"""Latent embedding via MOFA+ or PCA."""

from __future__ import annotations

import mudata as md
import numpy as np
from loguru import logger

from hbam.utils.logging import log_metric, log_step


def embed_mofa(
    mudata: md.MuData,
    n_factors: int = 15,
    seed: int = 42,
) -> md.MuData:
    """Run MOFA+ multi-omics factor analysis.

    Args:
        mudata: Input MuData (aligned and scaled).
        n_factors: Number of latent factors.
        seed: Random seed.

    Returns:
        MuData with .obsm["mofa_factors"] and .varm["mofa_weights"] per modality.
    """
    with log_step("embed_mofa", n_factors=n_factors):
        try:
            from muon import tl as mu_tl

            mu_tl.mofa(
                mudata,
                n_factors=n_factors,
                seed=seed,
                convergence_mode="slow",
                gpu_mode=False,
            )

            # Extract factors
            if "X_mofa" in mudata.obsm:
                mudata.obsm["mofa_factors"] = mudata.obsm["X_mofa"]

            log_metric("mofa_factors", n_factors)
            logger.info("MOFA+ embedding completed successfully")

        except ImportError:
            logger.warning("muon/mofapy2 not available. Use embed_pca as fallback.")
            raise

    return mudata


def embed_pca(
    mudata: md.MuData,
    n_components: int = 15,
    seed: int = 42,
) -> md.MuData:
    """Run PCA on concatenated modalities.

    Args:
        mudata: Input MuData (aligned and scaled).
        n_components: Number of PCA components.
        seed: Random seed.

    Returns:
        MuData with .obsm["pca_factors"] and variance explained.
    """
    with log_step("embed_pca", n_components=n_components):
        from sklearn.decomposition import PCA
        from sklearn.impute import SimpleImputer

        # Concatenate modalities for joint PCA
        # We need to handle the fact that different modalities may have different samples
        # For now, use the first modality's approach
        matrices = []
        mod_names = []

        for name, mod in mudata.mod.items():
            X = mod.X.copy()
            if np.isnan(X).any():
                imputer = SimpleImputer(strategy="mean")
                X = imputer.fit_transform(X)
            matrices.append(X)
            mod_names.append(name)

        # If modalities have same samples (vertical integration), concatenate features
        # If different samples (horizontal), use first modality
        if all(m.shape[0] == matrices[0].shape[0] for m in matrices):
            X_concat = np.hstack(matrices)
            logger.info(f"Concatenated {len(matrices)} modalities: {X_concat.shape}")
        else:
            # Use largest modality
            largest_idx = np.argmax([m.shape[0] for m in matrices])
            X_concat = matrices[largest_idx]
            logger.warning(f"Modalities have different sample counts. Using '{mod_names[largest_idx]}' only.")

        n_components = min(n_components, *X_concat.shape)

        pca = PCA(n_components=n_components, random_state=seed)
        factors = pca.fit_transform(X_concat)

        # Store results
        # mudata.obsm requires shape (mudata.n_obs, ...) which spans ALL modalities.
        # When modalities have different obs (horizontal integration), factors only
        # cover the concatenated obs of same-obs modalities, so store in uns instead.
        mudata.uns["pca_variance_explained"] = pca.explained_variance_ratio_.tolist()
        mudata.uns["pca_cumulative_variance"] = np.cumsum(pca.explained_variance_ratio_).tolist()
        if factors.shape[0] == mudata.n_obs:
            mudata.obsm["pca_factors"] = factors
        else:
            # Store per-modality in uns; also push into each same-obs modality's obsm
            mudata.uns["pca_factors"] = factors.tolist()
            for name, mod in mudata.mod.items():
                if mod.n_obs == factors.shape[0]:
                    mod.obsm["pca_factors"] = factors
                    break
            # Expose via mudata.obsm using the first modality's obs padded with NaN
            padded = np.full((mudata.n_obs, factors.shape[1]), np.nan)
            padded[: factors.shape[0]] = factors
            mudata.obsm["pca_factors"] = padded

        # Store loadings per modality
        loadings = pca.components_.T  # (n_features, n_components)
        offset = 0
        for name, mod in mudata.mod.items():
            n_feat = mod.n_vars
            if offset + n_feat <= loadings.shape[0]:
                mod.varm["pca_loadings"] = loadings[offset:offset + n_feat]
                offset += n_feat

        total_var = sum(pca.explained_variance_ratio_)
        log_metric("pca_components", n_components)
        log_metric("pca_total_variance_explained", f"{total_var:.4f}")

        return mudata


def auto_select_embedding(mudata: md.MuData) -> str:
    """Auto-select embedding method based on data characteristics.

    Uses MOFA+ if: multiple modalities with overlapping samples AND mofapy2 available.
    Falls back to PCA otherwise.

    Args:
        mudata: Input MuData.

    Returns:
        "mofa" or "pca".
    """
    # Check if MOFA+ is available
    mofa_available = False
    try:
        import mofapy2
        mofa_available = True
    except ImportError:
        pass

    try:
        import muon
        mofa_available = True
    except ImportError:
        pass

    n_modalities = len(mudata.mod)

    if mofa_available and n_modalities > 1:
        logger.info("Auto-selected MOFA+ (multiple modalities + mofapy2 available)")
        return "mofa"
    else:
        reason = "single modality" if n_modalities <= 1 else "mofapy2 not available"
        logger.info(f"Auto-selected PCA ({reason})")
        return "pca"


def run_embedding(
    mudata: md.MuData,
    method: str = "auto",
    n_factors: int = 15,
    seed: int = 42,
) -> md.MuData:
    """Run embedding with auto-selection and fallback.

    Args:
        mudata: Input MuData.
        method: "auto", "mofa", or "pca".
        n_factors: Number of factors/components.
        seed: Random seed.

    Returns:
        MuData with embedding results.
    """
    with log_step("run_embedding", method=method):
        if method == "auto":
            method = auto_select_embedding(mudata)

        if method == "mofa":
            try:
                result = embed_mofa(mudata, n_factors=n_factors, seed=seed)
            except (ImportError, Exception) as e:
                logger.warning(f"MOFA+ failed: {e}. Falling back to PCA.")
                result = embed_pca(mudata, n_components=n_factors, seed=seed)
        else:
            result = embed_pca(mudata, n_components=n_factors, seed=seed)

        # Silhouette score validation
        _validate_embedding_quality(result)

        return result


def _validate_embedding_quality(mudata: md.MuData) -> None:
    """Validate embedding quality using silhouette score on biological groups."""
    factor_key = None
    for key in ["pca_factors", "mofa_factors", "X_mofa"]:
        if key in mudata.obsm:
            factor_key = key
            break

    if factor_key is None:
        return

    factors = mudata.obsm[factor_key]
    first_mod = list(mudata.mod.keys())[0]
    mod = mudata.mod[first_mod]

    # Find a categorical column for clustering validation
    label_col = None
    for col in ["condition", "batch", "region"]:
        if col in mod.obs.columns and mod.obs[col].nunique() >= 2:
            label_col = col
            break

    if label_col is None or factors.shape[0] != mod.n_obs:
        logger.info("Silhouette validation skipped: no suitable labels or shape mismatch")
        return

    try:
        from sklearn.metrics import silhouette_score

        labels = mod.obs[label_col].values
        # Only compute if we have enough samples per group
        from collections import Counter
        counts = Counter(labels)
        if min(counts.values()) < 2:
            return

        score = silhouette_score(factors, labels)
        mudata.uns["silhouette_score"] = float(score)
        mudata.uns["silhouette_label"] = label_col

        if score > 0.3:
            logger.info(f"Silhouette score: {score:.4f} (PASS, threshold=0.3, labels={label_col})")
        else:
            logger.warning(f"Silhouette score: {score:.4f} (LOW, threshold=0.3, labels={label_col}). "
                          "Embedding may not capture biological structure well.")
    except Exception as e:
        logger.warning(f"Silhouette validation failed: {e}")
