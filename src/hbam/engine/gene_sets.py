"""Gene set management: literature, data-driven, and hybrid approaches."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from loguru import logger

from hbam.utils.logging import log_metric, log_step


# Default MSigDB hallmark sets for aging/muscle
DEFAULT_MATURATION_SETS = [
    "HALLMARK_MYOGENESIS",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_MTORC1_SIGNALING",
]

DEFAULT_DYSFUNCTION_SETS = [
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_P53_PATHWAY",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_APOPTOSIS",
]


def load_genage(path: Path | None = None) -> set[str]:
    """Load GenAge aging-associated gene set.

    Args:
        path: Path to GenAge CSV. If None, returns a curated subset.

    Returns:
        Set of human gene symbols.
    """
    with log_step("load_genage"):
        if path and Path(path).exists():
            df = pd.read_csv(path)
            gene_col = next((c for c in df.columns if "gene" in c.lower() or "symbol" in c.lower()), df.columns[0])
            genes = set(df[gene_col].dropna().astype(str).str.upper())
            log_metric("genage_genes", len(genes))
            return genes

        # Curated subset of well-known aging genes
        genes = {
            "TP53", "CDKN1A", "CDKN2A", "RB1", "TERT", "TERC",
            "SIRT1", "SIRT3", "SIRT6", "FOXO1", "FOXO3", "FOXO4",
            "MTOR", "RPTOR", "RICTOR", "AKT1", "IGF1", "IGF1R", "GHR",
            "LMNA", "PRKAA1", "PRKAA2", "ATM", "ATR", "BRCA1",
            "SOD1", "SOD2", "CAT", "GPX1", "NFE2L2",
            "IL6", "TNF", "IL1B", "CXCL8", "NFKB1", "RELA",
            "TGFB1", "BMP2", "WNT3A", "NOTCH1",
            "PPARGC1A", "TFAM", "NRF1", "POLG",
            "MSTN", "GDF11", "INHBA", "FST",
            "CTNNB1", "APC", "GSK3B",
            "HDAC1", "HDAC2", "EP300", "CREBBP",
            "MAP2K1", "MAPK1", "MAPK3", "JAK2", "STAT3",
        }
        log_metric("genage_genes_builtin", len(genes))
        return genes


def load_msigdb_hallmark(
    gene_set_names: list[str] | None = None,
    cache_dir: Path | None = None,
) -> dict[str, set[str]]:
    """Load MSigDB hallmark gene sets.

    Args:
        gene_set_names: Specific set names. If None, uses defaults.
        cache_dir: Directory for caching downloaded sets.

    Returns:
        Dict mapping set name -> set of gene symbols.
    """
    with log_step("load_msigdb_hallmark"):
        if gene_set_names is None:
            gene_set_names = DEFAULT_MATURATION_SETS + DEFAULT_DYSFUNCTION_SETS

        result = {}

        try:
            import gseapy as gp

            # Try to get hallmark gene sets
            hallmark = gp.get_library("MSigDB_Hallmark_2020")

            for name in gene_set_names:
                # Try exact match and variations
                for key in [name, name.replace("HALLMARK_", ""), name.replace("HALLMARK_", "").replace("_", " ").title()]:
                    if key in hallmark:
                        result[name] = set(hallmark[key])
                        break

                if name not in result:
                    logger.warning(f"Gene set '{name}' not found in MSigDB hallmark")

        except Exception as e:
            logger.warning(f"Failed to load MSigDB via gseapy: {e}. Using built-in fallback.")
            result = _builtin_hallmark_sets(gene_set_names)

        for name, genes in result.items():
            log_metric(f"msigdb_{name}", len(genes))

        return result


def _builtin_hallmark_sets(names: list[str]) -> dict[str, set[str]]:
    """Minimal built-in hallmark gene sets for offline use."""
    sets = {
        "HALLMARK_MYOGENESIS": {
            "MYOD1", "MYF5", "MYF6", "MYOG", "PAX7", "PAX3",
            "MEF2A", "MEF2C", "MEF2D", "DES", "TTN", "MYH1",
            "MYH2", "MYH3", "MYH7", "ACTA1", "ACTC1", "CKM",
            "MB", "TNNI1", "TNNI2", "TNNT1", "TNNT3", "TPM1",
        },
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION": {
            "NDUFA1", "NDUFA2", "NDUFB1", "NDUFS1", "NDUFS2",
            "SDHA", "SDHB", "UQCRC1", "UQCRC2", "COX4I1",
            "COX5A", "COX6A1", "ATP5F1A", "ATP5F1B", "ATP5MC1",
            "CS", "IDH1", "IDH2", "OGDH", "DLST", "FH",
        },
        "HALLMARK_MTORC1_SIGNALING": {
            "MTOR", "RPTOR", "RPS6KB1", "EIF4EBP1", "EIF4E",
            "RPS6", "ULK1", "ATG13", "DEPTOR", "MLST8",
            "AKT1", "TSC1", "TSC2", "RHEB", "PRAS40",
        },
        "HALLMARK_INFLAMMATORY_RESPONSE": {
            "IL6", "IL1B", "TNF", "CXCL8", "CCL2", "CCL5",
            "NFKB1", "RELA", "NFKBIA", "IKBKB", "TLR4",
            "MYD88", "IRAK1", "TRAF6", "CXCL1", "CXCL2",
            "IL1A", "IL18", "PTGS2", "NOS2", "ICAM1",
        },
        "HALLMARK_P53_PATHWAY": {
            "TP53", "MDM2", "CDKN1A", "BAX", "BBC3",
            "PMAIP1", "GADD45A", "GADD45B", "SESN1", "SESN2",
            "RRM2B", "DDB2", "XPC", "TIGAR", "SCO2",
            "PERP", "FAS", "TNFRSF10B", "PTEN", "TSC2",
        },
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB": {
            "NFKB1", "RELA", "NFKBIA", "TNFAIP3", "BIRC3",
            "TRAF1", "TRAF2", "IKBKB", "CHUK", "MAP3K7",
            "RIPK1", "TRADD", "FADD", "CASP8", "CFLAR",
            "BCL2L1", "XIAP", "CXCL1", "CCL20", "IL8",
        },
        "HALLMARK_APOPTOSIS": {
            "BAX", "BAK1", "BCL2", "BCL2L1", "MCL1",
            "BID", "BAD", "BIK", "BMF", "HRK",
            "CYCS", "APAF1", "CASP9", "CASP3", "CASP7",
            "DIABLO", "XIAP", "BIRC5", "FAS", "FASLG",
        },
    }

    return {n: sets.get(n, set()) for n in names if n in sets}


def load_custom_gene_set(path: Path) -> set[str]:
    """Load custom gene set from file.

    Supports: one gene per line (.txt) or GMT format (.gmt).

    Args:
        path: Path to gene set file.

    Returns:
        Set of gene symbols.
    """
    path = Path(path)

    if path.suffix == ".gmt":
        genes = set()
        with open(path) as f:
            for line in f:
                parts = line.strip().split("\t")
                genes.update(parts[2:])  # Skip name and description
        return genes
    else:
        with open(path) as f:
            return {line.strip().upper() for line in f if line.strip()}


def categorize_genes(
    hallmark_sets: dict[str, set[str]],
    genage_genes: set[str] | None = None,
    maturation_set_names: list[str] | None = None,
    dysfunction_set_names: list[str] | None = None,
) -> dict[str, set[str]]:
    """Categorize genes into maturation and dysfunction groups.

    Args:
        hallmark_sets: Dict of hallmark gene sets.
        genage_genes: GenAge aging genes (added to dysfunction).
        maturation_set_names: Which hallmark sets count as maturation.
        dysfunction_set_names: Which hallmark sets count as dysfunction.

    Returns:
        Dict with "maturation" and "dysfunction" gene sets.
    """
    if maturation_set_names is None:
        maturation_set_names = DEFAULT_MATURATION_SETS
    if dysfunction_set_names is None:
        dysfunction_set_names = DEFAULT_DYSFUNCTION_SETS

    maturation = set()
    for name in maturation_set_names:
        if name in hallmark_sets:
            maturation |= hallmark_sets[name]

    dysfunction = set()
    for name in dysfunction_set_names:
        if name in hallmark_sets:
            dysfunction |= hallmark_sets[name]

    if genage_genes:
        dysfunction |= genage_genes

    # Remove overlap (assign to dysfunction if in both)
    overlap = maturation & dysfunction
    if overlap:
        maturation -= overlap
        logger.info(f"Removed {len(overlap)} overlapping genes from maturation")

    log_metric("maturation_genes", len(maturation))
    log_metric("dysfunction_genes", len(dysfunction))

    return {"maturation": maturation, "dysfunction": dysfunction}


def refine_gene_sets(
    mudata: Any,
    literature_sets: dict[str, set[str]],
    method: str = "trajectory",
    factor_key: str = "pca_factors",
) -> dict[str, set[str]]:
    """Refine gene sets using data-driven signals.

    Args:
        mudata: MuData with latent factors.
        literature_sets: Initial literature-based gene sets.
        method: "trajectory" or "clustering".
        factor_key: Key in .obsm for latent factors.

    Returns:
        Refined gene sets with source tracking in metadata.
    """
    with log_step("refine_gene_sets", method=method):
        import mudata as md

        refined = {k: set(v) for k, v in literature_sets.items()}

        # Get first modality with gene expression data
        first_mod_name = list(mudata.mod.keys())[0]
        mod = mudata.mod[first_mod_name]

        if factor_key not in mudata.obsm and "mofa_factors" not in mudata.obsm:
            logger.warning("No latent factors found. Returning literature sets unchanged.")
            return refined

        factors = mudata.obsm.get(factor_key, mudata.obsm.get("mofa_factors"))

        if method == "trajectory":
            # Use first PC/factor as pseudo-trajectory
            trajectory = factors[:, 0]

            # Correlate each gene with trajectory
            from scipy.stats import spearmanr

            X = mod.X
            n_added_mat = 0
            n_added_dys = 0

            for j, gene in enumerate(mod.var_names):
                if gene in refined["maturation"] or gene in refined["dysfunction"]:
                    continue

                col = X[:, j]
                valid = ~np.isnan(col) if np.issubdtype(col.dtype, np.floating) else np.ones(len(col), dtype=bool)

                if valid.sum() < 5:
                    continue

                # Only correlate if we have matching sample count
                if valid.sum() == len(trajectory):
                    r, p = spearmanr(trajectory[valid], col[valid])
                elif valid.sum() < len(trajectory):
                    r, p = spearmanr(trajectory[:valid.sum()], col[valid])
                else:
                    continue

                if p < 0.01 and abs(r) > 0.5:
                    if r < 0:
                        refined["maturation"].add(gene)
                        n_added_mat += 1
                    else:
                        refined["dysfunction"].add(gene)
                        n_added_dys += 1

            log_metric("refined_added_maturation", n_added_mat)
            log_metric("refined_added_dysfunction", n_added_dys)

        elif method == "clustering":
            from sklearn.cluster import KMeans

            # Cluster genes by their loading patterns
            if "pca_loadings" in mod.varm:
                loadings = mod.varm["pca_loadings"]
            else:
                # Compute gene-gene correlations as features
                X = mod.X.copy()
                if np.isnan(X).any():
                    from sklearn.impute import SimpleImputer
                    X = SimpleImputer(strategy="mean").fit_transform(X)
                loadings = np.corrcoef(X.T)[:, :min(10, X.shape[0])]

            kmeans = KMeans(n_clusters=2, random_state=42, n_init=10)
            labels = kmeans.fit_predict(loadings)

            # Assign clusters based on overlap with literature sets
            for cluster_id in [0, 1]:
                cluster_genes = {mod.var_names[i] for i in range(len(labels)) if labels[i] == cluster_id}
                mat_overlap = len(cluster_genes & literature_sets.get("maturation", set()))
                dys_overlap = len(cluster_genes & literature_sets.get("dysfunction", set()))

                if mat_overlap > dys_overlap:
                    refined["maturation"] |= cluster_genes
                else:
                    refined["dysfunction"] |= cluster_genes

        # Remove overlap
        overlap = refined["maturation"] & refined["dysfunction"]
        if overlap:
            refined["maturation"] -= overlap

        log_metric("final_maturation_genes", len(refined["maturation"]))
        log_metric("final_dysfunction_genes", len(refined["dysfunction"]))

        return refined


def run_pathway_enrichment(
    gene_list: list[str],
    weights: pd.DataFrame | None = None,
    gene_sets_library: str = "MSigDB_Hallmark_2020",
    fdr_threshold: float = 0.05,
) -> pd.DataFrame:
    """Run pathway enrichment analysis using gseapy.

    Tries prerank (GSEA) if weights are provided, otherwise uses enrichr (ORA).

    Args:
        gene_list: List of gene symbols to test.
        weights: Optional DataFrame with gene and weight columns for ranked analysis.
        gene_sets_library: Gene set library name for gseapy.
        fdr_threshold: FDR threshold for significance.

    Returns:
        DataFrame with Term, P-value, Adjusted P-value, Genes columns.
    """
    with log_step("pathway_enrichment"):
        try:
            import gseapy as gp

            if weights is not None and len(weights) > 0:
                # Prerank GSEA using gene weights
                ranked = weights[["gene", "weight"]].copy()
                ranked = ranked.drop_duplicates(subset="gene")
                ranked = ranked.set_index("gene")["weight"].sort_values(ascending=False)

                try:
                    result = gp.prerank(
                        rnk=ranked,
                        gene_sets=gene_sets_library,
                        outdir=None,
                        seed=42,
                        verbose=False,
                        no_plot=True,
                    )
                    df = result.res2d
                    if len(df) > 0:
                        df = df.rename(columns={"FDR q-val": "Adjusted P-value", "NOM p-val": "P-value"})
                        df = df.sort_values("Adjusted P-value")
                        log_metric("enriched_pathways_fdr005", int((df["Adjusted P-value"] < fdr_threshold).sum()))
                        return df
                except Exception as e:
                    logger.warning(f"Prerank GSEA failed: {e}. Trying enrichr.")

            # Fallback: Over-Representation Analysis via enrichr
            result = gp.enrichr(
                gene_list=gene_list[:500],  # enrichr limit
                gene_sets=gene_sets_library,
                outdir=None,
                verbose=False,
                no_plot=True,
            )
            df = result.results
            if len(df) > 0:
                df = df.sort_values("Adjusted P-value")
                log_metric("enriched_pathways_fdr005", int((df["Adjusted P-value"] < fdr_threshold).sum()))
                return df

        except ImportError:
            logger.warning("gseapy not available. Skipping pathway enrichment.")
        except Exception as e:
            logger.warning(f"Pathway enrichment failed: {e}")

        # Return empty DataFrame with expected columns
        return pd.DataFrame(columns=["Term", "P-value", "Adjusted P-value", "Genes"])
