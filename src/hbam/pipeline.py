"""End-to-end HBAM pipeline orchestration."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import mudata as md
import numpy as np
import pandas as pd
from loguru import logger

from hbam.config import PipelineConfig
from hbam.utils.logging import log_step, setup_logging
from hbam.utils.validation import ValidationGate, validate_missingness


@dataclass
class PipelineResult:
    """Container for pipeline outputs."""
    mudata: md.MuData | None = None
    scores: pd.DataFrame | None = None
    weights: pd.DataFrame | None = None
    enrichment: pd.DataFrame | None = None
    figure_paths: list[Path] = field(default_factory=list)
    table_paths: list[Path] = field(default_factory=list)
    validation_report: dict[str, Any] = field(default_factory=dict)


def set_seeds(seed: int) -> None:
    """Set all random seeds for reproducibility."""
    import random
    random.seed(seed)
    np.random.seed(seed)
    logger.info(f"Random seeds set to {seed}")


def run_pipeline(config: PipelineConfig) -> PipelineResult:
    """Execute the full HBAM pipeline.

    Args:
        config: Pipeline configuration.

    Returns:
        PipelineResult with all outputs.
    """
    with log_step("full_pipeline"):
        setup_logging(log_level=config.log_level, log_file=config.log_file)
        set_seeds(config.seed)

        result = PipelineResult()
        output_dir = Path(config.data.output_dir)
        intermediate_dir = Path(config.data.intermediate_dir)
        intermediate_dir.mkdir(parents=True, exist_ok=True)

        # ===== PHASE 1: Data Loading =====
        modalities = _load_all_modalities(config)

        if not modalities:
            logger.warning("No modalities loaded. Generating synthetic data for demo.")
            modalities = _generate_demo_data(config.seed)

        # ===== PHASE 2: QC + Normalization + Imputation =====
        modalities = _process_modalities(modalities, config)

        # ===== PHASE 3: Modality-specific Analysis =====
        modalities = _analyze_modalities(modalities, config)

        # ===== PHASE 4: Integration =====
        mudata = _integrate_modalities(modalities, config)
        result.mudata = mudata

        # Save intermediate
        try:
            mudata.write(str(intermediate_dir / "integrated.h5mu"))
        except Exception as e:
            logger.warning(f"Could not save intermediate MuData: {e}")

        # ===== PHASE 5: HBAM Engine =====
        weights, gene_sets, enrichment = _run_hbam_engine(mudata, config)
        result.weights = weights
        result.enrichment = enrichment

        # ===== PHASE 6: Translation =====
        scores = _run_translation(mudata, weights, config)
        result.scores = scores

        # ===== PHASE 7: Output =====
        figure_paths, table_paths = _generate_outputs(
            mudata, scores, weights, enrichment, config
        )
        result.figure_paths = figure_paths
        result.table_paths = table_paths

        logger.info(f"Pipeline complete: {len(figure_paths)} figures, {len(table_paths)} tables")

        return result


def _is_diann_matrix(path: Path) -> bool:
    """Check if a file is a DIA-NN wide pg_matrix (not long report)."""
    import csv
    with open(path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
    return "Protein.Group" in header and "Run" not in header


def _load_all_modalities(config: PipelineConfig) -> dict[str, Any]:
    """Load all configured modalities."""
    from hbam.data.loaders import load_maxquant, load_diann, load_diann_matrix, load_matrix, load_stereo_gem

    modalities = {}

    for mod_input in config.data.modalities:
        with log_step("load_modality", modality_name=mod_input.name, format=mod_input.format):
            try:
                path = Path(mod_input.path)
                if not path.exists():
                    logger.warning(f"File not found: {path}. Skipping {mod_input.name}.")
                    continue

                loader_map = {
                    "maxquant": load_maxquant,
                    "diann": load_diann,
                    "diann_matrix": load_diann_matrix,
                    "csv": load_matrix,
                    "tsv": load_matrix,
                    "gem": load_stereo_gem,
                }

                # Auto-detect DIA-NN wide matrix format
                loader = loader_map.get(mod_input.format)
                if mod_input.format == "diann" and _is_diann_matrix(path):
                    loader = load_diann_matrix
                if loader is None:
                    logger.warning(f"Unknown format '{mod_input.format}'. Trying generic loader.")
                    loader = load_matrix

                adata = loader(path)

                # Store metadata
                adata.uns["species"] = mod_input.species
                adata.uns["modality_type"] = mod_input.modality_type

                # Inject sample metadata from companion file if available
                meta_path = Path(str(path).rsplit(".", 1)[0] + "_metadata.tsv")
                if not meta_path.exists():
                    # Try data/processed/ directory
                    meta_path = Path("data/processed") / f"{mod_input.name.replace(' ', '_')}_metadata.tsv"

                if meta_path.exists():
                    meta = pd.read_csv(meta_path, sep="\t")
                    # Match metadata to samples by sample_id
                    if "sample_id" in meta.columns:
                        meta = meta.set_index("sample_id")
                        shared = adata.obs_names.intersection(meta.index)
                        if len(shared) > 0:
                            for col in meta.columns:
                                adata.obs[col] = meta.loc[adata.obs_names, col].values if all(s in meta.index for s in adata.obs_names) else None
                            logger.info(f"Injected metadata from {meta_path}: {list(meta.columns)}")
                    else:
                        # Try row-order matching if same length
                        if len(meta) == adata.n_obs:
                            for col in meta.columns:
                                adata.obs[col] = meta[col].values
                            logger.info(f"Injected metadata (row-order) from {meta_path}")

                modalities[mod_input.name] = adata

            except Exception as e:
                logger.error(f"Failed to load {mod_input.name}: {e}")

    return modalities


def _generate_demo_data(seed: int = 42) -> dict[str, Any]:
    """Generate synthetic demo data when no real data is configured."""
    import anndata as ad

    rng = np.random.default_rng(seed)

    # Proteomics
    n_samples, n_proteins = 20, 200
    X_prot = rng.lognormal(mean=20, sigma=2, size=(n_samples, n_proteins))
    conditions = ["young"] * 10 + ["old"] * 10
    # Inject aging effect
    X_prot[10:, :30] *= rng.lognormal(0.5, 0.3, size=(10, 30))
    # Add missingness
    mask = rng.random(X_prot.shape) < 0.08
    X_prot[mask] = np.nan

    gene_names = [f"GENE_{i:04d}" for i in range(n_proteins)]
    real_genes = ["TP53", "MTOR", "CDKN1A", "IL6", "TNF", "MYOD1", "PAX7", "MYH1",
                  "SIRT1", "FOXO3", "IGF1", "MSTN", "DES", "TTN", "ACTA1",
                  "COL1A1", "NFKB1", "TGFB1", "SOD2", "CAT"]
    for i, g in enumerate(real_genes[:min(len(real_genes), n_proteins)]):
        gene_names[i] = g

    obs_prot = pd.DataFrame({
        "sample_id": [f"S{i:03d}" for i in range(n_samples)],
        "condition": conditions,
        "timepoint": [i // 5 + 1 for i in range(n_samples)],
    })
    obs_prot.index = obs_prot["sample_id"]
    var_prot = pd.DataFrame({"gene_names": gene_names}, index=gene_names)

    prot = ad.AnnData(X=X_prot, obs=obs_prot, var=var_prot, layers={"raw": X_prot.copy()})
    prot.uns["species"] = "human"
    prot.uns["modality_type"] = "temporal"

    # Spatial
    n_bins, n_genes_sp = 100, 150
    X_sp = rng.poisson(5, size=(n_bins, n_genes_sp)).astype(np.float32)
    zero_mask = rng.random(X_sp.shape) < 0.5
    X_sp[zero_mask] = 0
    coords = np.array([(i % 10, i // 10) for i in range(n_bins)], dtype=np.float32)

    sp_genes = gene_names[:n_genes_sp]  # Share some genes
    obs_sp = pd.DataFrame({"bin_id": [f"B{i:04d}" for i in range(n_bins)],
                           "region": [f"R{i%4}" for i in range(n_bins)]})
    obs_sp.index = obs_sp["bin_id"]
    var_sp = pd.DataFrame({"gene_names": sp_genes}, index=sp_genes)

    spatial = ad.AnnData(X=X_sp, obs=obs_sp, var=var_sp,
                         obsm={"spatial": coords}, layers={"raw": X_sp.copy()})
    spatial.uns["species"] = "human"
    spatial.uns["modality_type"] = "spatial"

    return {"proteomics": prot, "spatial": spatial}


def _process_modalities(modalities: dict, config: PipelineConfig) -> dict:
    """QC, normalize, and impute each modality."""
    from hbam.data.qc import filter_missingness, filter_variance
    from hbam.data.normalize import normalize_proteomics, normalize_spatial, normalize_generic
    from hbam.data.impute import impute_auto

    processed = {}

    for name, adata in modalities.items():
        with log_step("process_modality", modality_name=name):
            # QC
            adata = filter_missingness(adata, threshold=config.qc.min_presence, axis="var")
            adata = filter_variance(adata, percentile=config.qc.variance_percentile)

            # Normalize based on modality type
            mod_type = adata.uns.get("modality_type", "generic")
            if mod_type in ("temporal", "functional"):
                adata = normalize_proteomics(
                    adata,
                    method=config.normalization.proteomics_method,
                    pseudocount=config.normalization.pseudocount,
                )
            elif mod_type == "spatial":
                adata = normalize_spatial(adata, method=config.normalization.spatial_method)
            else:
                adata = normalize_generic(adata)

            # Impute
            if config.imputation.method != "none":
                adata = impute_auto(
                    adata,
                    seed=config.seed,
                    mnar_method=config.imputation.mnar_method,
                    mnar_quantile=config.imputation.mnar_quantile,
                    mcar_method=config.imputation.mcar_method,
                    mcar_k=config.imputation.knn_k,
                )

            processed[name] = adata

    return processed


def _analyze_modalities(modalities: dict, config: PipelineConfig) -> dict:
    """Run modality-specific analysis."""
    from hbam.modality.temporal import run_temporal_analysis
    from hbam.modality.spatial import run_spatial_analysis
    from hbam.modality.functional import run_functional_analysis

    analyzed = {}

    for name, adata in modalities.items():
        mod_type = adata.uns.get("modality_type", "generic")

        with log_step("analyze_modality", modality_name=name, type=mod_type):
            if mod_type == "temporal":
                adata = run_temporal_analysis(adata)
            elif mod_type == "spatial":
                adata = run_spatial_analysis(adata)
            elif mod_type == "functional":
                adata = run_functional_analysis(adata)
            else:
                logger.info(f"No specific analysis for modality type '{mod_type}'")

            analyzed[name] = adata

    return analyzed


def _integrate_modalities(modalities: dict, config: PipelineConfig) -> md.MuData:
    """Integrate modalities: orthologs, alignment, scaling, embedding."""
    from hbam.integration.orthologs import load_ortholog_table, harmonize_gene_names
    from hbam.integration.align import align_genes, scale_features
    from hbam.integration.embed import run_embedding

    with log_step("integration"):
        # Ortholog mapping
        ortho_table = load_ortholog_table(source="builtin")

        harmonized = {}
        for name, adata in modalities.items():
            species = adata.uns.get("species", "human")
            harmonized[name] = harmonize_gene_names(adata, species, ortho_table)

        # Create MuData
        mudata = md.MuData(harmonized)

        # Gene alignment
        try:
            mudata = align_genes(mudata, min_overlap=config.integration.min_gene_overlap)
        except Exception as e:
            logger.warning(f"Gene alignment failed: {e}. Using unaligned data.")

        # Feature scaling
        mudata = scale_features(mudata, method=config.integration.scaling_method)

        # Embedding
        mudata = run_embedding(
            mudata,
            method=config.integration.embedding_method,
            n_factors=config.integration.n_factors,
            seed=config.seed,
        )

        return mudata


def _run_hbam_engine(mudata: md.MuData, config: PipelineConfig) -> tuple:
    """Compute gene sets, weights, and HBAM scores."""
    from hbam.engine.gene_sets import (
        load_genage, load_msigdb_hallmark, categorize_genes, refine_gene_sets,
        run_pathway_enrichment,
    )
    from hbam.engine.weights import compute_weights_pca, normalize_weights
    from hbam.engine.score import compute_hbam_score, validate_score_direction, bootstrap_validate

    with log_step("hbam_engine"):
        # Gene sets
        genage = load_genage()
        hallmark = load_msigdb_hallmark()
        gene_sets = categorize_genes(hallmark, genage)

        # Data-driven refinement
        if config.engine.gene_set_source in ("hybrid", "data_driven"):
            gene_sets = refine_gene_sets(mudata, gene_sets, method=config.engine.refinement_method)

        # Weights
        weights = compute_weights_pca(mudata, gene_sets)

        if len(weights) == 0:
            logger.warning("No weights computed. Using uniform weights for available genes.")
            first_mod = list(mudata.mod.keys())[0]
            rows = []
            for gene in mudata.mod[first_mod].var_names[:50]:
                gene_upper = gene.upper()
                if gene_upper in gene_sets.get("dysfunction", set()):
                    rows.append({"gene": gene, "weight": 1.0, "category": "dysfunction", "modality": first_mod})
                elif gene_upper in gene_sets.get("maturation", set()):
                    rows.append({"gene": gene, "weight": 1.0, "category": "maturation", "modality": first_mod})
            weights = pd.DataFrame(rows) if rows else pd.DataFrame(columns=["gene", "weight", "category", "modality"])

        weights = normalize_weights(weights)

        # Scoring
        if len(weights) > 0:
            scores = compute_hbam_score(mudata, weights)

            # Direction validation
            first_mod = list(mudata.mod.keys())[0]
            validate_score_direction(scores, mudata.mod[first_mod].obs)

            # Bootstrap validation
            bootstrap_result = bootstrap_validate(
                mudata, gene_sets, compute_weights_pca,
                n_iter=config.engine.bootstrap_iterations,
                fraction=config.engine.bootstrap_fraction,
                threshold=config.engine.bootstrap_min_spearman,
                seed=config.seed,
            )
            logger.info(f"Bootstrap: mean_r={bootstrap_result.mean_correlation:.4f}, passed={bootstrap_result.passed}")

        # Pathway enrichment
        enrichment = None
        if len(weights) > 0:
            gene_list = weights["gene"].tolist()
            enrichment = run_pathway_enrichment(gene_list, weights)
            if enrichment is not None and len(enrichment) > 0:
                logger.info(f"Pathway enrichment: {len(enrichment)} terms found")

        return weights, gene_sets, enrichment


def _run_translation(mudata: md.MuData, weights: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    """Score samples and compute clinical metrics."""
    from hbam.translate.sample_score import score_samples, compare_conditions
    from hbam.translate.biomarkers import select_biomarkers
    from hbam.translate.clinical import compute_effect_sizes, discrimination_analysis, generate_clinical_summary

    with log_step("translation"):
        scores = score_samples(mudata, weights)

        # Statistical comparison
        if "condition" in scores.columns:
            comparison = compare_conditions(scores)
            logger.info(f"Condition comparison: p={comparison.get('p_value', 'N/A')}")

            effect_sizes = compute_effect_sizes(scores)

            try:
                disc = discrimination_analysis(scores)
                summary = generate_clinical_summary(scores, effect_sizes, disc)
                logger.info(f"Clinical summary generated. AUC={disc.get('auc', 'N/A')}")
            except Exception as e:
                logger.warning(f"Discrimination analysis failed: {e}")

        # Biomarker panel
        panel = select_biomarkers(weights, n=config.translation.top_biomarkers)
        logger.info(f"Biomarker panel: {len(panel)} genes selected")

        return scores


def _generate_outputs(
    mudata: md.MuData,
    scores: pd.DataFrame,
    weights: pd.DataFrame,
    enrichment: pd.DataFrame | None,
    config: PipelineConfig,
) -> tuple[list[Path], list[Path]]:
    """Generate all figures and tables."""
    from hbam.output.figures import (
        fig_pipeline_overview, fig_qc_summary, fig_latent_space,
        fig_weight_heatmap, fig_hbam_distribution, fig_pathway_enrichment,
        fig_volcano, fig_correlation_matrix, fig_biomarker_panel,
        fig_integration_diagnostic,
    )
    from hbam.output.tables import generate_summary_table, generate_gene_table, save_table
    from hbam.translate.biomarkers import select_biomarkers

    output_dir = Path(config.data.output_dir)
    fig_dir = output_dir / "figures"
    tab_dir = output_dir / "tables"

    figure_paths = []
    table_paths = []

    with log_step("generate_outputs"):
        # Figures
        figure_fns = [
            ("pipeline_overview", lambda: fig_pipeline_overview(fig_dir)),
            ("qc_summary", lambda: fig_qc_summary(mudata, fig_dir)),
            ("latent_space", lambda: fig_latent_space(mudata, fig_dir)),
            ("weight_heatmap", lambda: fig_weight_heatmap(weights, fig_dir, top_n=config.output.top_n_heatmap)),
            ("hbam_distribution", lambda: fig_hbam_distribution(scores, fig_dir)),
            ("pathway_enrichment", lambda: fig_pathway_enrichment(enrichment, fig_dir)),
            ("volcano", lambda: fig_volcano(None, fig_dir)),  # DE results from modality analysis
            ("correlation_matrix", lambda: fig_correlation_matrix(scores, fig_dir)),
            ("biomarker_panel", lambda: fig_biomarker_panel(select_biomarkers(weights, n=config.translation.top_biomarkers), mudata, fig_dir)),
            ("integration_diagnostic", lambda: fig_integration_diagnostic(mudata, fig_dir)),
        ]

        for name, fn in figure_fns:
            try:
                paths = fn()
                figure_paths.extend(paths)
            except Exception as e:
                logger.warning(f"Figure '{name}' failed: {e}")

        # Tables
        try:
            summary = generate_summary_table(mudata, scores, weights)
            table_paths.extend(save_table(summary, "summary", tab_dir))

            gene_table = generate_gene_table(weights)
            table_paths.extend(save_table(gene_table, "gene_weights", tab_dir))
        except Exception as e:
            logger.warning(f"Table generation failed: {e}")

    return figure_paths, table_paths
