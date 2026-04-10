# HBAM Multi-Omics Pipeline - Implementation Plan

**Created:** 2026-04-09
**Status:** PENDING USER APPROVAL
**Complexity:** HIGH (33 steps, ~60 files, greenfield project)

---

## Context

Build a reproducible Python pipeline for computing a High Burden Aging Muscle (HBAM) composite index from multi-omics data (mouse tissue proteomics + human iPSC-derived organoids + Stereo-seq spatial transcriptomics). The pipeline ingests heterogeneous data formats, performs modality-specific analysis, integrates across species via ortholog mapping, computes directional aging/dysfunction scores, and produces publication-ready figures.

## Work Objectives

1. Establish a modular, testable Python project with seed-controlled reproducibility
2. Build multi-format data ingestion with validation gates
3. Implement three parallelizable modality analysis sub-modules
4. Integrate modalities via cross-species ortholog mapping and latent embedding (MOFA+/PCA)
5. Compute the HBAM composite index with bootstrap validation
6. Generate a fixed catalog of 10 publication-ready figures
7. Assemble an end-to-end CLI pipeline

## Guardrails

### Must Have
- MuData as core data structure (scverse ecosystem)
- Seed control on every stochastic operation
- Validation gates between all modules (reject on >30% missing, gene mismatch, dimension inconsistency)
- Structured logging (loguru) at every processing step
- Each module independently importable and testable
- Colorblind-safe palette, vector PDF output for all figures
- Configuration-driven (YAML) with sensible defaults

### Must NOT Have
- Hard-coded file paths or magic numbers
- Silent data loss (every filter/drop must be logged with counts)
- Imputation of biological zeros in spatial data
- Batch correction when design matrix is fully confounded with biology
- Any dependency on interactive/GUI tools for core pipeline

---

## Phase 1: Foundation (Steps 1-5)

### Step 1: Project Scaffolding

**What:** Create the full directory tree, pyproject.toml, .gitignore, and empty `__init__.py` files for all packages.

**Files to create:**
- `pyproject.toml`
- `.gitignore`
- `README.md`
- `src/hbam/__init__.py` (with `__version__`)
- All sub-package `__init__.py` files (data/, modality/, integration/, engine/, translate/, output/, utils/)
- `tests/conftest.py`
- `config/default.yaml` (skeleton)
- Directory stubs: `data/raw/`, `data/processed/`, `data/external/`, `results/figures/`, `results/tables/`, `results/intermediate/`, `notebooks/`

**Key patterns:**
- `pyproject.toml` with `[project]` metadata, `[project.scripts]` for CLI entry point (`hbam = "hbam.cli:app"`)
- Dependencies: mudata, anndata, scanpy, mofapy2, pandas, numpy, scipy, scikit-learn, matplotlib, seaborn, gseapy, pybiomart, typer, pyyaml, loguru, pytest, pytest-cov
- Optional deps group `[project.optional-dependencies]` for dev tools (ruff, mypy)

**Acceptance criteria:**
- `pip install -e .` succeeds
- `python -c "import hbam; print(hbam.__version__)"` prints version string
- All directories exist

**Dependencies:** None (first step)

---

### Step 2: Configuration System

**What:** Build YAML-based configuration loading with validation and defaults.

**Files:**
- `src/hbam/config.py`
- `config/default.yaml`

**Key patterns:**
- Pydantic `BaseModel` (or dataclass) for typed config schema with nested sections: `DataConfig`, `QCConfig`, `NormConfig`, `ImputeConfig`, `IntegrationConfig`, `EngineConfig`, `TranslationConfig`, `OutputConfig`
- `load_config(path: Path) -> PipelineConfig` function that loads YAML, merges with defaults, validates
- All thresholds as config fields: `missingness_threshold: float = 0.3`, `variance_percentile: float = 10`, `min_genes_retained: int = 1000`, `bootstrap_iterations: int = 5`, `bootstrap_fraction: float = 0.8`, `bootstrap_min_spearman: float = 0.8`, `top_biomarkers: int = 20`
- `seed: int = 42` at top level

**Acceptance criteria:**
- Loading default.yaml returns fully populated PipelineConfig object
- Missing optional fields get defaults
- Invalid values (e.g., negative threshold) raise ValidationError
- Config is serializable back to YAML

**Dependencies:** Step 1

---

### Step 3: Logging Framework

**What:** Structured logging with loguru, configured from pipeline config.

**Files:**
- `src/hbam/utils/logging.py`

**Key patterns:**
- `setup_logging(config: PipelineConfig) -> None` that configures loguru sinks (file + stderr)
- Log file at `results/logs/pipeline_{timestamp}.log`
- Structured format: `{time} | {level} | {module}:{function}:{line} | {message}`
- Context manager `log_step(name: str)` that logs entry/exit with elapsed time
- Helper `log_filter(action: str, before: int, after: int, reason: str)` for data filtering events

**Acceptance criteria:**
- Calling `setup_logging()` creates log file
- `log_step` context manager records timing
- `log_filter` produces structured output like `FILTER | action=variance_filter | before=5000 | after=4500 | removed=500 | reason=below 10th percentile`

**Dependencies:** Step 2

---

### Step 4: Validation Framework

**What:** Gate validators that check data integrity between pipeline stages.

**Files:**
- `src/hbam/utils/validation.py`

**Key patterns:**
- `validate_missingness(adata: AnnData, threshold: float) -> ValidationResult` -- reject if fraction missing > threshold
- `validate_dimensions(mudata: MuData, expected: dict) -> ValidationResult` -- check n_obs, n_vars match expectations
- `validate_gene_overlap(genes_a: set, genes_b: set, min_overlap: int) -> ValidationResult` -- error if intersection < min_overlap
- `validate_normalization(adata: AnnData, fold_range: float = 1.5) -> ValidationResult` -- check median log2-intensity per sample within range
- `ValidationResult` dataclass with `passed: bool`, `message: str`, `metrics: dict`
- `ValidationGate` class that collects results and raises `PipelineValidationError` on failure

**Acceptance criteria:**
- Each validator returns correct pass/fail on synthetic data
- `ValidationGate` aggregates multiple checks
- Failed gate raises informative error with all failure details

**Dependencies:** Step 1

---

### Step 5: Test Fixtures and Conftest

**What:** Create synthetic test data fixtures for all pipeline stages.

**Files:**
- `tests/conftest.py`
- `tests/test_data/` (small fixture files)

**Key patterns:**
- `pytest` fixtures using `@pytest.fixture(scope="session")`
- `make_proteomics_adata(n_samples=10, n_proteins=500, seed=42)` -- returns AnnData with realistic structure
- `make_spatial_adata(n_bins=200, n_genes=300, seed=42)` -- returns AnnData with spatial coordinates in `.obsm["spatial"]`
- `make_mudata(modalities: dict)` -- wraps multiple AnnData into MuData
- `make_maxquant_file(tmp_path)` -- writes a minimal proteinGroups.txt
- `make_diann_file(tmp_path)` -- writes a minimal DIA-NN report.tsv
- `make_config()` -- returns PipelineConfig with test-appropriate settings
- Seed all random generators via `np.random.default_rng(seed)`

**Acceptance criteria:**
- All fixtures produce valid AnnData/MuData objects
- Fixture files parse without error
- Fixtures are deterministic (same seed = same data)

**Dependencies:** Steps 1-4

---

## Phase 2: Data Layer (Steps 6-10)

### Step 6: MaxQuant Proteomics Loader

**What:** Parse MaxQuant proteinGroups.txt into AnnData.

**Files:**
- `src/hbam/data/loaders.py`

**Key patterns:**
- `load_maxquant(path: Path, intensity_prefix: str = "LFQ intensity") -> AnnData`
- Filter rows: remove `Reverse`, `Potential contaminant`, `Only identified by site` entries (log counts)
- Extract intensity columns by prefix, transpose to samples x proteins
- Store protein IDs in `.var["protein_ids"]`, gene names in `.var["gene_names"]`
- Store raw values in `.layers["raw"]`, working matrix in `.X`
- Log: number of proteins loaded, samples detected, contaminants removed

**Acceptance criteria:**
- Parses fixture proteinGroups.txt correctly
- Contaminant/reverse rows excluded
- Dimensions match expected (samples x proteins)
- Gene names extracted and stored in `.var`

**Dependencies:** Step 5

**Tests:** `tests/test_loaders.py::test_load_maxquant`

---

### Step 7: DIA-NN and Generic CSV/TSV Loaders

**What:** Parse DIA-NN report.tsv (long format) and generic matrix files.

**Files:**
- `src/hbam/data/loaders.py` (extend)

**Key patterns:**
- `load_diann(path: Path, quantity_col: str = "PG.MaxLFQ") -> AnnData`
  - Pivot long-format: rows=Run, cols=ProteinGroup, values=quantity_col
  - Handle multiple charge states by taking max per protein group
- `load_matrix(path: Path, sep: str = "auto", gene_col: str = "auto") -> AnnData`
  - Auto-detect separator (comma vs tab)
  - Auto-detect gene name column (first non-numeric column)
  - Transpose if genes are rows (heuristic: more rows than columns = genes-as-rows)
- `load_stereo_gem(path: Path) -> AnnData`
  - Parse Stereo-seq GEM format (tab-separated: geneID, x, y, MIDcount)
  - Bin into spatial bins if needed (configurable bin size)
  - Store coordinates in `.obsm["spatial"]`

**Acceptance criteria:**
- DIA-NN long format pivots correctly to wide matrix
- CSV auto-detection works for comma and tab
- GEM parser creates AnnData with spatial coordinates
- All loaders produce consistent AnnData structure

**Dependencies:** Step 6

**Tests:** `tests/test_loaders.py::test_load_diann`, `test_load_matrix`, `test_load_gem`

---

### Step 8: Quality Control Module

**What:** Missingness filtering and variance filtering.

**Files:**
- `src/hbam/data/qc.py`

**Key patterns:**
- `filter_missingness(adata: AnnData, threshold: float = 0.7, axis: str = "var") -> AnnData`
  - Remove features present in fewer than `threshold` fraction of samples
  - Log: features before, after, removed count
- `filter_variance(adata: AnnData, percentile: float = 10) -> AnnData`
  - Remove features below the given percentile of variance (computed on non-NaN values)
  - Log: variance cutoff value, features removed
- `qc_report(adata: AnnData) -> dict`
  - Return summary: n_samples, n_features, missingness_fraction, median_values, variance_distribution
- All functions return new AnnData (no in-place mutation)

**Acceptance criteria:**
- 70% presence filter removes features missing in >30% of samples
- Variance filter removes bottom 10th percentile
- QC report returns correct statistics
- No silent data loss (all removals logged)
- Validation gate passes after QC on clean data

**Dependencies:** Steps 6-7

**Tests:** `tests/test_qc.py`

---

### Step 9: Normalization Module

**What:** Modality-appropriate normalization.

**Files:**
- `src/hbam/data/normalize.py`

**Key patterns:**
- `normalize_proteomics(adata: AnnData, method: str = "median", pseudocount: float = 1.0) -> AnnData`
  - Median centering: subtract per-sample median (in linear space), THEN log2(x + pseudocount)
  - Store normalized in `.X`, raw in `.layers["raw"]`
- `normalize_spatial(adata: AnnData, method: str = "sctransform") -> AnnData`
  - SCTransform via scanpy/scvi or scran normalization
  - Do NOT add pseudocount (SCTransform handles log internally)
  - Fallback to `scanpy.pp.normalize_total` + `scanpy.pp.log1p` if SCTransform unavailable
- `normalize_generic(adata: AnnData, method: str = "log2") -> AnnData`
  - log2(x + pseudocount) transformation
- After normalization, run `validate_normalization()` gate

**Acceptance criteria:**
- After median centering, per-sample medians within 1.5-fold range
- Log2 transform applied correctly (verify by exponentiating back)
- Raw values preserved in `.layers["raw"]`
- Validation gate passes post-normalization

**Dependencies:** Steps 4, 8

**Tests:** `tests/test_qc.py::test_normalize_proteomics`, `test_normalize_spatial`

---

### Step 10: Missing Value Imputation

**What:** MNAR vs MCAR detection and appropriate imputation.

**Files:**
- `src/hbam/data/impute.py`

**Key patterns:**
- `detect_missingness_type(adata: AnnData) -> str`
  - Heuristic: if missing values correlate with low intensity (left-censored pattern), classify as MNAR
  - Use Little's MCAR test approximation or intensity-dependence test
  - Return "MNAR", "MCAR", or "MIXED"
- `impute_mnar(adata: AnnData, method: str = "minprob", q: float = 0.01) -> AnnData`
  - MinProb: draw from N(q-th quantile of observed, sd=observed_sd * 0.3)
  - Seed-controlled via config
- `impute_mcar(adata: AnnData, method: str = "knn", k: int = 5) -> AnnData`
  - KNN imputation using scikit-learn KNNImputer
- `impute_auto(adata: AnnData, config: ImputeConfig) -> AnnData`
  - Detect type, route to appropriate method, log decision
- Store imputation mask in `.layers["imputed_mask"]` (boolean: True where imputed)

**Acceptance criteria:**
- MNAR detection correctly identifies left-censored patterns in synthetic data
- MinProb imputed values are in the low-intensity range
- KNN imputed values are reasonable (within observed range)
- Imputation mask correctly tracks which values were imputed
- No NaN values remain after imputation

**Dependencies:** Step 9

**Tests:** `tests/test_qc.py::test_impute_mnar`, `test_impute_mcar`, `test_impute_auto`

---

## Phase 3: Modality Analysis (Steps 11-13)

### Step 11: Temporal Agent (Liver Proteomics)

**What:** Time-series analysis and differential expression across timepoints for liver proteomics.

**Files:**
- `src/hbam/modality/temporal.py`

**Key patterns:**
- `run_temporal_analysis(adata: AnnData, time_col: str, config: PipelineConfig) -> AnnData`
- Differential expression: `scipy.stats.kruskal` across timepoints (non-parametric, small sample safe)
- Multiple testing correction: `statsmodels.stats.multitest.multipletests` with BH-FDR
- Store results in `.var`: `temporal_pvalue`, `temporal_fdr`, `temporal_fold_change`, `temporal_trend` (up/down/flat)
- Trend detection: Spearman correlation of protein abundance vs time
- Log: number of significant proteins at FDR < 0.05

**Acceptance criteria:**
- Differential expression runs on fixture data without error
- FDR correction applied correctly (FDR values >= raw p-values)
- Trend direction matches synthetic data pattern (injected trend should be detected)
- Results stored in `.var` columns

**Dependencies:** Steps 6-10

**Tests:** `tests/test_modality.py::test_temporal_analysis`

---

### Step 12: Spatial Agent (Stereo-seq)

**What:** Spatial QC, filtering, and pseudo-bulk or bin-level analysis for Stereo-seq data.

**Files:**
- `src/hbam/modality/spatial.py`

**Key patterns:**
- `run_spatial_analysis(adata: AnnData, config: PipelineConfig) -> AnnData`
- Spatial QC: filter bins by total counts (min_counts), genes detected (min_genes), mitochondrial fraction
- Zero-inflation handling: use NB-based methods; do NOT impute biological zeros
- Pseudo-bulk aggregation: `aggregate_pseudobulk(adata, group_col)` -- sum counts per group, return new AnnData
- Spatial variable gene detection: Moran's I or SpatialDE (if available), fallback to HVG via scanpy
- Store: `.var["spatially_variable"]`, `.var["morans_i"]`, `.obsm["spatial"]` preserved

**Acceptance criteria:**
- Spatial QC filters low-quality bins (logged)
- Pseudo-bulk produces correct aggregated counts (sum matches)
- Spatial coordinates preserved through all operations
- No imputation of zero-count entries
- Spatially variable genes identified

**Dependencies:** Steps 7-10

**Tests:** `tests/test_modality.py::test_spatial_analysis`

---

### Step 13: Functional Agent (Muscle Proteomics)

**What:** Functional annotation and protein group resolution for muscle proteomics.

**Files:**
- `src/hbam/modality/functional.py`

**Key patterns:**
- `run_functional_analysis(adata: AnnData, config: PipelineConfig) -> AnnData`
- Protein group resolution: when multiple gene names per protein group, select first (canonical) or split into separate entries (configurable)
- Functional annotation: map gene names to GO terms using gseapy or preloaded annotation table
- Differential expression: Wilcoxon rank-sum between conditions (young vs old, or treatment vs control)
- Store: `.var["go_terms"]`, `.var["functional_category"]`, `.var["de_pvalue"]`, `.var["de_fdr"]`, `.var["log2fc"]`

**Acceptance criteria:**
- Protein group resolution handles semicolon-separated gene names
- GO annotation maps genes correctly (at least partial coverage)
- Differential expression p-values are valid (0-1 range, FDR >= p)
- Results stored in `.var`

**Dependencies:** Steps 6-10

**Tests:** `tests/test_modality.py::test_functional_analysis`

---

## Phase 4: Integration (Steps 14-17)

### Step 14: Ortholog Mapping

**What:** Cross-species gene mapping between human (HGNC) and mouse (MGI).

**Files:**
- `src/hbam/integration/orthologs.py`

**Key patterns:**
- `load_ortholog_table(source: str = "biomart", cache_path: Path = None) -> pd.DataFrame`
  - Query Ensembl BioMart via pybiomart for human-mouse orthologs
  - Cache result to `data/external/ortholog_table.tsv`
  - Columns: human_gene, mouse_gene, orthology_type, sequence_identity
- `map_orthologs(genes: list[str], source_species: str, target_species: str, table: pd.DataFrame, strategy: str = "best") -> dict`
  - One-to-one: direct mapping
  - One-to-many: select best by sequence identity (strategy="best") or keep all (strategy="all")
  - Log: mapped count, unmapped count, one-to-many count
- `harmonize_gene_names(adata: AnnData, species: str, table: pd.DataFrame) -> AnnData`
  - Convert all gene names to a common namespace (default: human HGNC symbols)
  - Store original names in `.var["original_gene_name"]`, species in `.var["source_species"]`

**Acceptance criteria:**
- BioMart query returns >15,000 ortholog pairs (or cached table loads)
- One-to-many resolution picks highest sequence identity
- Unmapped genes logged but not silently dropped
- Harmonized gene names are consistent across modalities

**Dependencies:** Steps 11-13

**Tests:** `tests/test_integration.py::test_ortholog_mapping`

---

### Step 15: Gene Alignment

**What:** Compute gene intersection across modalities, validate overlap.

**Files:**
- `src/hbam/integration/align.py`

**Key patterns:**
- `align_genes(mudata: MuData, min_overlap: int = 1000) -> MuData`
  - Compute intersection of `.var_names` across all modalities
  - Log feature counts at each step: per-modality count, intersection count, retention fraction
  - Subset each modality to shared genes
  - Run `validate_gene_overlap()` gate -- error if < min_overlap
- `report_alignment(mudata: MuData) -> dict`
  - Return: per-modality gene counts, intersection size, Jaccard similarity, genes unique to each modality

**Acceptance criteria:**
- Intersection computed correctly (verified on known test sets)
- MuData modalities all have identical `.var_names` after alignment
- Error raised if intersection < 1000 genes
- Alignment report contains correct statistics
- Retains >= 50% of smaller modality's genes (logged)

**Dependencies:** Step 14

**Tests:** `tests/test_integration.py::test_gene_alignment`

---

### Step 16: Feature Scaling

**What:** Standardize features across modalities to comparable ranges.

**Files:**
- `src/hbam/integration/align.py` (extend)

**Key patterns:**
- `scale_features(mudata: MuData, method: str = "zscore") -> MuData`
  - Z-score: per-feature mean=0, std=1 within each modality
  - Min-max: scale to [0, 1] per feature (alternative)
  - Store unscaled in `.layers["unscaled"]`
- Apply after gene alignment, before embedding

**Acceptance criteria:**
- After z-score scaling, each feature has mean ~0 and std ~1
- Unscaled values preserved in layers
- No NaN introduced by scaling (check for zero-variance features first, remove or warn)

**Dependencies:** Step 15

**Tests:** `tests/test_integration.py::test_feature_scaling`

---

### Step 17: Latent Embedding (MOFA+ / PCA)

**What:** Multi-omics factor analysis for joint dimensionality reduction.

**Files:**
- `src/hbam/integration/embed.py`

**Key patterns:**
- `embed_mofa(mudata: MuData, n_factors: int = 15, seed: int = 42) -> MuData`
  - Configure MOFA+ via mofapy2: set number of factors, convergence criteria, seed
  - Fit model, store factor values in `.obsm["mofa_factors"]`
  - Store factor loadings (weights) in `.varm["mofa_weights"]` per modality
  - Store variance explained in `.uns["mofa_variance_explained"]`
- `embed_pca(mudata: MuData, n_components: int = 15, seed: int = 42) -> MuData`
  - Concatenate modalities, run PCA via sklearn
  - Store in `.obsm["pca_factors"]`, loadings in `.varm["pca_loadings"]`
- `auto_select_embedding(mudata: MuData, config: IntegrationConfig) -> str`
  - Use MOFA+ if: >1 modality with overlapping samples AND mofapy2 importable
  - Fallback to PCA otherwise
  - Log selection reason
- `run_embedding(mudata: MuData, config: IntegrationConfig) -> MuData`
  - Dispatch to MOFA+ or PCA based on auto-selection or config override

**Acceptance criteria:**
- MOFA+ produces n_factors factors (or PCA n_components components)
- Factor values stored correctly in `.obsm`
- Loadings stored in `.varm` with correct dimensions
- Silhouette score on biological groups > 0.3 (on test data with injected structure)
- PCA fallback works when MOFA+ unavailable
- Seed produces identical results across runs

**Dependencies:** Step 16

**Tests:** `tests/test_integration.py::test_embed_mofa`, `test_embed_pca`, `test_auto_select`

---

## Phase 5: HBAM Engine (Steps 18-22)

### Step 18: Literature Gene Set Loading

**What:** Load and parse curated aging/muscle gene sets.

**Files:**
- `src/hbam/engine/gene_sets.py`

**Key patterns:**
- `load_genage(path: Path = None) -> set[str]`
  - Parse GenAge database CSV (download or from data/external/)
  - Return set of human gene symbols associated with aging
- `load_msigdb_hallmark(collection: str, gene_set_names: list[str] = None) -> dict[str, set[str]]`
  - Use gseapy to fetch MSigDB hallmark gene sets
  - Default sets: HALLMARK_MYOGENESIS, HALLMARK_OXIDATIVE_PHOSPHORYLATION, HALLMARK_INFLAMMATORY_RESPONSE, HALLMARK_P53_PATHWAY, HALLMARK_MTORC1_SIGNALING
  - Cache locally in data/external/
- `load_custom_gene_set(path: Path) -> set[str]`
  - Load from text file (one gene per line) or GMT format
- `categorize_genes(gene_sets: dict) -> dict[str, set[str]]`
  - Group into "maturation" (developmental, growth) and "dysfunction" (senescence, inflammation, atrophy)
  - Configurable category assignments

**Acceptance criteria:**
- GenAge loads >300 genes
- MSigDB hallmark sets load with expected gene counts
- Custom gene set loading works for both formats
- Categories contain non-overlapping gene sets (or overlap is logged)

**Dependencies:** Step 1

**Tests:** `tests/test_engine.py::test_load_gene_sets`

---

### Step 19: Data-Driven Gene Set Refinement

**What:** Use trajectory/clustering on latent space to refine gene set membership.

**Files:**
- `src/hbam/engine/gene_sets.py` (extend)

**Key patterns:**
- `refine_gene_sets(mudata: MuData, literature_sets: dict, method: str = "trajectory") -> dict[str, set[str]]`
- Trajectory approach: 
  - Use diffusion pseudotime or PCA-based ordering on latent factors
  - Correlate gene expression with pseudotime
  - Genes positively correlated with aging trajectory -> dysfunction
  - Genes negatively correlated -> maturation
- Clustering approach:
  - Cluster genes by expression pattern (k-means on loadings)
  - Assign clusters to maturation/dysfunction based on overlap with literature sets
- Hybrid merge: union of literature + data-driven, with source tracking
- Store: gene -> source ("literature", "data_driven", "both") mapping

**Acceptance criteria:**
- Refinement adds genes not in original literature sets
- Source tracking correctly labels each gene
- At least 50% of literature genes retained after refinement
- Refined sets are non-empty for both categories

**Dependencies:** Steps 17, 18

**Tests:** `tests/test_engine.py::test_refine_gene_sets`

---

### Step 20: Gene Weight Computation

**What:** Compute per-gene weights from latent factor loadings.

**Files:**
- `src/hbam/engine/weights.py`

**Key patterns:**
- `compute_weights_mofa(mudata: MuData, gene_sets: dict, factors: list[int] = None) -> pd.DataFrame`
  - Extract MOFA+ loadings for HBAM-relevant genes
  - Weight = loading magnitude on aging-associated factors
  - Auto-select factors: those with highest correlation to age/condition metadata
  - Return DataFrame: gene, weight, factor, category (maturation/dysfunction), modality
- `compute_weights_pca(mudata: MuData, gene_sets: dict, components: list[int] = None) -> pd.DataFrame`
  - Same logic using PCA loadings
- `normalize_weights(weights: pd.DataFrame) -> pd.DataFrame`
  - Scale weights to sum to 1 within each category
  - Preserve sign (positive for contributing, negative for anti-contributing)

**Acceptance criteria:**
- Weights DataFrame has expected columns and no NaN
- Weights sum to ~1 per category after normalization
- Factor selection correlates with biological metadata
- All genes in gene sets have assigned weights

**Dependencies:** Step 19

**Tests:** `tests/test_engine.py::test_compute_weights`

---

### Step 21: HBAM Scoring

**What:** Compute directional HBAM composite score per sample.

**Files:**
- `src/hbam/engine/score.py`

**Key patterns:**
- `compute_hbam_score(mudata: MuData, weights: pd.DataFrame) -> pd.Series`
  - Dysfunction score: weighted sum of dysfunction genes per sample
  - Maturation score: weighted sum of maturation genes per sample
  - HBAM = dysfunction_score - maturation_score
  - Higher HBAM = more aging burden
  - Store in `.obs["hbam_score"]`, `.obs["dysfunction_score"]`, `.obs["maturation_score"]`
- `validate_score_direction(scores: pd.Series, metadata: pd.DataFrame, age_col: str) -> bool`
  - Check that HBAM positively correlates with known aging markers/conditions
  - If negative correlation, offer sign-flip option (log warning)
- Return scores as pd.Series indexed by sample ID

**Acceptance criteria:**
- HBAM score computed for all samples (no NaN)
- Score direction matches expected biology (positive correlation with age/condition)
- Scores stored in `.obs`
- Wilcoxon p < 0.05 between conditions on test data with injected effect

**Dependencies:** Step 20

**Tests:** `tests/test_engine.py::test_hbam_scoring`

---

### Step 22: Bootstrap Validation

**What:** Validate HBAM weight stability via bootstrap resampling.

**Files:**
- `src/hbam/engine/score.py` (extend)

**Key patterns:**
- `bootstrap_validate(mudata: MuData, gene_sets: dict, n_iter: int = 5, fraction: float = 0.8, seed: int = 42) -> BootstrapResult`
  - For each iteration: subsample 80% of samples, recompute weights, recompute HBAM scores
  - Compute Spearman correlation between each bootstrap's weights and full-data weights
  - `BootstrapResult` dataclass: correlations (list[float]), mean_correlation, passed (bool), threshold (float)
  - Pass criterion: all Spearman correlations > 0.8
- `bootstrap_confidence_intervals(mudata: MuData, gene_sets: dict, n_iter: int = 100) -> pd.DataFrame`
  - Compute 95% CI for each gene weight
  - Return DataFrame with gene, weight_mean, weight_lower, weight_upper
- Seed control: `np.random.default_rng(seed + iteration)` for each bootstrap

**Acceptance criteria:**
- Bootstrap produces n_iter results
- Spearman correlations computed correctly (validated against manual calculation)
- Pass/fail correctly determined by threshold
- Confidence intervals bracket the full-data weight
- Deterministic results with same seed

**Dependencies:** Step 21

**Tests:** `tests/test_engine.py::test_bootstrap_validation`

---

## Phase 6: Translation (Steps 23-25)

### Step 23: Sample-Level HBAM Scoring

**What:** Apply HBAM weights to individual samples for downstream analysis.

**Files:**
- `src/hbam/translate/sample_score.py`

**Key patterns:**
- `score_samples(mudata: MuData, weights: pd.DataFrame) -> pd.DataFrame`
  - Weighted sum per sample using validated HBAM weights
  - Return DataFrame: sample_id, hbam_score, dysfunction_score, maturation_score, condition, timepoint (from metadata)
- `compare_conditions(scores: pd.DataFrame, group_col: str) -> dict`
  - Wilcoxon rank-sum between groups
  - Cohen's d effect size
  - Return: p_value, effect_size, group_means, group_medians
- Store results in `.uns["hbam_comparison"]`

**Acceptance criteria:**
- All samples scored (no NaN)
- Statistical comparison runs without error
- Cohen's d > 0.8 on test data with injected large effect
- Results stored in `.uns`

**Dependencies:** Steps 21-22

**Tests:** `tests/test_translation.py::test_score_samples`

---

### Step 24: Biomarker Panel Reduction

**What:** Select top N genes as a reduced biomarker panel.

**Files:**
- `src/hbam/translate/biomarkers.py`

**Key patterns:**
- `select_biomarkers(weights: pd.DataFrame, n: int = 20, method: str = "top_weight") -> pd.DataFrame`
  - Rank genes by absolute HBAM weight
  - Select top N
  - Return DataFrame: gene, weight, rank, category, modality
- `evaluate_panel(mudata: MuData, panel_genes: list[str], group_col: str) -> dict`
  - Recompute HBAM using only panel genes
  - Compare discrimination power (AUC, sensitivity, specificity) vs full gene set
  - Return: auc, sensitivity, specificity, panel_correlation_with_full (Spearman)
- `optimize_panel_size(mudata: MuData, weights: pd.DataFrame, sizes: list[int] = [5, 10, 15, 20, 30, 50]) -> pd.DataFrame`
  - Evaluate each panel size, return performance vs size tradeoff

**Acceptance criteria:**
- Panel contains exactly N genes
- Panel genes are from top of ranked list
- Panel HBAM correlates > 0.9 with full HBAM (on test data)
- Optimization returns results for all requested sizes

**Dependencies:** Step 23

**Tests:** `tests/test_translation.py::test_biomarker_panel`

---

### Step 25: Clinical Metrics Simulation

**What:** Simulate clinical-relevant metrics from HBAM scores.

**Files:**
- `src/hbam/translate/clinical.py`

**Key patterns:**
- `compute_effect_sizes(scores: pd.DataFrame, comparisons: list[tuple]) -> pd.DataFrame`
  - For each pairwise comparison: Cohen's d, Hedge's g, Glass's delta
  - Return DataFrame: comparison, cohens_d, hedges_g, p_value, n1, n2
- `discrimination_analysis(scores: pd.DataFrame, group_col: str, positive_label: str) -> dict`
  - ROC-AUC (sklearn), optimal threshold (Youden's J)
  - Sensitivity, specificity at optimal threshold
- `generate_clinical_summary(scores: pd.DataFrame, effect_sizes: pd.DataFrame, discrimination: dict) -> str`
  - Formatted text summary suitable for manuscript methods/results

**Acceptance criteria:**
- Effect sizes computed correctly (validated against manual calculation for simple case)
- AUC between 0.5 and 1.0
- Clinical summary is human-readable text

**Dependencies:** Step 23

**Tests:** `tests/test_translation.py::test_clinical_metrics`

---

## Phase 7: Output (Steps 26-30)

### Step 26: Publication Style System

**What:** Consistent styling for all figures.

**Files:**
- `src/hbam/output/style.py`

**Key patterns:**
- `set_publication_style() -> None`
  - Set matplotlib rcParams: font sizes (8pt labels, 10pt titles), font family (Arial/Helvetica), line widths
  - Figure dimensions: single column (3.5"), double column (7"), full page (7" x 9")
- `PALETTE: dict` -- colorblind-safe palette (based on Wong 2011 or Okabe-Ito)
  - Named colors: "aging", "young", "dysfunction", "maturation", "modality_1", "modality_2", "modality_3"
- `save_figure(fig, name: str, output_dir: Path, formats: list = ["pdf", "png"], dpi: int = 300)`
  - Save in both PDF (vector) and PNG (raster) formats
  - Tight layout applied automatically
- `get_condition_colors(conditions: list) -> dict` -- map conditions to palette colors

**Acceptance criteria:**
- Style applies consistently across matplotlib and seaborn
- Palette passes colorblind simulation (deuteranopia check)
- PDF output is vector (text selectable)
- PNG output at 300 DPI

**Dependencies:** Step 1

**Tests:** `tests/test_output.py::test_style_system`

---

### Step 27: QC and Diagnostic Figures (Figures 2, 10)

**What:** QC summary and cross-modality diagnostic figures.

**Files:**
- `src/hbam/output/figures.py`

**Key patterns:**
- `fig_qc_summary(mudata: MuData, output_dir: Path) -> Path`
  - Multi-panel figure: (A) missingness heatmap per modality, (B) intensity distribution boxplots per sample, (C) sample correlation heatmap, (D) PCA of raw data colored by batch/condition
- `fig_integration_diagnostic(mudata: MuData, output_dir: Path) -> Path`
  - Multi-panel: (A) variance explained per factor/component, (B) factor-modality association heatmap, (C) sample mixing across modalities in latent space, (D) gene loading distributions per modality

**Acceptance criteria:**
- Figures render without error on test data
- Multi-panel layout has no overlapping text
- Saved as both PDF and PNG
- Figure dimensions appropriate for journal submission

**Dependencies:** Steps 8-10, 17, 26

**Tests:** `tests/test_output.py::test_fig_qc`, `test_fig_diagnostic`

---

### Step 28: Latent Space and Weight Figures (Figures 3, 4)

**What:** PCA/UMAP visualization of latent space and gene weight heatmap.

**Files:**
- `src/hbam/output/figures.py` (extend)

**Key patterns:**
- `fig_latent_space(mudata: MuData, output_dir: Path, method: str = "umap") -> Path`
  - Compute UMAP from latent factors (scanpy.tl.umap on `.obsm` factors)
  - Three panels colored by: (A) condition/age, (B) modality, (C) HBAM score (continuous colormap)
- `fig_weight_heatmap(weights: pd.DataFrame, output_dir: Path, top_n: int = 50) -> Path`
  - Clustered heatmap of top N gene weights across factors
  - Row annotation: category (maturation/dysfunction), modality
  - Column annotation: factor importance
  - Use seaborn.clustermap with diverging colormap (blue-white-red)

**Acceptance criteria:**
- UMAP computes from latent factors without error
- Heatmap shows clear clustering of maturation vs dysfunction genes (on structured test data)
- Color legends present and readable

**Dependencies:** Steps 17, 20, 26

**Tests:** `tests/test_output.py::test_fig_latent`, `test_fig_heatmap`

---

### Step 29: Score and Enrichment Figures (Figures 5, 6, 7)

**What:** HBAM score distribution, pathway enrichment, and volcano plot.

**Files:**
- `src/hbam/output/figures.py` (extend)

**Key patterns:**
- `fig_hbam_distribution(scores: pd.DataFrame, output_dir: Path) -> Path`
  - Violin/box plot of HBAM scores by condition
  - Annotate with p-value and effect size
  - Individual data points overlaid (jittered strip plot)
- `fig_pathway_enrichment(enrichment_results: pd.DataFrame, output_dir: Path) -> Path`
  - Dot plot: x=gene ratio, y=pathway name, size=gene count, color=-log10(FDR)
  - Top 20 pathways by significance
  - Run gseapy.enrichr or gseapy.prerank on HBAM-weighted gene list
- `fig_volcano(de_results: pd.DataFrame, output_dir: Path) -> Path`
  - log2FC vs -log10(p-value)
  - Color by significance: non-significant (gray), up (red), down (blue)
  - Label top 10 genes by significance

**Acceptance criteria:**
- Violin plot shows condition separation with statistics annotation
- Pathway enrichment includes >= 3 aging/muscle pathways at FDR < 0.05 (on realistic test data)
- Volcano plot labels do not overlap (adjust text)

**Dependencies:** Steps 11-13, 21, 26

**Tests:** `tests/test_output.py::test_fig_scores`, `test_fig_enrichment`, `test_fig_volcano`

---

### Step 30: Correlation and Biomarker Figures (Figures 1, 8, 9) and Tables

**What:** Pipeline overview, correlation matrix, biomarker panel, and summary tables.

**Files:**
- `src/hbam/output/figures.py` (extend)
- `src/hbam/output/tables.py`

**Key patterns:**
- `fig_pipeline_overview(output_dir: Path) -> Path`
  - Schematic created with matplotlib patches/arrows showing data flow
  - Modalities -> Integration -> HBAM Engine -> Scores
  - Annotate with sample counts and gene counts at each stage
- `fig_correlation_matrix(mudata: MuData, scores: pd.DataFrame, output_dir: Path) -> Path`
  - Correlation between HBAM components (dysfunction, maturation, individual factors)
  - Annotated heatmap with hierarchical clustering
- `fig_biomarker_panel(panel: pd.DataFrame, mudata: MuData, output_dir: Path) -> Path`
  - Heatmap of top N biomarker genes across samples
  - Samples ordered by HBAM score, annotated by condition
- Tables:
  - `generate_summary_table(mudata, scores, weights) -> pd.DataFrame` -- pipeline summary statistics
  - `generate_gene_table(weights) -> pd.DataFrame` -- full gene weight table for supplementary
  - Save as both CSV and formatted LaTeX

**Acceptance criteria:**
- All 10 figures render and save as PDF
- Tables contain expected columns and no NaN
- Pipeline overview is readable and correctly shows data flow
- LaTeX table compiles without error

**Dependencies:** Steps 21-24, 26

**Tests:** `tests/test_output.py::test_fig_overview`, `test_fig_correlation`, `test_fig_biomarker`, `test_tables`

---

## Phase 8: Pipeline Assembly (Steps 31-33)

### Step 31: CLI Entry Point

**What:** Typer-based CLI for running the full pipeline or individual modules.

**Files:**
- `src/hbam/cli.py`

**Key patterns:**
- `app = typer.Typer()`
- `@app.command("run")` -- full pipeline: load config, execute all steps, save outputs
  - `--config PATH` (required): path to YAML config
  - `--output-dir PATH` (default: results/): output directory
  - `--seed INT` (default from config): override seed
  - `--dry-run` flag: validate config and data, report what would run
- `@app.command("qc")` -- run only data loading + QC
- `@app.command("integrate")` -- run only integration (assumes processed data exists)
- `@app.command("score")` -- run only HBAM scoring (assumes integration done)
- `@app.command("figures")` -- regenerate figures from saved intermediate results
- Each command: setup logging, set seeds (numpy, random, torch if available), run steps, save results

**Acceptance criteria:**
- `python -m hbam run --config config/default.yaml` executes full pipeline (on test data)
- `--dry-run` validates without executing
- Individual sub-commands work independently
- Exit code 0 on success, non-zero on failure with informative error

**Dependencies:** Steps 1-30

**Tests:** `tests/test_cli.py`

---

### Step 32: End-to-End Pipeline Orchestration

**What:** Wire all modules together in the correct execution order with validation gates.

**Files:**
- `src/hbam/pipeline.py`

**Key patterns:**
- `run_pipeline(config: PipelineConfig) -> PipelineResult`
  - Step 1: Load data (dispatch to appropriate loader per modality config)
  - Step 2: QC each modality -> validation gate
  - Step 3: Normalize each modality -> validation gate
  - Step 4: Impute missing values
  - Step 5: Run modality-specific analysis (parallelizable with concurrent.futures)
  - Step 6: Ortholog mapping + gene alignment -> validation gate (min 1000 genes)
  - Step 7: Feature scaling
  - Step 8: Latent embedding -> validation (silhouette > 0.3)
  - Step 9: Load + refine gene sets
  - Step 10: Compute weights + HBAM scores
  - Step 11: Bootstrap validation -> validation gate (Spearman > 0.8)
  - Step 12: Translation (sample scores, biomarkers, clinical metrics)
  - Step 13: Generate all figures and tables
  - Step 14: Save intermediate results (MuData to .h5mu, scores to CSV)
- `PipelineResult` dataclass: mudata, scores, weights, figures_paths, validation_report
- Save intermediate MuData after each major phase (data/, integration/, engine/)
- Total elapsed time logging

**Acceptance criteria:**
- Pipeline runs end-to-end on test data without error
- All validation gates checked
- Intermediate files saved at each phase
- Total runtime logged
- PipelineResult contains all expected fields

**Dependencies:** Steps 1-31

**Tests:** `tests/test_pipeline.py::test_end_to_end`

---

### Step 33: Integration Tests and Final Validation

**What:** End-to-end tests, edge case tests, and documentation.

**Files:**
- `tests/test_pipeline.py`
- `tests/test_edge_cases.py`
- `docs/pipeline.md`

**Key patterns:**
- Integration tests:
  - `test_end_to_end_small()` -- 10 samples, 500 genes, verify all outputs exist
  - `test_end_to_end_minimal()` -- 3 samples, 100 genes, verify graceful handling
  - `test_reproducibility()` -- run twice with same seed, assert identical outputs
  - `test_pca_fallback()` -- mock MOFA+ import failure, verify PCA path works
- Edge case tests:
  - `test_high_missingness()` -- >50% missing data, verify rejection
  - `test_single_modality()` -- only one data type, verify degraded but functional
  - `test_no_ortholog_overlap()` -- disjoint gene sets, verify informative error
  - `test_confounded_batch()` -- batch = condition, verify warning
- Documentation:
  - Pipeline overview with flow diagram
  - Installation instructions (conda environment recommended for mofapy2)
  - Configuration reference (all YAML fields)
  - Example usage
  - Troubleshooting (MOFA+ on Windows, memory issues)

**Acceptance criteria:**
- All tests pass
- >= 80% code coverage
- Pipeline runs end-to-end: `python -m hbam run --config config/default.yaml`
- All 10 figures rendered as valid PDFs
- Documentation covers installation, configuration, and usage

**Dependencies:** Steps 1-32

**Tests:** `tests/test_pipeline.py`, `tests/test_edge_cases.py`

---

## Success Criteria (Global)

1. Each input file loads without error, matrix dimensions match metadata
2. After normalization, median log2-intensity per sample within 1.5-fold range
3. Gene alignment retains >= 50% of smaller modality's genes (logged)
4. Latent embedding: silhouette score on biological groups > 0.3
5. HBAM weights: bootstrap Spearman > 0.8 across 5 iterations
6. HBAM scores: Wilcoxon p < 0.05 between known conditions, Cohen's d > 0.8
7. Pathway enrichment: >= 3 aging/muscle-related pathways at FDR < 0.05
8. All 10 figures render as valid PDFs with colorblind-safe palette
9. Pipeline runs end-to-end with `python -m hbam run --config config/default.yaml`
10. All tests pass, >= 80% coverage on new code

## Risks and Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| Ortholog mapping ambiguity (one-to-many) | Medium | Best-ortholog by sequence identity; log all mappings; configurable strategy |
| Batch-biology confounding | High | Check design matrix before batch correction; warn and skip if fully confounded |
| MOFA+ installation issues on Windows | Medium | PCA fallback; document conda environment setup; auto-select embedding |
| Insufficient gene overlap (<1000) | High | Log aggressively; diagnostic figure; warn user; configurable minimum |
| Zero-inflated spatial data mishandled | High | NB-based methods; never impute biological zeros; separate MNAR/MCAR paths |
| HBAM score direction arbitrary | Medium | Validate against known aging markers; sign-flip option; direction check in validation |

## Dependency Graph (Condensed)

```
Phase 1: [1] -> [2] -> [3,4] -> [5]
Phase 2: [5] -> [6] -> [7] -> [8] -> [9] -> [10]
Phase 3: [10] -> [11,12,13]  (parallel)
Phase 4: [11,12,13] -> [14] -> [15] -> [16] -> [17]
Phase 5: [1] -> [18], [17,18] -> [19] -> [20] -> [21] -> [22]
Phase 6: [22] -> [23] -> [24,25]  (parallel)
Phase 7: [1] -> [26], [various] -> [27,28,29,30]  (parallel after deps met)
Phase 8: [all] -> [31] -> [32] -> [33]
```

## Estimated File Count

- Source files: ~30 (src/hbam/**)
- Test files: ~10 (tests/**)
- Config files: ~3 (pyproject.toml, config/default.yaml, .gitignore)
- Documentation: ~2 (README.md, docs/pipeline.md)
- **Total: ~45 files**
