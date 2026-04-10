# Design: Spatial Data Integration for Muscle Aging

**Feature**: spatial-data-integration
**Created**: 2026-04-10
**Status**: ACTIVE
**PDCA Phase**: Design
**References**: [Plan](../../01-plan/features/spatial-data-integration.plan.md)

---

## 1. Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                   EXISTING PIPELINE (unchanged)                  │
│  load_data → QC → normalize → impute → modality_analysis        │
│                                                                  │
│  MODALITY 1: muscle proteomics (PXD047296, 40 samples)          │
│  MODALITY 2: muscle spatial    (Walter 2024, NEW)     ◄── NEW   │
│                                                                  │
│              ┌──────────────────────┐                            │
│              │   INTEGRATION LAYER  │                            │
│              │  gene alignment      │  ← already implemented    │
│              │  feature scaling     │  ← already implemented    │
│              │  MOFA+ / PCA embed   │  ← already implemented    │
│              └──────────┬───────────┘                            │
│                         ▼                                        │
│              ┌──────────────────────┐                            │
│              │   HBAM ENGINE        │                            │
│              │  multi-modal weights │  ← ENHANCEMENT NEEDED     │
│              │  per-spot scoring    │  ← NEW                    │
│              └──────────┬───────────┘                            │
│                         ▼                                        │
│              ┌──────────────────────┐                            │
│              │   OUTPUT LAYER       │                            │
│              │  spatial HBAM map    │  ← NEW FIGURE             │
│              │  region biomarkers   │  ← NEW TABLE              │
│              └──────────────────────┘                            │
└─────────────────────────────────────────────────────────────────┘
```

### What Already Works (No Changes Needed)
- `load_stereo_gem()` — Stereo-seq GEM loader → works for GEM format
- `load_matrix()` — generic CSV/TSV loader → works for count matrices
- `run_spatial_analysis()` — spatial QC, HVG, Moran's I
- `aggregate_pseudobulk()` — collapse spots to regions
- `align_genes()` — gene intersection across modalities
- `scale_features()` — z-score/minmax normalization
- `embed_pca()` — PCA fallback for joint embedding
- `normalize_spatial()` — scran/total normalization

### What Needs Enhancement
- `run_embedding()` — handle different sample counts across modalities
- `compute_hbam_score()` — support per-spot scoring on spatial data
- `pipeline.py` — wire multi-modality config properly

### What's New
- `scripts/download_spatial.py` — download and prepare Walter 2024 data
- `src/hbam/output/spatial_figures.py` — spatial HBAM map figure
- `config/muscle_spatial.yaml` — multi-modal config
- `load_h5ad()` — h5ad loader (if data is in AnnData format)

## 2. Data Flow Design

### Input Data Formats

| Modality | Source | Format | Loader | Samples | Genes |
|----------|--------|--------|--------|---------|-------|
| Proteomics | PXD047296 | DIA-NN pg_matrix | `load_diann_matrix()` | 40 | ~4,300 |
| Spatial | Walter 2024 | h5ad or count matrix | `load_h5ad()` or `load_matrix()` | 3+ sections | ~15,000+ |

### Processing Pipeline

```
Step 1: Load both modalities
  proteomics → AnnData (40 samples × 4,299 genes)
  spatial    → AnnData (N spots × M genes, .obsm["spatial"])

Step 2: QC each independently
  proteomics → filter_missingness → filter_variance → 4,299 genes
  spatial    → run_spatial_analysis (min_counts, min_genes) → filtered spots

Step 3: Normalize each independently
  proteomics → normalize_proteomics (median + log2)
  spatial    → normalize_spatial (scran/total + log1p)

Step 4: Impute (proteomics only — spatial zeros are biological)
  proteomics → impute_auto (MNAR detected)
  spatial    → NO imputation

Step 5: Modality analysis
  proteomics → run_temporal_analysis (Kruskal-Wallis across timepoints)
  spatial    → run_spatial_analysis (HVG, Moran's I)

Step 6: Integration
  6a. Pseudo-bulk spatial data by region/condition for sample-level alignment
  6b. Gene alignment (intersection ≥ 1,000 genes)
  6c. Feature scaling (z-score)
  6d. Joint embedding (MOFA+ or PCA)

Step 7: HBAM scoring
  7a. Multi-modal gene weights from joint embedding
  7b. Aggregate HBAM scores (weighted by modality contribution)
  7c. Per-spot HBAM scores on spatial data (project weights back)

Step 8: Output
  8a. All existing figures (updated with multi-modal data)
  8b. NEW: Spatial HBAM map (per-spot scores on tissue coordinates)
  8c. NEW: Region-specific aging biomarker table
```

## 3. Component Design

### 3.1 New Loader: `load_h5ad()` in `src/hbam/data/loaders.py`

```python
def load_h5ad(path: Path) -> ad.AnnData:
    """Load pre-processed AnnData h5ad file.
    
    Preserves spatial coordinates in .obsm["spatial"] if present.
    Stores raw counts in .layers["raw"].
    """
```

**Implementation**: ~15 lines. Read with `ad.read_h5ad()`, ensure `.layers["raw"]` exists, log metrics.

### 3.2 New Figure: `fig_spatial_hbam_map()` in `src/hbam/output/spatial_figures.py`

```python
def fig_spatial_hbam_map(
    adata: ad.AnnData,
    score_col: str = "hbam_score",
    output_dir: Path,
) -> list[Path]:
    """Spatial scatter plot of HBAM scores on tissue coordinates.
    
    Panels:
    A) HBAM score (diverging colormap, blue=young, red=aged)
    B) Dysfunction score
    C) Maturation score
    D) Top dysfunction gene expression
    """
```

**Implementation**: ~80 lines. Use `ax.scatter(coords[:, 0], coords[:, 1], c=scores, cmap="RdBu_r")`.

### 3.3 New Figure: `fig_spatial_gene_overlay()` in `src/hbam/output/spatial_figures.py`

```python
def fig_spatial_gene_overlay(
    adata: ad.AnnData,
    genes: list[str],
    output_dir: Path,
) -> list[Path]:
    """Overlay top aging gene expression on spatial coordinates.
    
    2x3 grid showing top 6 genes from HBAM weights on tissue.
    """
```

**Implementation**: ~50 lines. Grid of scatter plots, one per gene.

### 3.4 Enhanced: Per-Spot HBAM Scoring in `src/hbam/engine/score.py`

```python
def compute_spatial_hbam_score(
    adata: ad.AnnData,
    weights: pd.DataFrame,
) -> ad.AnnData:
    """Compute HBAM score for each spatial spot.
    
    Unlike sample-level scoring, this applies gene weights to each spot
    in the spatial data, creating a continuous spatial aging map.
    
    Stores: .obs["hbam_score"], .obs["dysfunction_score"], .obs["maturation_score"]
    """
```

**Implementation**: ~40 lines. Same weighted-sum logic as `compute_hbam_score()` but applied to spatial AnnData directly.

### 3.5 New: Region Biomarker Analysis in `src/hbam/translate/spatial_biomarkers.py`

```python
def identify_region_biomarkers(
    adata: ad.AnnData,
    score_col: str = "hbam_score",
    n_regions: int = 4,
    top_n: int = 10,
) -> pd.DataFrame:
    """Identify genes enriched in high-HBAM vs low-HBAM spatial regions.
    
    1. Bin spots into high/low HBAM regions (quartiles)
    2. DE between high-HBAM and low-HBAM spots
    3. Return top genes driving spatial aging pattern
    """
```

**Implementation**: ~60 lines. Quartile split + Mann-Whitney per gene.

### 3.6 Enhanced Config: `config/muscle_spatial.yaml`

```yaml
data:
  modalities:
    - name: muscle_proteomics
      path: data/raw/PXD047296/MuscleDIANNStep3report.pg_matrix.tsv
      format: diann
      species: mouse
      modality_type: temporal
      condition_col: condition
    - name: muscle_spatial
      path: data/raw/walter2024/spatial_muscle.h5ad  # or .gem or .csv
      format: h5ad
      species: mouse
      modality_type: spatial

integration:
  min_gene_overlap: 500  # Lower for cross-platform
  embedding_method: auto  # MOFA+ if available
  n_factors: 10
```

### 3.7 Pipeline Enhancement: `src/hbam/pipeline.py`

Changes needed in `_load_all_modalities()`:
- Add `"h5ad": load_h5ad` to loader_map

Changes needed in `_generate_outputs()`:
- After standard figures, if any modality has `.obsm["spatial"]`:
  - Compute per-spot HBAM scores
  - Generate `fig_spatial_hbam_map()`
  - Generate `fig_spatial_gene_overlay()`
  - Generate region biomarker table

## 4. File Changes Summary

| File | Change Type | Lines (~) | Description |
|------|-------------|-----------|-------------|
| `src/hbam/data/loaders.py` | Add function | +20 | `load_h5ad()` |
| `src/hbam/engine/score.py` | Add function | +45 | `compute_spatial_hbam_score()` |
| `src/hbam/output/spatial_figures.py` | **New file** | +150 | `fig_spatial_hbam_map()`, `fig_spatial_gene_overlay()` |
| `src/hbam/translate/spatial_biomarkers.py` | **New file** | +80 | `identify_region_biomarkers()` |
| `src/hbam/pipeline.py` | Modify | +30 | Add h5ad loader, spatial outputs |
| `config/muscle_spatial.yaml` | **New file** | +50 | Multi-modal config |
| `scripts/download_spatial.py` | **New file** | +100 | Download Walter 2024 data |
| `tests/test_spatial_integration.py` | **New file** | +80 | Integration tests |
| **Total** | | ~555 | |

## 5. Implementation Order

```
1. scripts/download_spatial.py          (data acquisition - independent)
2. src/hbam/data/loaders.py             (add load_h5ad - independent)
3. src/hbam/engine/score.py             (add spatial scoring - depends on 2)
4. src/hbam/output/spatial_figures.py    (new spatial figures - depends on 3)
5. src/hbam/translate/spatial_biomarkers.py  (region analysis - depends on 3)
6. src/hbam/pipeline.py                 (wire everything - depends on 2,3,4,5)
7. config/muscle_spatial.yaml           (config - depends on 6)
8. tests/test_spatial_integration.py    (tests - depends on all)
```

**Parallelizable**: Steps 1-2 can run in parallel. Steps 3-5 can run in parallel after step 2.

## 6. Acceptance Criteria (Testable)

| # | Criterion | Test Method |
|---|-----------|-------------|
| AC1 | `load_h5ad()` loads AnnData with `.obsm["spatial"]` preserved | Unit test with fixture |
| AC2 | Spatial QC retains >50% spots, >200 genes/spot | Assertion in test |
| AC3 | Gene overlap between proteomics and spatial ≥ 500 | `validate_gene_overlap()` gate |
| AC4 | Joint embedding produces ≥5 factors with >50% variance | Check `.uns["pca_variance_explained"]` |
| AC5 | Per-spot HBAM scores computed, no NaN in `.obs["hbam_score"]` | Assertion |
| AC6 | `fig_spatial_hbam_map` renders valid PDF with spatial pattern | File exists + non-empty |
| AC7 | Region biomarkers table has ≥10 genes with FDR < 0.05 | DataFrame check |
| AC8 | Pipeline runs with multi-modal config in <5 minutes | Timing assertion |
| AC9 | All 52 existing tests still pass | `pytest tests/` regression |
| AC10 | At least 2 new spatial-specific tests pass | New test file |

## 7. Risk Mitigations

| Risk | Design Mitigation |
|------|-------------------|
| Walter data unavailable | `load_h5ad()` also works with any h5ad; fallback to STDS0000247 liver |
| Different sample counts in MOFA+ | `aggregate_pseudobulk()` reduces spatial to sample-level first |
| Too many spots (>1M) | Configurable bin_size in spatial loader; subsample if needed |
| Gene overlap too low | Lower `min_gene_overlap` to 500; log diagnostic |
| Spatial coords missing | Check `.obsm["spatial"]` exists; skip spatial figures if absent |
