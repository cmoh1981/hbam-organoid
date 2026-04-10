# HBAM Pipeline Documentation

## Overview

The HBAM (High Burden Aging Muscle) pipeline is a reproducible Python tool for computing a composite aging index from multi-omics data. It integrates proteomics and spatial transcriptomics from mouse tissue and human iPSC-derived organoids.

## Installation

```bash
# Clone and install
git clone <repo-url>
cd hbam-organoid
pip install -e ".[all]"

# For MOFA+ support (recommended: use conda)
conda install -c bioconda mofapy2
```

### Requirements
- Python >= 3.10
- Key dependencies: mudata, anndata, scanpy, gseapy, scikit-learn, matplotlib, seaborn

## Quick Start

```bash
# Run with default config (generates demo data if no real data configured)
python -m hbam run --config config/default.yaml

# Validate config without running
python -m hbam run --config config/default.yaml --dry-run

# Run individual steps
python -m hbam qc --config config/default.yaml
python -m hbam integrate --config config/default.yaml
python -m hbam score --config config/default.yaml
```

## Configuration Reference

All parameters are in `config/default.yaml`. Key sections:

| Section | Description | Key Parameters |
|---------|-------------|----------------|
| `seed` | Random seed | Default: 42 |
| `data.modalities` | Input data files | name, path, format, species, modality_type |
| `qc` | Quality control | missingness_threshold (0.3), variance_percentile (10) |
| `normalization` | Normalization | proteomics_method (median), spatial_method (scran) |
| `imputation` | Missing values | method (auto), mnar_method (minprob), mcar_method (knn) |
| `integration` | Multi-omics | embedding_method (auto), n_factors (15), min_gene_overlap (1000) |
| `engine` | HBAM scoring | gene_set_source (hybrid), bootstrap_iterations (5) |
| `translation` | Clinical metrics | top_biomarkers (20) |
| `output` | Figures/tables | figure_formats ([pdf, png]), figure_dpi (300) |

### Configuring Data Sources

```yaml
data:
  modalities:
    - name: liver_proteomics
      path: data/raw/proteinGroups.txt
      format: maxquant          # maxquant, diann, csv, tsv, gem
      species: mouse            # mouse or human
      modality_type: temporal   # temporal, spatial, functional
      time_col: timepoint
      condition_col: condition
    - name: spatial_muscle
      path: data/raw/stereo_seq.gem
      format: gem
      species: mouse
      modality_type: spatial
```

## Pipeline Architecture

```
Data Layer -> Modality Analysis -> Integration -> HBAM Engine -> Translation -> Output
```

1. **Data Layer**: Multi-format ingestion, QC filtering, normalization, imputation
2. **Modality Analysis**: Temporal DE, spatial variable genes, functional annotation
3. **Integration**: Cross-species ortholog mapping, gene alignment, MOFA+/PCA embedding
4. **HBAM Engine**: Gene set loading, weight computation, HBAM = dysfunction - maturation
5. **Translation**: Sample scoring, biomarker panel, clinical effect sizes
6. **Output**: 10 publication-ready figures, summary tables

## HBAM Formula

```
HBAM_score = dysfunction_score - maturation_score
```

- Higher HBAM = greater aging burden
- Gene weights derived from latent factor loadings
- Bootstrap validated (Spearman > 0.8 required)

## Output Files

### Figures (results/figures/)
1. Pipeline overview schematic
2. QC summary (missingness, distributions)
3. Latent space (UMAP/PCA by condition)
4. Gene weight heatmap (top 50)
5. HBAM score distribution by condition
6. Pathway enrichment dot plot
7. Volcano plot
8. Score correlation matrix
9. Biomarker panel heatmap
10. Integration diagnostic

### Tables (results/tables/)
- `summary.csv` / `summary.tex` — Pipeline statistics
- `gene_weights.csv` / `gene_weights.tex` — Full gene weight table

## Troubleshooting

### MOFA+ not installing on Windows
Use PCA fallback (automatic). Or install via conda:
```bash
conda create -n hbam python=3.11
conda activate hbam
conda install -c bioconda mofapy2
pip install -e ".[dev]"
```

### Memory issues with large Stereo-seq data
Increase bin size in config or reduce spatial bins before loading.

### Low gene overlap warning
Check species annotation. Use `source: builtin` for ortholog table as fallback.

## Running Tests

```bash
python -m pytest tests/ -v
python -m pytest tests/ -v --cov=hbam --cov-report=term-missing
```
