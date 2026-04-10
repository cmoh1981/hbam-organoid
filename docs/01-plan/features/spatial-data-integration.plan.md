# Plan: Spatial Data Integration for Muscle Aging

**Feature**: spatial-data-integration
**Created**: 2026-04-10
**Status**: ACTIVE
**PDCA Phase**: Plan

---

## 1. Objective

Integrate spatial transcriptomics data into the HBAM pipeline to add tissue-level spatial resolution to the aging analysis. This enables answering WHERE in the tissue aging changes occur, not just WHAT changes.

## 2. Background

The HBAM pipeline currently runs on single-modality proteomics (PXD047296, mouse muscle, 40 samples). The initial analysis identified 2,452 age-changed proteins and key therapeutic targets (C3, ACLY, FASN). Adding spatial data will:

- Reveal spatial patterns of aging (e.g., are dysfunction genes enriched in specific muscle regions?)
- Enable true multi-omics integration (the pipeline's core purpose)
- Validate proteomic findings with transcriptomic evidence
- Identify cell-type-specific aging signatures

## 3. Data Sources

### Track A: Stereo-seq Liver Aging (STDS0000247)
- **Paper**: Ma et al. "Spatial transcriptomic landscape unveils immunoglobin-associated senescence" (Cell, 2024)
- **Technology**: Stereo-seq (BGI)
- **Tissues**: 9 organs including **liver** (no skeletal muscle)
- **Ages**: 2M (young), 4M, 13M (middle), 19M, 25M (old) — C57BL/6J male mice
- **Scale**: 1,535,191 spots, ~1,450 genes/spot
- **Format**: h5ad (processed), GEF
- **Download**: STOmicsDB (requires CNGB login) — `https://db.cngb.org/stomics/datasets/STDS0000247`
- **Note**: Liver spatial data can integrate with future liver proteomics from PXD047296

### Track B: Curio Seeker Muscle Aging (Walter et al. 2024) — PRIMARY
- **Paper**: Walter et al. "Transcriptomic analysis of skeletal muscle regeneration across mouse lifespan" (Nature Aging, 2024)
- **Technology**: Curio Seeker (Slide-seq-based, high resolution)
- **Tissue**: Tibialis anterior (TA) skeletal muscle — **matches our proteomics tissue**
- **Ages**: Young (~5 mo), Old (~20 mo), Geriatric (~26 mo) — C57BL/6J mice
- **Scale**: 273,923 single-cell spatial transcriptomes
- **Format**: Expected h5ad or count matrices
- **PMID**: 39578558, PMC11645289
- **Key finding**: Senescent-like muscle stem cell subsets elevated in aged injury zones

### Data Source Decision

| Criterion | Track A (Liver) | Track B (Muscle) |
|-----------|:-:|:-:|
| Tissue match with existing proteomics | No (liver vs muscle) | **Yes** (both muscle) |
| Same organism/strain | Yes (C57BL/6J) | Yes (C57BL/6J) |
| Age overlap with PXD047296 | Partial (4M overlap) | **Good** (5M≈4M, 20M=20M) |
| Spatial resolution | Sub-cellular | High (Slide-seq) |
| Data accessibility | CNGB login required | GEO expected |

**Decision: Prioritize Track B** (muscle spatial) for direct integration with muscle proteomics. Track A (liver) deferred to future cycle.

## 4. Integration Strategy

### 4.1 Cross-Modal Integration Type
- **Horizontal integration** (different samples from same tissue/organism)
- Proteomics: 40 samples at 4 timepoints (4/8/12/20 months)
- Spatial: young (~5mo) vs old (~20mo) vs geriatric (~26mo)
- Common axis: **age/condition** (young vs old), not sample-level

### 4.2 Gene Alignment
- Both datasets use mouse gene symbols — no ortholog mapping needed within this modality pair
- Existing ortholog step still runs for human gene set comparison
- Expected gene overlap: proteomics covers ~4,300 genes, spatial may cover 15,000+ → overlap likely 3,000-4,000

### 4.3 Updated HBAM Formula
- Current: `HBAM = dysfunction_score - maturation_score` (proteomics weights only)
- Updated: `HBAM = w_prot * HBAM_prot + w_spatial * HBAM_spatial` (multi-modal weighted)
- Spatial contributes: spatially variable aging genes, cell-type-specific weights
- Weight calibration: variance-explained from MOFA+ factors

## 5. Implementation Plan

### Phase 1: Data Acquisition (Steps 1-2)
1. **Locate and download Walter et al. 2024 spatial data** from GEO/paper supplements
2. **Inspect and prepare** data: convert to AnnData, validate format, create metadata

### Phase 2: Spatial Processing (Steps 3-5)
3. **Spatial QC**: filter low-quality spots, compute per-spot metrics
4. **Spatial normalization**: SCTransform or scran (already implemented)
5. **Spatially variable gene detection**: Moran's I, SpatialDE, or HVG (already implemented)

### Phase 3: Integration (Steps 6-8)
6. **Gene alignment**: intersect proteomics and spatial gene sets (no ortholog needed)
7. **Multi-modal embedding**: Run MOFA+ with both modalities (this is why we built MOFA+ support)
8. **Silhouette validation**: verify biological signal in joint latent space

### Phase 4: Enhanced HBAM (Steps 9-11)
9. **Multi-modal gene weights**: extract weights from MOFA+ factors per modality
10. **Spatial HBAM scoring**: compute per-spot HBAM scores in spatial data
11. **Spatial pattern analysis**: identify tissue regions with highest/lowest aging burden

### Phase 5: Outputs (Steps 12-14)
12. **Spatial HBAM map**: publication figure showing HBAM scores overlaid on tissue
13. **Region-specific biomarkers**: genes driving aging in specific muscle regions
14. **Updated pipeline config**: multi-modality YAML config

## 6. Acceptance Criteria

1. Spatial data loads and passes QC (>50% spots retained, >200 genes/spot)
2. Gene overlap between proteomics and spatial ≥ 1,000 genes
3. MOFA+ runs with both modalities and produces ≥5 interpretable factors
4. Per-spot HBAM scores computed for spatial data (no NaN)
5. Spatial HBAM map figure renders showing spatial aging pattern
6. Top 10 dysfunction genes from proteomics show expected spatial localization
7. Pipeline runs end-to-end with multi-modal config in <5 minutes
8. All existing tests still pass (regression-free)

## 7. Risks & Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| Walter data not on GEO | High | Contact authors, use supplement tables, or use STDS0000247 liver instead |
| Low gene overlap (<1000) | Medium | Lower min_overlap threshold, use broader gene universe |
| MOFA+ fails on mixed sample counts | Medium | PCA fallback with pseudo-bulk aggregation of spatial data |
| Spatial resolution too high (millions of spots) | Medium | Bin to lower resolution, use pseudo-bulk per region |
| Batch effects dominate biology | Medium | ComBat/Harmony batch correction before MOFA+ |

## 8. Dependencies

- Existing HBAM pipeline (completed, 52 tests passing)
- PXD047296 muscle proteomics (downloaded, processed)
- Walter et al. 2024 spatial data (to be acquired)
- MOFA+/mofapy2 or PCA fallback (implemented)

## 9. Success Metrics

- Multi-modal HBAM explains more variance than single-modal (R² improvement)
- Spatial aging pattern is biologically interpretable (matches known muscle anatomy)
- At least 2 new therapeutic insights from spatial localization (e.g., "C3 is enriched in perivascular regions → local complement inhibition")
