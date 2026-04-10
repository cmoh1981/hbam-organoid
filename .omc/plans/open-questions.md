# Open Questions

## hbam-pipeline - 2026-04-09

- [ ] Which specific Stereo-seq datasets will be used (GEM vs GEF format)? — Determines whether GEF binary parser is needed in addition to GEM text parser
- [ ] Are sample metadata files (condition, timepoint, age) provided separately or embedded in data files? — Affects loader design and metadata join logic
- [ ] What is the target journal for publication figures? — Figure dimensions and style requirements vary by journal (e.g., Nature vs Cell vs PNAS)
- [ ] Should the pipeline support resuming from intermediate checkpoints, or is full re-run acceptable? — Affects whether intermediate .h5mu files need versioned saving
- [ ] Is there a preferred Python version constraint (3.10+, 3.11+, 3.12+)? — mofapy2 compatibility varies across Python versions
- [ ] Should batch correction (e.g., ComBat, Harmony) be included as an optional step, or is the current design (warn-only on confounding) sufficient? — Relevant if data has known batch effects separable from biology
- [ ] What is the expected scale of spatial data (number of bins/cells)? — Determines whether out-of-core processing (backed AnnData) is needed for the spatial module
- [ ] Are there specific aging markers (e.g., p16, p21, SA-beta-gal targets) that should be used to validate HBAM score direction? — Currently using generic correlation with age/condition metadata
