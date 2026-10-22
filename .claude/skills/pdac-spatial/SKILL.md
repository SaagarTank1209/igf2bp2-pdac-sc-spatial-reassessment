---
name: pdac-spatial
description: Use for the BINF6420 project reassessing an IGF2BP2-high malignant state in pancreatic cancer with GEO GSE111672 single-cell and spatial data. Applies locked inputs, exclusions, statistics rules, paths, and response style.
---

# Project
Test whether the IGF2BP2-high malignant state (Chen et al. 2026, npj Precis Oncol) colocalizes with
immunosuppressive and lipid metabolism signatures after deconvolution and spatially aware statistics,
and whether it repeats across tumors. Results are descriptive. Key findings due Nov 20, 2026.

# Locked inputs (see docs/input_manifest.md)
- Single-cell: PDAC-A (1,926 cells), PDAC-B (1,733 cells), pooled filtered matrices, 19,737 genes.
- Spatial: 9 arrays, 5 tumors: A1-A3, B1-B3, D, E, G. Genes x spots layout.
- Excluded: PDAC-F (truncated), PDAC-C (no spatial), transposed unfiltered A1/B1 files, per-library files.
- PDAC-D: drop 29 duplicated spot columns (identical), keep first.
- Matched reference exists for A and B only. D, E, G use a cross-patient reference, lower confidence.

# Rules
- State is derived here from A and B malignant cells. Chen supplement has no gene list.
- Test held out: derive in A and test B, and the reverse. Apply to D, E, G.
- Include tumor-cell fraction as a covariate. Never interpret a raw tumor-signature map as niche evidence.
- Use a spatially aware null. Plain spot-level Spearman is only a baseline for comparison.
- Report per array and per tumor. Do not pool. State limits in every summary.
- Never invent numbers. If a value was not computed, say so.

# Paths
- Data: /scratch/$USER/pdac_spatial/{raw,counts}. Work and results on scratch, code in the repo.
- Run df -h on the target folder before any download or environment build.
- Heavy jobs go through the scheduler, never the login node.

# Response style
Concise and formal. Small numbered steps. Ask for compact output (head, wc, summaries), not full dumps.
