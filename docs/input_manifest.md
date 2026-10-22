# Input manifest (locked 2026-10-21)

Source: GEO GSE111672 (Moncada et al. 2020, Nat Biotechnol). Downloaded to
/scratch/$USER/pdac_spatial/{raw,counts}. All files are public. Data are not stored in this repo.

## Used

| Role | File(s) | Content |
|---|---|---|
| Single-cell reference, tumor A | GSE111672_PDAC-A-indrop-filtered-expMat.txt.gz | 1,926 cells x 19,737 genes (genes as rows) |
| Single-cell reference, tumor B | GSE111672_PDAC-B-indrop-filtered-expMat.txt.gz | 1,733 cells x 19,737 genes |
| ST, A1 | GSM3036911_PDAC-A-ST1-filtered.txt.gz | 428 spots |
| ST, A2, A3 | GSM4100721_PDAC-A-st2, GSM4100722_PDAC-A-st3 | 316, 288 spots |
| ST, B1 | GSM3405534_PDAC-B-ST1-filtered.txt.gz | 224 spots |
| ST, B2, B3 | GSM4100723_PDAC-B-st2, GSM4100724_PDAC-B-st3 | 321, 242 spots |
| ST, D, E, G | GSM4100725 (D), GSM4100726 (E), GSM4100728 (G) | 431 (after dedup), 359, 288 spots |

Genes x spots layout. Spot labels are array coordinates (for example 10x13). Gene lists are identical
between single-cell and spatial tables (verified for A1 by diff).

## Excluded, with reason

| File | Reason |
|---|---|
| PDAC-F ST (GSM4100727) | Truncated: 5,757 complete genes, last row cut mid-row (299 of 360 fields). IGF2BP2, SGMS2, KRT19 absent. |
| GSM3036911.tsv.gz, GSM3405534_PDAC-B-ST1.tsv.gz | Transposed unfiltered duplicates of A1 and B1. Filtered files used instead. |
| Per-library inDrop files (A, B, C) and GSM3036909, GSM3036910 | Unfiltered libraries. Pooled filtered matrices used. 909 and 910 unlabeled, treated as unused. |
| PDAC-C | Single-cell only, no spatial data. |
| H&E and Cy3 images | Not needed for the main analysis. |

## Cleaning rules

- PDAC-D: 29 duplicated spot labels, all with identical columns. Keep the first, drop the rest.
- Unfiltered arrays (A2, A3, B2, B3, D, E, G) hold 242 to 460 spots, so they appear tissue-only. This is an inference and should be checked against coordinates in the QC notebook.
- Gene symbols `1-Mar` and `2-Mar` each appear twice in every table (a spreadsheet date-conversion artifact). The two rows behind each name have different counts, so they are different genes. They are renamed `1-Mar__1`, `1-Mar__2`, `2-Mar__1`, `2-Mar__2`, are not summed or dropped, and are excluded from signatures and interpretation.

## Design decisions

1. IGF2BP2-high state: Chen et al. (2026) supplementary material contains no State 3 gene list. The state is derived by this project from malignant cells of PDAC-A and PDAC-B, following the published procedure (subcluster malignant ductal cells into three states, State 3 is the IGF2BP2-enriched one, signature genes with expression proportion >= 25 percent, weighted score).
2. Held-out testing: derive in A, test on B arrays, and the reverse. Apply the A and B signature to D, E, G as independent tumors.
3. Reference for D, E, G: cross-patient (A and B). Results there are labeled lower confidence.
4. Statistics: tumor-cell fraction from deconvolution is a covariate. Significance uses a spatially aware null, not plain spot-level Spearman. Results are reported per array and per tumor, not pooled, and treated as descriptive.

## Expected outputs (draft)

- Per-array table of deconvolved cell-type fractions per spot.
- Per-array state score and immunosuppressive and lipid signature scores per spot.
- Per-array statistics table: naive correlation, composition-adjusted association, spatial-null p-value, Moran's I of each score.
- Summary figure per tumor and a cross-tumor replication figure.
- Short report stating what held, what did not, and the limits.
