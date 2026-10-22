# Input manifest (locked 2026-10-21)

Source: GEO GSE111672 (Moncada et al. 2020, Nat Biotechnol). Downloaded to
/scratch/$USER/pdac_spatial/{raw,counts}. All files are public. Data are not stored in this repo.

## Used

| Role | File(s) | Content |
|---|---|---|
| Single-cell reference, tumor A | GSE111672_PDAC-A-indrop-filtered-expMat.txt.gz | 1,926 cells x 19,738 genes (genes as rows) |
| Single-cell reference, tumor B | GSE111672_PDAC-B-indrop-filtered-expMat.txt.gz | 1,733 cells x 19,738 genes |
| ST, A1 | GSM3036911_PDAC-A-ST1-filtered.txt.gz | 428 spots |
| ST, A2, A3 | GSM4100721_PDAC-A-st2, GSM4100722_PDAC-A-st3 | 316, 288 spots |
| ST, B1 | GSM3405534_PDAC-B-ST1-filtered.txt.gz | 224 spots |
| ST, B2, B3 | GSM4100723_PDAC-B-st2, GSM4100724_PDAC-B-st3 | 321, 242 spots |
| ST, D, E, G | GSM4100725 (D), GSM4100726 (E), GSM4100728 (G) | 429(after dedup), 359, 288 spots |

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

- PDAC-D: 29 duplicated spot labels (31 extra columns), all with identical counts. Keep the first of each, drop the rest, leaving 429 spots.
- Unfiltered arrays (A2, A3, B2, B3, D, E, G) hold 242 to 460 spots. The step 01 geometry check (see QC findings) supports treating them as tissue-only.
- Gene symbols damaged by spreadsheet date conversion: 27 rows in the PDAC-A single-cell table match a number-month pattern (1-Dec, 1-Mar to 11-Mar, 1-Sep to 14-Sep, with 13-Sep absent), and the gene lists are identical across tables. `1-Mar` and `2-Mar` each occur twice with different counts, so they are different genes. They are renamed `1-Mar__1`, `1-Mar__2`, `2-Mar__1`, `2-Mar__2` and nothing is summed or dropped. All 27 rows stay in the matrices but are excluded from signatures, marker lists, and gene-set lookups, because their labels cannot be matched to current gene symbols.

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

## QC findings (step 01, run 2026-10-22)

- Gene tables hold 19,738 genes, identical and in the same order in every file. Single-cell: 1,926 cells (A), 1,733 cells (B). Nine spatial arrays total 2,895 spots.
- Spatial depth varies about sevenfold: median UMI per spot is 403 (A3), 455 (B2), 522 (A2), 903 (B3), 1,167 (A1), 1,493 (G), 1,560 (B1), 2,356 (E), 2,994 (D).
- Spots form contiguous tissue-shaped regions (bounding-box fill 0.49 to 0.72, 99 to 100 percent of spots have a grid neighbour), so all nine arrays are treated as tissue-only. Not compared with the H&E images.
- Single-cell detection (fraction of cells with count above 0), A and B: IGF2BP2 0.31, 0.26; SGMS2 0.36, 0.34; CD274 0.015, 0.034; PTPRC 0.060, 0.038.
- Spatial detection (fraction of spots): IGF2BP2 0.04 to 0.16 (15 to 55 counts per array); SGMS2 0.00 to 0.10; CD274 0.00 to 0.01 (0 to 6 counts); PTPRC 0.00 to 0.06 (0 to 32 counts); CD68 0.15 to 0.73.
- Consequence: single-gene maps of IGF2BP2, SGMS2, CD274, and PTPRC are not interpretable. Use multi-gene signature scores and composition estimates. Immune-cell signals are close to undetectable, so T-cell exhaustion analyses are exploratory and underpowered.
- Single-cell tables contain no genes with prefix MT-, so mitochondrial percentage cannot be re-checked.

## Controls and success criterion (locked 2026-10-22, before any spatial result)

- Positive control: Moncada et al. reported colocalization of a cancer stress-response module with inflammatory fibroblasts across these arrays. The pipeline must recover that association before any result on the IGF2BP2 state is believed.
- Negative control: random gene sets matched for expression level, run through the identical statistics, to estimate the false-association rate.
- Independent units are tumors (5), not arrays (9) or spots. A1 to A3 come from one tumor, and B1 to B3 from another.
- Success: the same direction of effect in at least 4 of 5 tumors, and a pooled interval that excludes zero after adjustment for log total UMI and tumor-cell fraction. Anything weaker is reported as inconclusive.
- Order of work: crude marker-based tumor score and least-squares composition first, cell2location afterward as a sensitivity check.

## Annotation findings (step 02, 2026-10-22)

- Leiden clusters (resolution 0.8): 18 in tumor A, 16 in tumor B. Provisional labels are in config/cluster_labels.yaml.
- Epithelial cells: about 1,720 of 1,926 (A) and 1,493 of 1,733 (B).
- Candidate cancer cells: A5 (108 cells, TM4SF1 high) and A6 (168 cells, S100A4, KRT16, KRT17) match the two cancer clusters of Moncada et al. in this tumor. B0 (338 cells, TM4SF1, CDKN2A, KRT17) matches their single cancer cluster in tumor B. Marker agreement only, no copy-number evidence, so these are candidates. Clusters A8 and B10 may also contain cancer cells and are labeled epithelial.
- TM4SF1 alone does not mark cancer cells: endothelial clusters carry it too.
- IGF2BP2 detected in 54 percent of B0 cells against 16 to 29 percent of other epithelial clusters in B, and in 58 percent (A5) and 50 percent (A6) of cancer candidates against 24 to 43 percent of other epithelial clusters in A. It behaves partly as a cancer-lineage marker, so tumor-cell fraction must be a covariate in every spatial model.
- Reference limits: 7 fibroblasts (B only), 42 T and NK cells (A only, none in B). Stromal composition is estimated from spatial marker scores, and cell2location is a sensitivity check only.

## Gate for the IGF2BP2 program (locked before step 03 results)

Step 03 derives, within candidate cancer cells of each tumor, the genes that differ between IGF2BP2-detected and undetected cells (adjusted for log total UMI and a cell-cycle score) and takes the top 100. The null shuffles detection among cells of similar depth in both tumors (1,000 shuffles) and repeats the whole derivation.
- The program counts as reproducible only if both held-out AUCs (A to B and B to A) exceed 0.5 with p below 0.05 against the null, and the top-100 overlap between tumors has p below 0.05.
- If it passes: the spatial state score uses the program genes. For arrays of A and B, use the program derived in the other tumor. For D, E and G, use the union of both top-100 lists.
- If it fails: report that no reproducible IGF2BP2 program was found in these data. The spatial analysis then uses IGF2BP2 expression alone, adjusted for tumor-cell fraction and log total UMI, and is labeled exploratory.