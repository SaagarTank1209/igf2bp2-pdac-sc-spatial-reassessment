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

## Result of the reproducibility gate (step 03, run 2026-10-22)

Candidate cancer cells: A5 and A6 (276 cells), B0 (338 cells). The gene universe is about 3,900 genes, inferred from the chance overlap of the null.
- Top-100 overlap between tumors: 1 gene (DIP2B) against a null mean of 2.58 (p 0.935).
- AUC, program derived in A and tested in B: 0.459 (null mean 0.505, p 0.932).
- AUC, program derived in B and tested in A: 0.554 (null mean 0.500, p 0.048).
- Spearman correlation of the t statistics over all genes: 0.008 (p 0.351).
Verdict: the gate is not passed (the A to B AUC and the overlap fail). The single p of 0.048 is one of four tests and is not read as evidence. No reproducible IGF2BP2-detection program was found in these two tumors at this sample size. This is not a test of the Chen et al. State 3, whose gene list is unavailable and which was defined by subclustering and not by detection of one gene. Power of this test for realistic effect sizes is unknown.
Top genes: the tumor A list is dominated by epithelial junction and polarity genes (CLDN4, TJP3, F11R, ITGB4, EPPK1, TRIM29) and includes TM4SF1, the A5 marker, so it may partly reflect cancer-cluster identity (not checked). The tumor B list mixes signalling genes (MAPK1, YAP1, CTNNB1), RNA-processing and nucleolar genes, and a few matrix genes. Neither list is dominated by lipid-metabolism genes, though a few appear (PNPLA2 and LIPH in A, ADIPOR1 and SGPP2 in B). Enrichment was not tested.
- Diagnostics (not part of the gate): genes with t above +2 / below -2: 167 / 93 in tumor A and 135 / 102 in tumor B (about 89 expected per side under no signal, from a naive normal expectation). Median t is 0.11 (A) and 0.06 (B), mean 0.12 and 0.07. A shift of that size accounts for part of the tilt toward positive t, and the remainder fits slightly inflated t statistics. No sign of a gene program. A sensitivity run with an extra covariate for number of detected genes was considered and not done.

## Amended fallback for the spatial analysis (written before any spatial result)

Exposure: IGF2BP2 detection per spot (detected or not), modeled with log total UMI, tumor-cell fraction and a stromal score as covariates, and tested against a spatially aware null. IGF2BP2 is detected in only 4 to 16 percent of spots, so power is expected to be low and all results are exploratory. Visual single-gene maps of IGF2BP2 are still not interpreted. Only the adjusted statistical test is reported.

## Spatial analysis specification (locked 2026-10-22, before any spatial result)

Deviation from the earlier plan: the deconvolved tumor-cell fraction is replaced by a marker-based epithelial score, because the single-cell reference is too thin for deconvolution (7 fibroblasts) and cancer cells cannot be separated from normal duct cells by markers alone. The epithelial score measures epithelial content, not cancer fraction, so residual confounding by cancer fraction remains. Tumor G has little epithelium (KRT19 detected in 34 percent of spots).

Positive control, revised: Moncada et al.'s inflammatory fibroblast signature cannot be scored. IL6 is detected in at most 1 percent of spots, HAS1 in at most 4, CCL2 in at most 6, CXCL12 in at most 7 and LIF in at most 11 (CXCL1 and CXCL2 only in array D, at 63 and 29 percent). The substitute tests the stress module as listed in their Figure 5a (CEBPB, DNAJA1, FOS, FOSB, JUN, JUNB, KLF2, KLF4, MAP2K3) against the pan-stromal score (COL1A1, COL1A2, COL3A1, DCN, LUM). It checks that the pipeline can recover a spatial association and does not test the inflammatory-fibroblast claim. FOS, JUN and JUNB also respond to tissue handling, so a positive result can be technical.

Scores: mean z-score within each array of log-normalized expression (counts divided by the spot total, times 10,000, log1p) over the signature genes that vary. Epithelial covariate: KRT19, KRT7, KRT8, KRT18, MUC1, EPCAM.

Models, one per array: z-scored outcome on exposure, log total UMI, epithelial score and, for H1 and H2, the stromal score. Exposure for H1 and H2: IGF2BP2 detected in the spot (about 14 to 47 detected spots per array). Exposure for the control: stromal score, with covariates log total UMI and epithelial score.

Spatially aware null: the exposure is shuffled within 5 by 5 grid-unit blocks, 1,000 times, which keeps large-scale gradients and tests local association only. Negative control: 200 random gene sets per outcome, matched to the signature size and expression decile, go through the same procedure and give p_vs_random.

Tumors are the independent units (5). Tumor effect = mean over the tumor's arrays. An outcome meets the criterion if its effect has the same sign in at least 4 of 5 tumors, the 95 percent t interval over the 5 tumors (4 degrees of freedom) excludes zero, and the median p_vs_random over arrays is below 0.05. H1 and H2 are hypothesized positive. S1 (lipid-associated macrophages) and S2 (T-cell exhaustion) are exploratory. If the control fails, H1 and H2 are reported as unvalidated.

## Results of the spatial analysis (step 04, run 2026-10-22)

Nine arrays, five tumors. Tumor-level results (mean adjusted effect in SD of the outcome per unit of exposure, 95 percent t interval over 5 tumors):
- Control C1 (stress score against stromal score): 0 of 5 tumors positive, mean -0.046 (-0.077 to -0.015), median p_vs_random 0.592. The control did not recover a positive association, and the weak negative association it shows cannot be told apart from matched random gene sets. Control failed.
- H1 (lipid score, IGF2BP2 detected): 3 of 5 tumors positive, mean 0.068 (-0.108 to 0.245), median p_vs_random 0.562. Criterion not met.
- H2 (macrophage score): 1 of 5 tumors positive, mean -0.098 (-0.237 to 0.041), median p_vs_random 0.726. Criterion not met. The lean is negative, opposite to the hypothesis.
- S1 (lipid-associated macrophages): mean -0.026 (-0.249 to 0.196). S2 (T-cell exhaustion): mean 0.140 (-0.160 to 0.440). Exploratory, no random-set comparison.
Verdict under the locked rules: the control failed, so H1 and H2 are reported as unvalidated, and neither meets its criterion. The data are inconclusive for H1 and H2 at this sample size. They are compatible with no effect and with effects up to about a quarter of an SD.

Array level: 4 of 45 block-permutation tests are below 0.05 (A1 H1 0.012, A3 H1 0.038, B3 H1 0.037, A2 S2 0.041), and A1 H2 is at 0.050 with a negative effect. p_vs_random for the H1 and H2 hits: A1 H1 0.065, A3 H1 0.592, B3 H1 0.080, A1 H2 0.090. Three arrays show H1 effects of 0.38 to 0.47 and B2 shows -0.34.

Calibration on real data (share of matched random gene sets with block p below 0.05, expected about 0.05): control model 6 to 21 percent per array (mean 13), H1 model 4 to 28 percent (mean 9), H2 model 3 to 29 percent (mean 9). Array A3, the lowest-depth array, is the outlier for H1 and H2 (28 and 29 percent, with 62 and 46 percent of random sets at |t| above 2). The block-permutation p values are anti-conservative on these data, and p_vs_random is the stricter comparison. Possible causes (not verified): a nonlinear depth effect in A3 that the linear covariate misses, and a compositional effect from the very abundant stromal genes. By hand calculation from the array table, excluding A3 changes the H1 mean to about 0.048 and the H2 mean to about -0.087 with no change in the number of positive tumors.

Any analysis added after this point is post-hoc and labeled as such.

## Post-hoc diagnostics (plan written before step 05 was run)

All results of step 05 are post-hoc. They do not change the locked criterion or the verdict of step 04.
1. Spike-in: a shift of 0 to 0.75 SD is added to the lipid score of IGF2BP2-detected spots and the locked criterion is re-evaluated. We report the smallest shift in the grid at which the criterion is met, or "not detectable up to 0.75 SD". This is conditional on the observed noise (one realization), not an expected-power calculation. An approximate minimum detectable effect from the interval of the unshifted analysis is also reported.
2. Split-half plumbing control: the genes of a signature are split into two random halves (30 splits per array, median taken) and one half is tested as the exposure for the other, with the same covariates as the corresponding locked model. Co-variation counts as recovered if the tumor-level effect is positive in at least 4 of 5 tumors with the interval above zero. If the stress signature fails this, the stress score is not a coherent module in these data and the failed control says little about the pipeline.
3. Signed random-set statistics: a model has a material bias if the mean t of its matched random gene sets, averaged over arrays, differs from zero by more than 0.5. If so, z_vs_random (real t minus mean random t, divided by the random SD) is reported as a supplementary bias-corrected metric, labeled post-hoc, and is never substituted for the locked criterion.

## Results of the post-hoc diagnostics (step 05, run 2026-10-22; post-hoc, not part of the verdict)

Check: the delta = 0 row reproduces the step 04 H1 result exactly (3 of 5 tumors positive, mean 0.068, interval -0.108 to 0.245).

1. Spike-in (shift added to the lipid score of IGF2BP2-detected spots, in SD of the score). The locked criterion is first met at 0.4 SD (mean 0.458, interval 0.287 to 0.629, median p_vs_random 0.040). At 0.3 SD the interval excludes zero (0.190 to 0.535) but median p_vs_random is 0.114. At 0.2 SD the interval excludes zero (0.091 to 0.440) and median p_vs_random is 0.338. At 0.1 SD or less the interval includes zero. Approximate minimum detectable tumor-level effect at 80 percent power from the interval of the unshifted analysis (t test over 5 tumors, ignoring the random-set condition): 0.24 SD. This is conditional on the observed noise, one realization. The data rule out a tumor-level lipid shift much above 0.25 SD (upper bound of the interval) and cannot distinguish zero from smaller effects.

2. Split-half plumbing control (tumor-level effect of one random half of a signature on the other, with the covariates of the locked models): stromal 5 of 5 tumors positive, mean 0.388 (0.213 to 0.564), recovered. Macrophage 5 of 5, mean 0.175 (0.052 to 0.298), recovered. Lipid 4 of 5, mean 0.049 (-0.033 to 0.131), not recovered. Stress 4 of 5, mean 0.055 (-0.001 to 0.111), not recovered. Under the rule written beforehand, the stress score is not a coherent module in these data, so the failed control of step 04 says little about the pipeline. The pipeline recovers real gene-gene co-variation for the stromal and macrophage signatures. The lipid score also shows almost no internal co-variation and is a weak outcome measure. Rough reliability of the full score by the Spearman-Brown formula from the half-half association (itself adjusted, so only an indicator): lipid 0.09, stress 0.10, macrophage 0.30, stromal 0.56.

3. Signed random-set statistics. Mean t of matched random gene sets averaged over arrays: control 0.107, H1 0.391, H2 0.251. All are below the 0.5 threshold written beforehand, so by that rule no model has material bias on average. Array A3 is a clear exception: mean random t 2.34 (H1) and 1.80 (H2), with 100 and 97.5 percent of random sets positive. Random genes score higher in IGF2BP2-detected spots of A3, consistent with a nonlinear depth or detection effect that the linear covariate does not remove (not verified). In A3 the real H1 t (2.05) is below the random mean, so the A3 H1 effect is explained by this artifact. Supplementary bias-corrected metric z_vs_random (tumor level, post-hoc): H1 mean -0.012 (-0.694 to 0.670), 4 of 5 tumors positive; per-array values for the arrays with apparent positive H1 effects are A1 1.86, A3 -0.33 and B3 1.51. H2 mean -0.595 (-1.512 to 0.321), 1 of 5 positive. Control mean -0.638 (-1.036 to -0.241), 0 of 5 positive, interpreted with caution because the stress score is not coherent. The random-set t distribution of the control model is overdispersed (mean SD 1.37 against about 1.0 for H1 and H2), consistent with the anti-conservative block p values found in step 04.

Conclusion: the step 04 verdict stands (control failed, H1 and H2 inconclusive and unvalidated). The diagnostics show that the lipid outcome is a weak measure, that the control was built on an incoherent score, and that the apparent positive H1 effects shrink after bias correction (corrected tumor-level mean about zero). There is no evidence for a lipid association. The data rule out a large one and are uninformative about a small one.

## Software environment and deviations from the project summary (recorded 2026-10-22)

- Containers were not used. The pipeline runs in conda environments located by path in nextflow.config: the analysis environment (environment.yml, with exact versions in environment.lock.yml) and the cell2location environment (scripts/build_c2l_env.sbatch, versions in environment.cell2loc.lock.txt). The project summary planned containers.
- cell2location was installed and verified (imports work and a V100 GPU is visible) but was not used, because the single-cell reference is too thin for deconvolution (7 fibroblasts, no T or NK cells in tumor B). Cell-type composition is represented by marker-based scores.
- Data are fetched from GEO by scripts/download_data.sh. Nothing from GSE111672 is stored in this repository.

## External check on GSE202051 (plan locked 2026-10-22, before any IGF2BP2 association was computed in this dataset)

Data: GSE202051_totaldata-final-toshare.h5ad (Hwang et al. 2022), single-nucleus RNA-seq, 224,988 nuclei, 22,164 genes, 43 patients (pid) in 44 samples. Only the structure of the file has been inspected (columns, sizes, gene presence). No IGF2BP2 detection rate or association has been looked at.

Unit and cells: the unit is the patient (pid). Cancer nuclei are those with Level 1 Annotation equal to "Epithelial (malignant)" (64,538 nuclei, the same count as in the two other label columns). The labels were assigned by the dataset authors, apparently with copy-number inference (suggested by column names, not verified). Expression is the log-normalized X, detection means counts above zero in the counts layer, and depth is log of total_counts.
Primary population: patients with treatment_status "Untreated". Secondary, reported separately: all patients.
Eligibility: a patient is eligible if, among their cancer nuclei, at least 30 have IGF2BP2 counts and at least 30 do not. The eligibility step reads only how many nuclei have IGF2BP2 counts. If fewer than 8 untreated patients are eligible, the check is declared infeasible and no association is computed.

E1, reproducibility of an IGF2BP2-linked program across tumors (leave one tumor out). Per patient and gene: regress log-normalized expression on IGF2BP2 detection, log total_counts and a cell-cycle score (the same 19 genes as step 03), and keep the t statistic for detection. Gene universe: genes detected in at least 5 percent of the pooled cancer nuclei of the eligible patients, excluding IGF2BP2. For each eligible patient, the program is the 100 genes with the largest mean t over the other eligible patients. The held-out cancer nuclei are scored with the mean z-scored expression of those genes, the score is residualized on log depth and cycle score, and the AUC for IGF2BP2-detected versus undetected nuclei is computed. Statistic: mean held-out AUC over patients, with a 95 percent t interval over patients and the share of patients with AUC above 0.5. Null: IGF2BP2 detection shuffled among nuclei of similar depth (10 depth strata) within each patient, with the whole leave-one-out derivation repeated 200 times. E1 passes if the mean AUC is above 0.5, the one-sided empirical p is at most 0.05, and the AUC is above 0.5 in at least 60 percent of eligible patients.

E2, lipid score in cancer nuclei (cell-intrinsic version of H1). Per patient: the z-scored lipid signature (the 15 genes locked above) regressed on IGF2BP2 detection with log depth and cycle score as covariates, effect in SD of the score. Random-set comparison: 200 matched random gene sets per patient, as in step 04, giving p_vs_random per patient. E2 meets the criterion if the effect is positive in at least 70 percent of eligible patients, the 95 percent t interval over patients excludes zero, and the median p_vs_random over patients is below 0.05.

Both tests are reported whatever the outcome, and they do not change the verdicts of steps 03 and 04. A pass of E1 would show that an IGF2BP2-linked program reproduces in nuclei across tumors. It would not show that this is the State 3 of Chen et al. (defined by clustering). A failure would not refute our result in cells, because nuclei and cells differ. This dataset was also used by Chen et al. for their treatment-response comparison, so it is independent of our data but not of their paper. Anything beyond E1 and E2 is post-hoc.

## External check: eligibility result (step 07, run 2026-10-22)

The three cancer-label columns agree for all 224,988 nuclei (64,538 cancer nuclei). The counts layer is integer-valued (maximum 1,156 in a sample of 2,000 nuclei). 43 patients, none with more than one treatment status. Untreated: 18 patients, 52,999 cancer nuclei, all 18 eligible (threshold 8), so the check is feasible. Treated groups: CRT 14 patients (6 eligible), CRTl 5 (0), CRTx 2 (1), RT 1 (1), CRTln 1 (0), CRTn 1 (0), GART 1 (0). Only the number of cancer nuclei with IGF2BP2 counts was read. No association has been computed.
IGF2BP2 is detected in 28 to 90 percent of cancer nuclei per untreated patient (above 75 percent in 12 of 18), much more often than in the spatial spots. The primary contrast of detected versus undetected therefore compares most nuclei with a minority of low-depth or dropout nuclei.

Amendment A1 (written before any association was computed in this dataset): a secondary exposure is added to E1 and E2, namely IGF2BP2 log-normalized expression strictly above the patient's median among cancer nuclei (equal to detection when the median is zero). It goes through exactly the same procedures and is reported beside the primary result, never in place of it. If the two disagree, both are reported and the disagreement is stated.

## Reproducibility check (clean-room run, 2026-10-22)

A fresh clone of the repository at commit dc82ec8 was run in a new folder with data downloaded again by scripts/download_data.sh (archive size verified, 25 count files extracted), using one command: nextflow run main.nf -profile explorer --step all. Six tasks completed in 7 minutes 17 seconds. Ten key result tables are byte-identical to those of the original step-by-step runs: spatial_qc, sc_qc, clusters_A, clusters_B, reproducibility, array_tests, tumor_summary, spike_in_summary, split_half_summary and random_set_bias_summary. Limits: same cluster and software; the conda environments were reused by path and not rebuilt from the lock files; figures, h5ad files and the external steps 07 and 08 (run separately with sbatch) were not compared.