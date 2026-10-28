# IGF2BP2 in pancreatic cancer: a spatially aware reassessment

BINF6420 final project, Northeastern University, Fall 2026. Author: Sagar Shantilal Tank.

## Summary

Chen et al. (2026) report that pancreatic cancer cells with high IGF2BP2 form a distinct state that is linked to lipid metabolism and immune suppression. This project tests part of that claim with public data, using rules written down before the results, statistics that count patients (not cells or spots) as the units, and negative controls.

Findings in short:

1. **Spatial test (GSE111672, 9 tissue arrays from 5 tumors): inconclusive.** IGF2BP2 detected spots show no clear link to a lipid score (mean 0.068 SD, interval -0.108 to 0.245, 3 of 5 tumors positive) or a macrophage score (-0.098 SD, interval -0.237 to 0.041). A positive control failed, and the data could only have detected effects of about 0.4 SD or larger.
2. **Reproducibility in two tumors: not shown.** Genes linked to IGF2BP2 in tumor A and in tumor B overlapped by 1 gene against 2.6 expected by chance.
3. **External check (GSE202051, 43 patients, 64,538 cancer nuclei): a reproducible association, but not specific to IGF2BP2.** Across 18 untreated patients a program derived from the other patients separated IGF2BP2 detected from undetected nuclei (mean AUC 0.566), but 49 of 100 random genes detected at similar rates did as well or better. No lipid link was found (mean effect -0.042 SD, interval -0.101 to 0.017).
4. **State extension in two datasets: the group repeats, but IGF2BP2 is not special.** Cancer cells were clustered after patient alignment and tested on held out patients in 20 random splits, first in the nuclei of GSE202051 (18 untreated patients) and then in the source dataset of the paper (Peng et al.; whole cells, 23 untreated tumors, the original authors' cancer cell labels). Both gave outcome B: a cancer cell cluster with more IGF2BP2 recurs across patients, but IGF2BP2 did not stand out from its 100 matched comparison genes. It beat 95 of them in enrichment in 0 of 20 splits in both datasets, and in signature in 0 of 20 splits (nuclei) and 3 of 20 splits (whole cells); 14 were needed.
5. **Overall:** no evidence here for an IGF2BP2 specific malignant state with lipid or immune links. This does not refute the paper: its state gene list and code are not public, and every test here uses IGF2BP2 detection or our own clustering as a stand in.

## What was and was not tested

Tested: IGF2BP2 detection against lipid, macrophage, and exhaustion-related scores in tissue maps; reproducibility of an IGF2BP2-linked program across tumors; the lipid link inside cancer nuclei.
Not tested: the paper's State 3 as defined by clustering, T-cell exhaustion in space (the genes are almost undetectable), survival or immunotherapy prediction, the SGMS2 and PD-L1 mechanism, and the drug (montelukast) experiments.

## Data

Public data only. No data are stored in this repository.

- GSE111672 (Moncada et al. 2020): single-cell data for tumors A (1,926 cells) and B (1,733 cells), and 9 spatial arrays from 5 tumors. Exclusions and cleaning rules are in `docs/input_manifest.md`.
- GSE202051 (Hwang et al. 2022): single-nucleus data, 224,988 nuclei from 43 patients, with the authors' cell labels and treatment status.

## Pipeline

Nextflow on Slurm (Northeastern Explorer). Steps are numbered scripts in `notebooks/`.

| Step | Script | What it does |
|---|---|---|
| 01 | `01_load_qc.py` | Loads and cleans the data, writes QC tables |
| 02 | `02_sc_annotation.py` | Clusters single cells and writes marker tables for labeling |
| 03 | `03_state_program.py` | Reproducibility gate: IGF2BP2 linked genes in one tumor tested in the other |
| 04 | `04_spatial_tests.py` | Spatial tests with covariates, a block shuffle null, and matched random gene sets |
| 05 | `05_posthoc_diagnostics.py` | Post hoc: detectability, signature coherence, random set bias |
| 06 | `06_figures.py` | Figures 1 to 5 |
| 07 | `07_external_extract.py` | Extracts cancer nuclei from GSE202051 and checks eligibility |
| 08 | `08_external_analysis.py` | External tests E1 (program reproducibility) and E2 (lipid link) |
| 09 | `09_external_robustness.py` | Post hoc: technical covariates and matched random exposures |
| 10 | `10_external_figure.py` | Figure 6 |
| 11 | `11_state_map.py` | Comparison genes and the patient aligned map of cancer cells |
| 12 | `12_state_test.py` | Held out test of the IGF2BP2 enriched cluster (criteria C0 to C4) |
| 13 | `13_state_figure.py` | Figures 7 and 8 |
| 14 | `14_leave_one_out.py` | Post hoc: spatial result without array A3 and with each tumor dropped |
| 15 | `15_peng_extract.py` | Extracts the cancer cells of the source dataset (Peng et al.) and checks eligibility |

Steps 01 to 06 run in one command. Steps 07 to 09, 11 to 12 and 15 run through the separate `scripts/run_*.sbatch` files; steps 10, 13 and 14 are plain Python scripts that read the saved tables.

## Reproduce

1. Clone the repository on Explorer, for example into `/scratch/$USER`.
2. Build the environment (`scripts/build_env.sbatch`). The conda environment is located by path in `nextflow.config`. Exact versions are in `environment.lock.yml`; `environment.extension.txt` lists harmonypy 2.1.0, installed afterwards with `pip install --no-deps`. The cell2location environment (`scripts/build_c2l_env.sbatch`) was built and tested but is not used by the analysis.
3. Download the data: `bash scripts/download_data.sh /scratch/$USER/pdac_spatial`
4. Run the main pipeline: `sbatch scripts/run_pipeline.sbatch --step all --datadir /scratch/$USER/pdac_spatial --outdir /scratch/$USER/pdac_spatial/results`
5. External check, in order: `sbatch scripts/run_external_extract.sbatch`, then `scripts/run_external_analysis.sbatch`, then `scripts/run_external_robustness.sbatch`.
6. State extension in the nuclei: a smoke test first with `NSPLITS=3 sbatch scripts/run_state_extension.sbatch`, then the full run on the map already built with `SKIP_MAP=1 sbatch scripts/run_state_extension.sbatch`. The secondary run uses `POP=all_eligible sbatch scripts/run_state_extension.sbatch`.
7. Source dataset of the paper: download `StdWf1_PRJCA001063_CRC_besca2.raw.h5ad` from the Zenodo record 10.5281/zenodo.3969339 to `/scratch/$USER/pdac_spatial/external/Peng2019/raw.h5ad` (md5 78c75a92b45724e20700ae76f51819f1), run `sbatch scripts/run_peng_extract.sbatch`, then the state extension with the cache and output folder set: `CACHE=/scratch/$USER/pdac_spatial/external/Peng2019/Peng2019_cancer_cells.h5ad OUT=/scratch/$USER/pdac_spatial/results/12_peng NSPLITS=3 sbatch scripts/run_state_extension.sbatch`, and the full run with `SKIP_MAP=1` and without `NSPLITS`.

A clean room run of the main pipeline from a fresh clone (commit dc82ec8) finished six tasks in 7 minutes 17 seconds, and ten key result tables were byte identical to the original runs. Steps 07 to 15 were not part of that check. Paths and Slurm settings are specific to Explorer. Small result tables and figures are saved in `docs/results` and `docs/figures`.

## Repository layout

- `docs/input_manifest.md`: the project record. It holds the locked inputs, each plan written before its analysis, and every result with its verdict.
- `docs/results`, `docs/figures`: result tables and figures.
- `notebooks/`, `modules/local/`, `main.nf`, `nextflow.config`: analysis scripts and workflow.
- `scripts/`: environment builds, data download, space check, Slurm job files.
- `config/cluster_labels.yaml`: provisional cell-type labels for the single-cell clusters (marker-based, no copy-number evidence).
- `.claude/`: project rules for the AI assistant (see below).

## Pre-specification and provenance

The rules for each analysis were committed before its results existed (all dates are in 2026). The external check plan is commit 6f419e3 (22 October). Amendment A1, a secondary exposure, is commit dc82ec8, written after detection rates were seen and before any association was computed. The analysis script and implementation clarifications are commit 804a97b, and the robustness plan and script 09 are commit 1ae8e72. The state extension plan is commit c38d789 (22 October, 23:36). Its text was accidentally overwritten in commit e9a0f9b (28 October) and restored unchanged from c38d789 in commit 3f2b2c0, before any state analysis was run; the clarifications and scripts are in e9a0f9b, and a 3 split smoke test preceded the full run. The plan for the second dataset (Peng et al.) is commit daef1c8 (28 October); the file structure had been inspected and the average IGF2BP2 expression by cell type had been seen on the TISCH2 portal, but no association was computed, and the extraction result and the rule for a short comparison list were recorded before the map was built. Steps 05 and 09, the leave one out tables and the bias corrected values are labeled post hoc in `docs/input_manifest.md`.

Two design points are recorded openly. The external pass rule for E1 tested for a reproducible association and not for specificity to IGF2BP2, so E1 passed as written, and the matched random exposures of step 09 were added because of that gap. The state extensions then required specificity (beating 95 of 100 matched genes) in their rules from the start.

## Limitations

- Five tumors in the spatial analysis, and old arrays with 20 to 70 cells per spot, shallow depth, and almost no detectable immune genes.
- IGF2BP2 detection is a noisy stand in for the paper's state, whose gene list is not public. The state extensions used our own clustering, so a negative result is not a failure to reproduce the paper.
- The copy of the source dataset (Roche reprocessing of the data of Peng et al.) is not necessarily what Chen et al. processed. Its cancer cell labels are the original authors' and rest on copy number variation.
- The state extensions depend on the clustering resolution and the maps were built on all patients, including those of the validation half. In the whole cells the group repeats in 15 of 20 splits (14 needed) and the chosen cluster changes with the discovery half.
- The comparison genes are matched on detection and expression level but not on gene length. Whether the shared axis is technical or biological is not determined.
- Cancer cell calls in GSE111672 rest on marker genes only. The GSE202051 labels are the authors' and may contain errors, and the stricter label sensitivity was not run.
- The positive control failed because the stress score is not a coherent module in these data.
- The lipid and macrophage gene lists were chosen for this project and are not the paper's scores.
- GSE202051 was also used by Chen et al. for other analyses.

## AI assistance

Claude (Anthropic) was used as a coding and planning assistant. It wrote the scripts, tested their logic on synthetic data in a sandbox without access to the real data, and drafted documentation. I ran all analyses on the cluster and reviewed the outputs. No model runs inside the pipeline. The `.claude` folder holds project rules for the assistant.

## References

- Chen et al. (2026), npj Precision Oncology. 
- Moncada R, et al. Integrating microarray-based spatial transcriptomics and single-cell RNA-seq reveals tissue architecture in pancreatic ductal adenocarcinomas. Nat Biotechnol (2020).
- Hwang WL, et al. Single-nucleus and spatial transcriptome profiling of pancreatic cancer identifies multicellular dynamics associated with neoadjuvant treatment. Nat Genet (2022).

## License

MIT.