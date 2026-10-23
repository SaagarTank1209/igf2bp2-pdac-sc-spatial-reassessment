# IGF2BP2 in pancreatic cancer: a spatially aware reassessment

BINF6420 final project, Northeastern University, Fall 2026. Author: Sagar Shantilal Tank.

## Summary

Chen et al. (2026) report that pancreatic cancer cells with high IGF2BP2 form a distinct state that is linked to lipid metabolism and immune suppression. This project tests part of that claim with public data. It uses rules written down before the results, statistics that count patients (not cells or spots) as the units, and negative controls.

Findings in short:

1. **Spatial test (GSE111672, 9 tissue arrays from 5 tumors): inconclusive.** IGF2BP2-detected spots show no clear link to a lipid score (mean 0.068 SD, interval -0.108 to 0.245, 3 of 5 tumors positive) or a macrophage score (-0.098 SD, interval -0.237 to 0.041). A positive control failed, and the data could only have detected effects of about 0.4 SD or larger.
2. **Reproducibility in two tumors: not shown.** Genes linked to IGF2BP2 in tumor A and in tumor B overlapped by 1 gene against 2.6 expected by chance.
3. **External check (GSE202051, 43 patients, 64,538 cancer nuclei): a reproducible association, but not specific to IGF2BP2.** Across 18 untreated patients a program derived from the other patients separated IGF2BP2-detected from undetected nuclei (mean AUC 0.566). But 49 of 100 randomly chosen genes detected at similar rates did as well or better, so this reflects a shared axis of variation and not an IGF2BP2-specific program. No lipid link was found (mean effect -0.042 SD, interval -0.101 to 0.017).
4. **Overall:** no evidence here for an IGF2BP2-specific malignant program with lipid or immune links. This does not refute the paper. Every test uses IGF2BP2 detection as a stand-in for a cluster-defined state whose gene list is not public, and the data are small or non-specific.

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
| 03 | `03_state_program.py` | Reproducibility gate: IGF2BP2-linked genes in one tumor tested in the other |
| 04 | `04_spatial_tests.py` | Spatial tests with covariates, a block-shuffle null, and matched random gene sets |
| 05 | `05_posthoc_diagnostics.py` | Post-hoc: detectability, signature coherence, random-set bias |
| 06 | `06_figures.py` | Figures |
| 07 | `07_external_extract.py` | Extracts cancer nuclei from GSE202051 and checks eligibility |
| 08 | `08_external_analysis.py` | External tests E1 (program reproducibility) and E2 (lipid link) |
| 09 | `09_external_robustness.py` | Post-hoc: technical covariates and matched random exposures |

Steps 01 to 06 run in one command. Steps 07 to 09 run through the separate `scripts/run_external_*.sbatch` files.

## Reproduce

1. Clone the repository on Explorer, for example into `/scratch/$USER`.
2. Build the environment (`scripts/build_env.sbatch`). The conda environment is located by path in `nextflow.config`. Exact versions are in `environment.lock.yml`. The cell2location environment (`scripts/build_c2l_env.sbatch`) was built and tested but is not used by the analysis.
3. Download the data: `bash scripts/download_data.sh /scratch/$USER/pdac_spatial`
4. Run the main pipeline:

```
sbatch scripts/run_pipeline.sbatch --step all --datadir /scratch/$USER/pdac_spatial --outdir /scratch/$USER/pdac_spatial/results
```

5. External check, in order: `sbatch scripts/run_external_extract.sbatch`, then `scripts/run_external_analysis.sbatch`, then `scripts/run_external_robustness.sbatch`.

A clean-room run of the main pipeline from a fresh clone (commit dc82ec8) finished six tasks in 7 minutes 17 seconds, and ten key result tables were byte-identical to the original runs. Paths and Slurm settings are specific to Explorer. Small result tables and figures are saved in `docs/results` and `docs/figures`.

## Repository layout

- `docs/input_manifest.md`: the project record. It holds the locked inputs, each plan written before its analysis, and every result with its verdict.
- `docs/results`, `docs/figures`: result tables and figures.
- `notebooks/`, `modules/local/`, `main.nf`, `nextflow.config`: analysis scripts and workflow.
- `scripts/`: environment builds, data download, space check, Slurm job files.
- `config/cluster_labels.yaml`: provisional cell-type labels for the single-cell clusters (marker-based, no copy-number evidence).
- `.claude/`: project rules for the AI assistant (see below).

## Pre-specification and provenance

The rules for each analysis were committed before its results existed (dates of 22 October 2026). The external check plan is commit 6f419e3. Amendment A1, a secondary exposure, is commit dc82ec8, written after detection rates were seen and before any association was computed. The analysis script and implementation clarifications are commit 804a97b. The robustness plan and script 09 are commit 1ae8e72. Steps 05 and 09 and the bias-corrected values are labeled post-hoc in `docs/input_manifest.md`.

One design flaw is recorded openly: the external pass rule for E1 tested for a reproducible association and not for specificity to IGF2BP2, so E1 passed as written. The matched random exposures in step 09 were added because of that gap.

## Limitations

- Five tumors in the spatial analysis, and old arrays with 20 to 70 cells per spot, shallow depth, and almost no detectable immune genes.
- IGF2BP2 detection is a noisy stand-in for the paper's state. The paper's gene list is not in its supplement.
- The cancer-cell calls in GSE111672 rest on marker genes only. The GSE202051 labels come from the dataset authors and may contain errors.
- The positive control failed because the stress score is not a coherent module in these data.
- The lipid and macrophage gene lists were chosen for this project and are not the paper's scores.
- GSE202051 is nuclei, not whole cells, and the paper also used this dataset for other analyses.
- Matched random exposures were matched loosely, and gene length was not controlled.

## AI assistance

Claude (Anthropic) was used as a coding and planning assistant. It wrote the scripts, tested their logic on synthetic data in a sandbox without access to the real data, and drafted documentation. I ran all analyses on the cluster and reviewed the outputs. No model runs inside the pipeline. The `.claude` folder holds project rules for the assistant.

## References

- Chen et al. (2026), npj Precision Oncology. 
- Moncada R, et al. Integrating microarray-based spatial transcriptomics and single-cell RNA-seq reveals tissue architecture in pancreatic ductal adenocarcinomas. Nat Biotechnol (2020).
- Hwang WL, et al. Single-nucleus and spatial transcriptome profiling of pancreatic cancer identifies multicellular dynamics associated with neoadjuvant treatment. Nat Genet (2022).

## License

MIT.