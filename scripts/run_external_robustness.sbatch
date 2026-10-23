#!/usr/bin/env bash
#SBATCH --job-name=ext_robust
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/%u/pdac_spatial/logs/ext_robust_%j.out
# Post-hoc robustness and specificity checks of E1 on the cancer nuclei of GSE202051.
# Submit from the repository root:  sbatch scripts/run_external_robustness.sbatch
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
env PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4 /scratch/$USER/envs/pdac-spatial/bin/python notebooks/09_external_robustness.py \
    --cache /scratch/$USER/pdac_spatial/external/GSE202051/GSE202051_cancer_nuclei.h5ad \
    --outdir /scratch/$USER/pdac_spatial/results/09_external