#!/usr/bin/env bash
#SBATCH --job-name=ext_analysis
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/%u/pdac_spatial/logs/ext_analysis_%j.out
# Runs the locked E1 and E2 analyses on the cancer nuclei extracted by run_external_extract.sbatch.
# Submit from the repository root:  sbatch scripts/run_external_analysis.sbatch
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
env PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4 /scratch/$USER/envs/pdac-spatial/bin/python notebooks/08_external_analysis.py \
    --cache /scratch/$USER/pdac_spatial/external/GSE202051/GSE202051_cancer_nuclei.h5ad \
    --outdir /scratch/$USER/pdac_spatial/results/08_external