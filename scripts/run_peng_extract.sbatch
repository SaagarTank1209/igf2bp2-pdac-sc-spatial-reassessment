#!/usr/bin/env bash
#SBATCH --job-name=peng_extract
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/%u/pdac_spatial/logs/peng_extract_%j.out
# Extracts the cancer cells of the Peng et al. dataset (Zenodo copy 10.5281/zenodo.3969339, raw.h5ad) and checks eligibility.
# Submit from the repository root:  sbatch scripts/run_peng_extract.sbatch
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
D=/scratch/$USER/pdac_spatial/external/Peng2019
env PYTHONNOUSERSITE=1 /scratch/$USER/envs/pdac-spatial/bin/python notebooks/15_peng_extract.py \
    --h5ad "$D/raw.h5ad" --outdir /scratch/$USER/pdac_spatial/results/15_peng --cache_dir "$D"
