#!/usr/bin/env bash
#SBATCH --job-name=build_env
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/%u/pdac_spatial/logs/build_env_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
ENV_DIR=/scratch/$USER/envs/pdac-spatial
bash scripts/check_space.sh /scratch/$USER/envs 15
if [ -d "$ENV_DIR" ]; then echo "Environment already exists at $ENV_DIR"; exit 1; fi
conda env create -p "$ENV_DIR" -f environment.yml
conda list -p "$ENV_DIR" | wc -l
du -sh "$ENV_DIR"
echo "BUILD DONE"