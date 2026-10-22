#!/bin/bash -l
#SBATCH --job-name=pdac_nf
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/%u/pdac_spatial/logs/nf_%j.out
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
source scripts/explorer_env.sh
bash scripts/check_space.sh /scratch/$USER/pdac_spatial 10
nextflow run main.nf -profile explorer "$@"