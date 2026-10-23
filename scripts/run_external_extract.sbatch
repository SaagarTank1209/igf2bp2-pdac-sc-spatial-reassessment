#!/usr/bin/env bash
#SBATCH --job-name=ext_extract
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/%u/pdac_spatial/logs/ext_extract_%j.out
# Downloads the GSE202051 single-nucleus file if needed, then runs notebooks/07_external_extract.py.
# Submit from the repository root:  sbatch scripts/run_external_extract.sbatch
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
EXT=/scratch/$USER/pdac_spatial/external/GSE202051
F=GSE202051_totaldata-final-toshare.h5ad
bash scripts/check_space.sh "$EXT" 15
if [ ! -s "$EXT/$F" ]; then
  [ -s "$EXT/$F.gz" ] || (cd "$EXT" && curl -fO "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202051/suppl/$F.gz")
  (cd "$EXT" && gunzip -k "$F.gz")
fi
env PYTHONNOUSERSITE=1 /scratch/$USER/envs/pdac-spatial/bin/python notebooks/07_external_extract.py \
    --h5ad "$EXT/$F" --outdir /scratch/$USER/pdac_spatial/results/07_external --cache_dir "$EXT"