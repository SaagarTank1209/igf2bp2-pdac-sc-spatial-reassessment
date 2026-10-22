#!/usr/bin/env bash
#SBATCH --job-name=build_c2l
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/%u/pdac_spatial/logs/build_c2l_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export PYTHONNOUSERSITE=1
export PIP_CACHE_DIR=/scratch/$USER/pip_cache
mkdir -p "$PIP_CACHE_DIR"
ENV_DIR=/scratch/$USER/envs/cell2loc
bash scripts/check_space.sh /scratch/$USER/envs 20
if [ -d "$ENV_DIR" ]; then echo "Environment already exists at $ENV_DIR"; exit 1; fi
conda create -y -p "$ENV_DIR" --override-channels -c conda-forge python=3.10
"$ENV_DIR/bin/python" -m pip install "cell2location[tutorials]"
"$ENV_DIR/bin/python" -m pip freeze > /scratch/$USER/pdac_spatial/cell2loc.lock.txt
"$ENV_DIR/bin/python" -c "import importlib.metadata as m; [print(p, m.version(p)) for p in ('cell2location','scvi-tools','torch','numpy','scanpy')]"
du -sh "$ENV_DIR"
echo "BUILD DONE"