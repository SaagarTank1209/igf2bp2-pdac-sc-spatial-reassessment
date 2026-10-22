#!/usr/bin/env bash
# CUDA 12.1 build for Explorer GPU drivers (CUDA 12.3). Replaces the earlier CUDA 13 build.
#SBATCH --job-name=build_c2l_cu121
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/%u/pdac_spatial/logs/build_c2l_cu121_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export PYTHONNOUSERSITE=1
export PIP_CACHE_DIR=/scratch/$USER/pip_cache
mkdir -p "$PIP_CACHE_DIR"
ENV_DIR=/scratch/$USER/envs/cell2loc_cu121
bash scripts/check_space.sh /scratch/$USER/envs 20
if [ -d "$ENV_DIR" ]; then echo "Environment already exists at $ENV_DIR"; exit 1; fi
conda create -y -p "$ENV_DIR" --override-channels -c conda-forge python=3.10
PY="$ENV_DIR/bin/python"
"$PY" -m pip install torch==2.5.1 --index-url https://download.pytorch.org/whl/cu121
"$PY" -m pip install "cell2location[tutorials]"
"$PY" -m pip freeze > /scratch/$USER/pdac_spatial/cell2loc_cu121.lock.txt
"$PY" -c "import importlib.metadata as m; [print(p, m.version(p)) for p in ('cell2location','scvi-tools','torch','numpy','scanpy')]"
"$PY" -c "import torch, cell2location, scvi; print('imports OK'); print('torch', torch.__version__, 'cuda build', torch.version.cuda)"
du -sh "$ENV_DIR"
echo "BUILD DONE"