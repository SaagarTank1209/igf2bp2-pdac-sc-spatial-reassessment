#!/usr/bin/env bash
#SBATCH --job-name=state_ext
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/%u/pdac_spatial/logs/state_ext_%j.out
# State extension on the cancer nuclei of GSE202051: builds the patient-corrected map, then runs the held-out
# tests at Leiden resolutions 0.3 (primary), 0.2 and 0.5, and the three-cluster variant k3 (S4).
# Defaults: untreated patients, main variant, 20 splits.
#   Smoke test first (3 splits, writes files tagged _n3):   NSPLITS=3 sbatch scripts/run_state_extension.sbatch
#   Full run on the map already built (skips the map):      SKIP_MAP=1 sbatch scripts/run_state_extension.sbatch
#   Secondary populations and variants:  POP=all_eligible sbatch ...      VAR=S1 sbatch ...
#   Another dataset (cache and output folder overridden):  CACHE=/path/cache.h5ad OUT=/path/results sbatch ...
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
POP=${POP:-untreated}
VAR=${VAR:-main}
NSPLITS=${NSPLITS:-20}
PY="env PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4 /scratch/$USER/envs/pdac-spatial/bin/python"
CACHE=${CACHE:-/scratch/$USER/pdac_spatial/external/GSE202051/GSE202051_cancer_nuclei.h5ad}
OUT=${OUT:-/scratch/$USER/pdac_spatial/results/11_state}
if [ "${SKIP_MAP:-0}" = "1" ] && [ -s "$OUT/clusters_${POP}_${VAR}.tsv" ]; then
  echo "skipping the map: using $OUT/clusters_${POP}_${VAR}.tsv"
else
  $PY notebooks/11_state_map.py --cache "$CACHE" --outdir "$OUT" --population "$POP" --variant "$VAR"
fi
for R in r0.3 r0.2 r0.5 k3; do
  $PY notebooks/12_state_test.py --cache "$CACHE" --outdir "$OUT" --population "$POP" --variant "$VAR" --resolution "$R" --n_splits "$NSPLITS"
done
