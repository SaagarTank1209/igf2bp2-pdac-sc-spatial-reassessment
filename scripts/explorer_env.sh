#!/usr/bin/env bash
# Run at the start of every Explorer session:  source scripts/explorer_env.sh
module load nextflow/25.04.6
export NXF_HOME=/scratch/$USER/nxf_home
export NXF_APPTAINER_CACHEDIR=/scratch/$USER/apptainer_images
export APPTAINER_CACHEDIR=/scratch/$USER/apptainer_cache
export APPTAINER_TMPDIR=/scratch/$USER/apptainer_tmp
mkdir -p "$NXF_HOME" "$NXF_APPTAINER_CACHEDIR" "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
echo "Nextflow environment ready"