#!/usr/bin/env bash
# Download GEO GSE111672 (Moncada et al. 2020) and extract the count tables, without the H&E images.
# Usage: bash scripts/download_data.sh TARGET_DIR      (needs about 1 GB; checks free space first)
set -euo pipefail
target="${1:?Usage: download_data.sh TARGET_DIR}"
base="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111672/suppl"
bash "$(dirname "$0")/check_space.sh" "$target" 5
mkdir -p "$target/raw" "$target/counts"
cd "$target/raw"
for f in GSE111672_RAW.tar GSE111672_PDAC-A-indrop-filtered-expMat.txt.gz GSE111672_PDAC-B-indrop-filtered-expMat.txt.gz; do
  [ -s "$f" ] || curl -fO "$base/$f"
done
size=$(stat -c %s GSE111672_RAW.tar)
if [ "$size" != "749291520" ]; then
  echo "GSE111672_RAW.tar has $size bytes, expected 749291520" >&2
  exit 1
fi
tar -xf GSE111672_RAW.tar -C ../counts --wildcards '*.tsv.gz' '*filtered.txt.gz'
echo "count files extracted: $(ls ../counts | wc -l) (expected 25)"