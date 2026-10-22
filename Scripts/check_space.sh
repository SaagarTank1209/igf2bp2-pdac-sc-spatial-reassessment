#!/usr/bin/env bash
# Usage: bash scripts/check_space.sh TARGET_DIR NEEDED_GB
# Stops with an error if the filesystem holding TARGET_DIR has less than NEEDED_GB free.
set -euo pipefail
target="${1:?Usage: check_space.sh TARGET_DIR NEEDED_GB}"
need_gb="${2:?Usage: check_space.sh TARGET_DIR NEEDED_GB}"
mkdir -p "$target"
avail_kb=$(df -Pk "$target" | awk 'NR==2 {print $4}')
avail_gb=$(( avail_kb / 1024 / 1024 ))
used=$(du -sh "$target" 2>/dev/null | cut -f1 || true)
echo "Target: $target"
echo "Free on filesystem: ${avail_gb} GB (shared, not your personal quota)"
echo "Already used in target: ${used:-unknown}"
echo "Needed: ${need_gb} GB"
if [ "$avail_gb" -lt "$need_gb" ]; then
  echo "FAIL: not enough free space" >&2
  exit 1
fi
echo "OK"