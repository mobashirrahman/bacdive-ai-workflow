#!/usr/bin/env bash
# Stage the offline fixture into a run directory.
#
# Usage: benchmark/.test/stage.sh RUN_DIR
#
# Files are touched in dependency order (genomes, proteins, annotations) so
# every staged output is newer than its input. Otherwise Snakemake may decide a
# staged annotation is stale and run the real InterProScan on the dummy proteins.
set -euo pipefail

run_dir=${1:?usage: stage.sh RUN_DIR}
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

mkdir -p "$run_dir"
cp -r "$here"/staged/* "$run_dir"/
for stage in genomes proteins interpro; do
    find "$run_dir/results/$stage" -type f -exec touch {} +
    sleep 1
done
