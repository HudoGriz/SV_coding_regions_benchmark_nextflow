#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 ANALYSIS_SIF TRANSITIONS_TSV OUTPUT_DIR" >&2
    exit 2
fi

analysis_sif=$1
transitions=$2
output_dir=$3

apptainer exec "$analysis_sif" plot-target-boundary-mechanisms \
    --transitions "$transitions" \
    --output-dir "$output_dir"

