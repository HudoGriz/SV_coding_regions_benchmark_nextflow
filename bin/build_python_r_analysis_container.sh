#!/usr/bin/env bash
set -euo pipefail

if [[ $# -gt 1 ]]; then
    echo "Usage: $0 [OUTPUT_SIF]" >&2
    exit 2
fi

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_sif=${1:-"$repo_root/containers/python-r-analysis.sif"}

mkdir -p "$(dirname "$output_sif")"
cd "$repo_root"
apptainer build --force "$output_sif" containers/Singularity.python-r-analysis
apptainer test "$output_sif"
echo "Built and tested: $output_sif"

