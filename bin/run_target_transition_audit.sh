#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 || $# -gt 5 ]]; then
    echo "Usage: $0 RESULTS ASSEMBLY ANALYSIS_SIF OUTPUT_DIR [PIPELINE_LABEL]" >&2
    exit 2
fi

results=$1
assembly=$2
analysis_sif=$3
output_dir=$4
pipeline_filter=${5:-}
batch_size=20

mkdir -p "$output_dir/batches"

pipelines=(
    "Illumina_WGS:Manta:Manta"
    "ONT:CuteSV:ONT_CuteSV"
    "ONT:Sniffles:ONT_Sniffles"
    "PacBio:CuteSV:PacBio_CuteSV"
    "PacBio:Pbsv:PacBio_Pbsv"
)

for spec in "${pipelines[@]}"; do
    IFS=: read -r tech caller label <<< "$spec"
    if [[ -n "$pipeline_filter" && "$label" != "$pipeline_filter" ]]; then
        continue
    fi
    for start in $(seq 1 "$batch_size" 500); do
        prefix="$output_dir/batches/${assembly}_${label}_${start}"
        if [[ -s "${prefix}.metadata.json" ]]; then
            continue
        fi
        lock="${prefix}.lock"
        if ! mkdir "$lock" 2>/dev/null; then
            continue
        fi
        trap 'rmdir "$lock" 2>/dev/null || true' EXIT
        apptainer exec "$analysis_sif" target-transition-audit \
            --results "$results" \
            --assembly "$assembly" \
            --simulations \
            --simulation-start "$start" \
            --simulation-limit "$batch_size" \
            --pipeline "$tech:$caller" \
            --prefix "$prefix"
        rmdir "$lock"
        trap - EXIT
    done
done

if [[ -z "$pipeline_filter" ]]; then
    apptainer exec "$analysis_sif" merge-transition-audits \
        --input-dir "$output_dir/batches" \
        --prefix "$output_dir/${assembly}_simulations"
fi
