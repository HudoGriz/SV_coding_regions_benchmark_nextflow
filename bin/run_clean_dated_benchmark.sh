#!/usr/bin/env bash
set -euo pipefail

run_date=${1:-$(date +%F)}
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
data_root=/home/45483vrhovsek/genome_in_bottle/SV_coding_regions_benchmark_nextflow_run_with_prep
rerun_root="$data_root/clean-rerun-$run_date"
hpc_config=/home/Software/configs/nextflow_local/configs/conf/kisld_hpc.config
analysis_sif="$repo_root/containers/python-r-analysis.sif"
truvari_sif=/home/45483vrhovsek/sif/truvari-5.4.0-bench-overlaps-numeric.sif

if [[ -e "$rerun_root" ]]; then
    echo "Refusing to overwrite existing rerun: $rerun_root" >&2
    exit 1
fi
for required in "$analysis_sif" "$truvari_sif" "$hpc_config"; do
    [[ -e "$required" ]] || { echo "Missing required file: $required" >&2; exit 1; }
done

mkdir -p "$rerun_root/logs" "$rerun_root/work"
{
    echo "run_date=$run_date"
    echo "started_at=$(date --iso-8601=seconds)"
    echo "pipeline=$repo_root"
    echo "analysis_container=$analysis_sif"
    echo "truvari_container=$truvari_sif"
    git -C "$repo_root" rev-parse HEAD 2>/dev/null | sed 's/^/git_head=/' || true
    git -C "$repo_root" status --porcelain=v1 | sha256sum | sed 's/  -$/  git_status_porcelain/' || true
    git -C "$repo_root" diff --binary | sha256sum | sed 's/  -$/  git_tracked_diff/' || true
    sha256sum "$analysis_sif" "$truvari_sif"
} > "$rerun_root/RUN_MANIFEST.txt"

git -C "$repo_root" status --porcelain=v1 > "$rerun_root/SOURCE_STATUS.txt" || true
git -C "$repo_root" diff --binary > "$rerun_root/SOURCE_TRACKED.diff" || true

module load anaconda
set +u
conda activate nf-core
set -u

overall_status=0
for assembly in GRCh37 GRCh38; do
    params_file="$data_root/$assembly/params_$assembly.yaml"
    assembly_work="$rerun_root/work/$assembly"
    assembly_out="$rerun_root/results-$assembly"
    log_file="$rerun_root/logs/nextflow.$assembly.log"
    stdout_file="$rerun_root/logs/run_$assembly.out"
    mkdir -p "$assembly_work"

    echo "[$(date --iso-8601=seconds)] Starting clean $assembly run" | tee -a "$rerun_root/RUN_STATUS.log"
    if nextflow -log "$log_file" run "$repo_root" \
        -params-file "$params_file" \
        -c "$hpc_config" \
        -profile cpu \
        -work-dir "$assembly_work" \
        --outdir "$assembly_out" \
        --run_name "${assembly}_clean_${run_date}" \
        --reference_assembly "$assembly" \
        --num_simulations 500 \
        --simulate_targets true \
        --gather_statistics true \
        --generate_transition_evidence true \
        --analysis_container "$analysis_sif" \
        --truvari_container "$truvari_sif" \
        -ansi-log false > "$stdout_file" 2>&1; then
        echo "[$(date --iso-8601=seconds)] Completed $assembly" | tee -a "$rerun_root/RUN_STATUS.log"
    else
        status=$?
        overall_status=$status
        echo "[$(date --iso-8601=seconds)] FAILED $assembly (exit $status)" | tee -a "$rerun_root/RUN_STATUS.log"
        tail -50 "$stdout_file" || true
    fi
done

echo "finished_at=$(date --iso-8601=seconds)" >> "$rerun_root/RUN_MANIFEST.txt"
exit "$overall_status"
