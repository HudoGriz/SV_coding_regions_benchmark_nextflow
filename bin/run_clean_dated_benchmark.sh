#!/usr/bin/env bash
set -euo pipefail

run_date=${1:-$(date +%F)}
# Locate the repository. Under sbatch this runs from a spool copy, so
# BASH_SOURCE does not point into the repo; SV_REPO_ROOT and the submit
# directory cover that case.
repo_root=""
for _candidate in "${SV_REPO_ROOT:-}" \
                  "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)" \
                  "${SLURM_SUBMIT_DIR:-}" "$PWD"; do
    if [[ -n "$_candidate" && -f "$_candidate/bin/common.sh" ]]; then
        repo_root=$_candidate; break
    fi
done
[[ -n "$repo_root" ]] || { echo "ERROR: cannot locate the repository; set SV_REPO_ROOT" >&2; exit 1; }
source "$repo_root/bin/common.sh"

# Everything site-specific is an environment variable. SV_DATA_ROOT is the only
# one without a sensible default, since it is where the prepared data lives.
data_root=${SV_DATA_ROOT:?set SV_DATA_ROOT to the directory holding the prepared per-assembly data}
rerun_root="$data_root/clean-rerun-$run_date"

# An extra Nextflow config for the local cluster. Optional: unset means the
# pipeline runs with its own configuration only.
hpc_config=${SV_HPC_CONFIG:-}

# Empty means "use the published defaults in nextflow.config", which are
# registry URIs Nextflow pulls itself. Set either to a local .sif to override.
analysis_sif=${ANALYSIS_SIF:-}
truvari_sif=${TRUVARI_SIF:-}

if [[ -e "$rerun_root" ]]; then
    echo "Refusing to overwrite existing rerun: $rerun_root" >&2
    exit 1
fi
[[ -d "$data_root" ]] || { echo "SV_DATA_ROOT is not a directory: $data_root" >&2; exit 1; }
for required in "$hpc_config" "$analysis_sif" "$truvari_sif"; do
    [[ -z "$required" || -e "$required" ]] || {
        echo "Missing required file: $required" >&2; exit 1; }
done

mkdir -p "$rerun_root/logs" "$rerun_root/work"
{
    echo "run_date=$run_date"
    echo "started_at=$(date --iso-8601=seconds)"
    echo "pipeline=$repo_root"
    echo "analysis_container=${analysis_sif:-<pipeline default>}"
    echo "truvari_container=${truvari_sif:-<pipeline default>}"
    git -C "$repo_root" rev-parse HEAD 2>/dev/null | sed 's/^/git_head=/' || true
    git -C "$repo_root" status --porcelain=v1 | sha256sum | sed 's/  -$/  git_status_porcelain/' || true
    git -C "$repo_root" diff --binary | sha256sum | sed 's/  -$/  git_tracked_diff/' || true
    # Checksum only the overrides. When a container comes from the pipeline
    # default it is a registry URI pinned to an immutable tag, which identifies
    # the image on its own.
    for image in "$analysis_sif" "$truvari_sif"; do
        [[ -n "$image" ]] && image_checksum "$image"
    done
} > "$rerun_root/RUN_MANIFEST.txt"

git -C "$repo_root" status --porcelain=v1 > "$rerun_root/SOURCE_STATUS.txt" || true
git -C "$repo_root" diff --binary > "$rerun_root/SOURCE_TRACKED.diff" || true

# Bring Nextflow onto PATH. Both steps are site-specific and skipped when the
# variables are unset, which is the right behaviour anywhere nextflow is already
# installed. SV_ENV_MODULE loads an environment module, SV_CONDA_ENV activates a
# conda environment.
if [[ -n "${SV_ENV_MODULE:-}" ]]; then
    module load "$SV_ENV_MODULE"
fi
if [[ -n "${SV_CONDA_ENV:-}" ]]; then
    set +u
    conda activate "$SV_CONDA_ENV"
    set -u
fi
command -v nextflow >/dev/null 2>&1 || {
    echo "ERROR: nextflow is not on PATH. Install it, or set SV_ENV_MODULE / SV_CONDA_ENV." >&2
    exit 1
}

overall_status=0
for assembly in GRCh37 GRCh38; do
    params_file="$data_root/$assembly/params_$assembly.yaml"
    assembly_work="$rerun_root/work/$assembly"
    assembly_out="$rerun_root/results-$assembly"
    log_file="$rerun_root/logs/nextflow.$assembly.log"
    stdout_file="$rerun_root/logs/run_$assembly.out"
    mkdir -p "$assembly_work"

    # Only pass the optional switches that are actually configured, so the
    # pipeline defaults apply everywhere else.
    nf_args=(-params-file "$params_file" -profile "${SV_PROFILE:-singularity}")
    [[ -n "$hpc_config" ]] && nf_args+=(-c "$hpc_config")
    [[ -n "$analysis_sif" ]] && nf_args+=(--analysis_container "$analysis_sif")
    [[ -n "$truvari_sif" ]] && nf_args+=(--truvari_container "$truvari_sif")

    echo "[$(date --iso-8601=seconds)] Starting clean $assembly run" | tee -a "$rerun_root/RUN_STATUS.log"
    if nextflow -log "$log_file" run "$repo_root" \
        "${nf_args[@]}" \
        -work-dir "$assembly_work" \
        --outdir "$assembly_out" \
        --run_name "${assembly}_clean_${run_date}" \
        --reference_assembly "$assembly" \
        --num_simulations 500 \
        --simulate_targets true \
        --gather_statistics true \
        --generate_transition_evidence true \
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
