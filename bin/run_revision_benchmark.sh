#!/usr/bin/env bash
#SBATCH --job-name=sv_revision
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00
#
# Run the full pipeline for one assembly into a new run directory:
#
#   bin/run_revision_benchmark.sh <run_label> <assembly>
#
# Results go to $SV_DATA_ROOT/<run_label>/results-<assembly>. Each assembly is
# launched from its own directory there, so the two assemblies can run at the
# same time and each keeps its own Nextflow cache for -resume.
#
# The script refuses to write into an existing results directory. Set
# SV_RESUME=1 to continue an interrupted run of the same label and assembly;
# nothing outside that run directory is ever written.
#
# Environment (all optional except SV_DATA_ROOT):
#   SV_DATA_ROOT     directory holding the prepared per-assembly data
#   SV_PARAMS_FILE   params file (default $SV_DATA_ROOT/<assembly>/params_<assembly>.yaml)
#   SV_HPC_CONFIG    extra Nextflow config for the local cluster
#   SV_PROFILE       Nextflow profile (default singularity)
#   SV_ENV_MODULE    environment module that provides conda or nextflow
#   SV_CONDA_ENV     conda environment with nextflow
#   SV_MAX_TIME      ceiling on any task's time request (default 120.h)
#   SV_RESUME        1 to resume this run instead of refusing
#   ANALYSIS_SIF, TRUVARI_SIF   local images overriding the published defaults
set -euo pipefail

run_label=${1:?usage: run_revision_benchmark.sh <run_label> <assembly>}
assembly=${2:?usage: run_revision_benchmark.sh <run_label> <assembly>}
[[ "$assembly" == GRCh37 || "$assembly" == GRCh38 ]] || {
    echo "ERROR: assembly must be GRCh37 or GRCh38, got $assembly" >&2; exit 1; }

# Under sbatch this runs from a spool copy, so BASH_SOURCE does not point into
# the repo; SV_REPO_ROOT and the submit directory cover that case.
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

data_root=${SV_DATA_ROOT:?set SV_DATA_ROOT to the directory holding the prepared per-assembly data}
run_root="$data_root/$run_label"
results="$run_root/results-$assembly"
launch_dir="$run_root/launch-$assembly"
work_dir="$run_root/work/$assembly"
params_file=${SV_PARAMS_FILE:-$data_root/$assembly/params_$assembly.yaml}
hpc_config=${SV_HPC_CONFIG:-}
analysis_sif=${ANALYSIS_SIF:-}
truvari_sif=${TRUVARI_SIF:-}
resume=${SV_RESUME:-0}

if [[ -e "$results" && "$resume" != 1 ]]; then
    echo "Refusing to overwrite existing results: $results (set SV_RESUME=1 to resume this run)" >&2
    exit 1
fi
for required in "$params_file" "$hpc_config" "$analysis_sif" "$truvari_sif"; do
    [[ -z "$required" || -e "$required" ]] || { echo "Missing required file: $required" >&2; exit 1; }
done

mkdir -p "$run_root/logs" "$launch_dir" "$work_dir"
manifest="$run_root/RUN_MANIFEST.$assembly.txt"
stamp=$(date +%Y%m%dT%H%M%S)

# Compute nodes may have no git, so the pipeline code is also identified by a
# hash of the files that define it. The per-file listing lets any commit be
# matched to the run later: hash the same paths in a checkout of that commit.
pipeline_files="$run_root/logs/pipeline_files.$assembly.$stamp.sha256"
(
    cd "$repo_root"
    find main.nf nextflow.config nextflow_schema.json modules.json conf modules workflows bin preparation \
        -type f ! -name '*.pyc' ! -path '*/__pycache__/*' -print0 | sort -z | xargs -0 sha256sum
) > "$pipeline_files"
{
    echo "run_label=$run_label"
    echo "assembly=$assembly"
    echo "started_at=$(date --iso-8601=seconds)"
    echo "resume=$resume"
    echo "slurm_job_id=${SLURM_JOB_ID:-}"
    echo "pipeline=$repo_root"
    echo "pipeline_tree_sha256=$(sha256sum < "$pipeline_files" | cut -d' ' -f1)  ($(wc -l < "$pipeline_files") files, listed in $pipeline_files)"
    if command -v git >/dev/null 2>&1; then
        git -C "$repo_root" rev-parse --abbrev-ref HEAD 2>/dev/null | sed 's/^/git_branch=/' || true
        git -C "$repo_root" rev-parse HEAD 2>/dev/null | sed 's/^/git_head=/' || true
        git -C "$repo_root" status --porcelain=v1 | sha256sum | sed 's/  -$/  git_status_porcelain/' || true
        git -C "$repo_root" diff --binary | sha256sum | sed 's/  -$/  git_tracked_diff/' || true
    else
        echo "git=unavailable on $(hostname); identify the commit with pipeline_tree_sha256"
    fi
    echo "params_file=$params_file"
    sha256sum "$params_file"
    echo "hpc_config=${hpc_config:-<none>}"
    [[ -n "$hpc_config" ]] && sha256sum "$hpc_config"
    echo "analysis_container=${analysis_sif:-<pipeline default>}"
    echo "truvari_container=${truvari_sif:-<pipeline default>}"
    for image in "$analysis_sif" "$truvari_sif"; do
        [[ -n "$image" ]] && image_checksum "$image"
    done
    echo "----"
} >> "$manifest"
# Snapshot the exact inputs the run was started with, one copy per start.
cp "$params_file" "$run_root/logs/params.$assembly.$stamp.yaml"
if command -v git >/dev/null 2>&1; then
    git -C "$repo_root" status --porcelain=v1 > "$run_root/logs/source_status.$assembly.$stamp.txt" || true
    git -C "$repo_root" diff --binary > "$run_root/logs/source_tracked.$assembly.$stamp.diff" || true
fi

# Bring Nextflow onto PATH; both steps are skipped when unset.
if [[ -n "${SV_ENV_MODULE:-}" ]]; then
    module load "$SV_ENV_MODULE"
fi
if [[ -n "${SV_CONDA_ENV:-}" ]]; then
    set +u
    eval "$(conda shell.bash hook)"
    conda activate "$SV_CONDA_ENV"
    set -u
fi
command -v nextflow >/dev/null 2>&1 || {
    echo "ERROR: nextflow is not on PATH. Install it, or set SV_ENV_MODULE / SV_CONDA_ENV." >&2
    exit 1
}
export NXF_OPTS=${NXF_OPTS:--Xms1g -Xmx6g}

nf_args=(-params-file "$params_file" -profile "${SV_PROFILE:-singularity}")
[[ -n "$hpc_config" ]] && nf_args+=(-c "$hpc_config")
[[ -n "$analysis_sif" ]] && nf_args+=(--analysis_container "$analysis_sif")
[[ -n "$truvari_sif" ]] && nf_args+=(--truvari_container "$truvari_sif")
[[ "$resume" == 1 ]] && nf_args+=(-resume)

cd "$launch_dir"
echo "[$(date --iso-8601=seconds)] Starting $assembly (resume=$resume)" | tee -a "$run_root/RUN_STATUS.log"
status=0
nextflow -log "$run_root/logs/nextflow.$assembly.$stamp.log" run "$repo_root" \
    "${nf_args[@]}" \
    -work-dir "$work_dir" \
    --outdir "$results" \
    --run_name "${assembly}_${run_label}" \
    --reference_assembly "$assembly" \
    --max_time "${SV_MAX_TIME:-120.h}" \
    --num_simulations 500 \
    --simulate_targets true \
    --gather_statistics true \
    --generate_transition_evidence true \
    --sensitivity_benchmarks true \
    -ansi-log false > "$run_root/logs/run_$assembly.$stamp.out" 2>&1 || status=$?

if [[ "$status" == 0 ]]; then
    echo "[$(date --iso-8601=seconds)] Completed $assembly" | tee -a "$run_root/RUN_STATUS.log"
else
    echo "[$(date --iso-8601=seconds)] FAILED $assembly (exit $status)" | tee -a "$run_root/RUN_STATUS.log"
    tail -50 "$run_root/logs/run_$assembly.$stamp.out" || true
fi
echo "finished_at=$(date --iso-8601=seconds) exit=$status" >> "$manifest"
exit "$status"
