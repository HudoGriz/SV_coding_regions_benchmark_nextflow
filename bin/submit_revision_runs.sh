#!/usr/bin/env bash
# Submit bin/run_revision_benchmark.sh for one or both assemblies as Slurm jobs.
#
#   bin/submit_revision_runs.sh <run_label> [GRCh37] [GRCh38]
#
# Each job first writes its params file with preparation/generate_params.sh, as
# $SV_DATA_ROOT/<assembly>/params_<assembly>.<run_label>.yaml, unless that file
# already exists. Generating it inside the job lets a GRCh38 job wait for the
# BAM rebuild (SV_DEPENDENCY_GRCh38=afterok:<jobid>:...) and still pick up the
# rebuilt BAMs, which generate_params.sh only finds once they exist.
#
# The job runs a frozen copy of run_revision_benchmark.sh, taken from the
# committed version (git HEAD) at submission and named after the commit. Bash
# reads a script from disk while it runs it, so a job pointed at the working
# copy would pick up any later edit mid-run, from its old byte offset; the copy
# is never edited, and its name records the code the run used.
#
# Environment: everything run_revision_benchmark.sh reads, plus
#   SV_PARTITION              Slurm partition for the Nextflow head job (default cpu)
#   SV_HEAD_TIME              head job time limit (default 7-00:00:00, the QOS ceiling).
#                             A shorter limit backfills sooner on a busy cluster; a run
#                             that hits it is continued with SV_RESUME=1.
#   SV_DEPENDENCY_<assembly>  optional sbatch --dependency value for that assembly
set -euo pipefail

run_label=${1:?usage: submit_revision_runs.sh <run_label> [GRCh37] [GRCh38]}
shift
assemblies=("$@")
[[ ${#assemblies[@]} -gt 0 ]] || assemblies=(GRCh37 GRCh38)

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
data_root=${SV_DATA_ROOT:?set SV_DATA_ROOT to the directory holding the prepared per-assembly data}
log_dir="$data_root/$run_label/logs"
mkdir -p "$log_dir"
commit=$(git -C "$repo_root" rev-parse --short=12 HEAD)

for assembly in "${assemblies[@]}"; do
    [[ "$assembly" == GRCh37 || "$assembly" == GRCh38 ]] || { echo "Unknown assembly: $assembly" >&2; exit 1; }
    params="$data_root/$assembly/params_$assembly.$run_label.yaml"
    dependency_var="SV_DEPENDENCY_$assembly"
    dependency=${!dependency_var:-}

    sbatch_args=(
        --parsable
        --job-name="sv_${run_label}_${assembly}"
        --partition="${SV_PARTITION:-cpu}"
        --cpus-per-task=2 --mem=8G --time="${SV_HEAD_TIME:-7-00:00:00}"
        --output="$log_dir/head_${assembly}.%j.out"
        --export=ALL,SV_PARAMS_FILE="$params",SV_REPO_ROOT="$repo_root"
    )
    [[ -n "$dependency" ]] && sbatch_args+=(--dependency="$dependency")

    driver="$log_dir/run_revision_benchmark.$assembly.$commit.$(date +%Y%m%dT%H%M%S).sh"
    git -C "$repo_root" show "HEAD:bin/run_revision_benchmark.sh" > "$driver"
    chmod 0555 "$driver"

    job_id=$(sbatch "${sbatch_args[@]}" --wrap "set -euo pipefail
if [[ ! -e '$params' ]]; then
    bash '$repo_root/preparation/generate_params.sh' --genome '$assembly' \
        --datadir '$data_root/$assembly/data' --outfile '$params'
fi
exec bash '$driver' '$run_label' '$assembly'")
    echo "$assembly job $job_id${dependency:+ (dependency $dependency)} driver $driver (commit $commit)" \
        | tee -a "$log_dir/submitted_jobs.txt"
done
