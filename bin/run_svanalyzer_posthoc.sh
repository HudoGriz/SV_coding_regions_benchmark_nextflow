#!/usr/bin/env bash
#
# Second-comparator check, outside the pipeline:
#
#   bin/run_svanalyzer_posthoc.sh <run_root> <assembly> [TECH:CALLER ...]
#
# For each WGS pipeline of a finished run, the truth and candidate
# VCFs are pre-filtered independently to HCI and to the boundary target with
# Truvari's own membership and filters (prefilter_vcf.py), benchmarked with
# `svanalyzer benchmark`, and every HCI-TP to target-FN truth variant is traced
# to its HCI partners (svanalyzer_eligibility.py). Writes
# <run_root>/svanalyzer/<assembly>/ and refuses to overwrite it.
#
# Environment:
#   SVANALYZER_SIF   image built from containers/Singularity.svanalyzer (required)
#   TRUVARI_SIF, ANALYSIS_SIF   local images instead of the published URIs
#   SV_IMAGE_CACHE               where published images are pulled to
#   SV_POSTHOC_JOBS              parallel pipelines (default: all online CPUs)
#   SVA_MAXDIST (500), SVA_NORMSHIFT (1.0), SVA_NORMSIZEDIFF (0.3), SVA_NORMDIST (1.0)
#     maxdist mirrors Truvari's refdist, normsizediff 0.3 its pctsize 0.7, and
#     normdist 1.0 its pctseq 0 (no sequence requirement); normshift keeps the
#     SVanalyzer default.
set -euo pipefail

run_root=${1:?usage: run_svanalyzer_posthoc.sh <run_root> <assembly> [TECH:CALLER ...]}
assembly=${2:?usage: run_svanalyzer_posthoc.sh <run_root> <assembly> [TECH:CALLER ...]}
shift 2
run_root=$(cd "$run_root" && pwd)

repo_root=""
for _candidate in "${SV_REPO_ROOT:-}" \
                  "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)" \
                  "$PWD"; do
    if [[ -n "$_candidate" && -f "$_candidate/bin/common.sh" ]]; then
        repo_root=$_candidate; break
    fi
done
[[ -n "$repo_root" ]] || { echo "ERROR: cannot locate the repository; set SV_REPO_ROOT" >&2; exit 1; }
source "$repo_root/bin/common.sh"

results="$run_root/results-$assembly"
out="$run_root/svanalyzer/$assembly"
[[ -d "$results" ]] || { echo "ERROR: no results at $results" >&2; exit 1; }
[[ ! -e "$out" ]] || { echo "Refusing to overwrite $out" >&2; exit 1; }
svanalyzer_sif=${SVANALYZER_SIF:?set SVANALYZER_SIF to an image built from containers/Singularity.svanalyzer}
[[ -f "$svanalyzer_sif" ]] || { echo "ERROR: $svanalyzer_sif not found" >&2; exit 1; }

params_file=$(awk -F= '/^params_file=/{print $2; exit}' "$run_root/RUN_MANIFEST.$assembly.txt")
param() { awk -F"'" -v key="$1" '$0 ~ "^"key":" {print $2; exit}' "$params_file"; }
hci_bed=$(param high_confidence_targets)
target_bed=$(param wes_utr_targets)
truth_vcf=$(param benchmark_vcf)
reference=$(param fasta)

image_cache=${SV_IMAGE_CACHE:-$HOME/.cache/sv-benchmark-images}
analysis_sif=$(resolve_image "${ANALYSIS_SIF:-$ANALYSIS_IMAGE_DEFAULT}" "$image_cache")
truvari_sif=$(resolve_image "${TRUVARI_SIF:-$TRUVARI_IMAGE_DEFAULT}" "$image_cache")
engine=$(container_engine)

maxdist=${SVA_MAXDIST:-500}
normshift=${SVA_NORMSHIFT:-1.0}
normsizediff=${SVA_NORMSIZEDIFF:-0.3}
normdist=${SVA_NORMDIST:-1.0}

declare -A raw_vcf=(
    [Illumina_WGS:Manta]="Illumina_WGS/Manta/*.diploid_sv.vcf.gz"
    [Illumina_WGS:Delly]="Illumina_WGS/Delly/*.vcf.gz"
    [ONT:CuteSV]="ONT/CuteSV/*.vcf.gz"
    [ONT:Sniffles]="ONT/Sniffles/*.vcf.gz"
    [PacBio:CuteSV]="PacBio/CuteSV/*.vcf.gz"
    [PacBio:Pbsv]="PacBio/PBSV/*.vcf.gz"
)
# Runs with the 0/0 filter keep the VCF each benchmark actually scored in
# benchmarked_calls/<technology>/<caller>/; use it when it is there.
if [[ -d "$results/benchmarked_calls" ]]; then
    for spec in "${!raw_vcf[@]}"; do
        raw_vcf[$spec]="../benchmarked_calls/${spec%%:*}/${spec##*:}/*.benchmarked.vcf.gz"
    done
fi
pipelines=("$@")
if [[ ${#pipelines[@]} -eq 0 ]]; then
    for spec in "${!raw_vcf[@]}"; do
        compgen -G "$results/sv_calls/${raw_vcf[$spec]}" >/dev/null && pipelines+=("$spec")
    done
fi

mkdir -p "$out"/{inputs,runs,logs,reference}
# SVanalyzer (GTB::FASTA) reads <fasta>.fai and runs `samtools faidx` to rebuild
# it when it is missing or older than the FASTA. The FASTA is symlinked and its
# index copied, never linked: the copy is newer, so nothing is rebuilt, and a
# rebuild could never write through a link into the shared reference directory.
ln -s "$(readlink -f "$reference")" "$out/reference/reference.fasta"
cp "$(readlink -f "$reference").fai" "$out/reference/reference.fasta.fai"
{
    echo "run_root=$run_root"
    echo "assembly=$assembly"
    echo "started_at=$(date --iso-8601=seconds)"
    echo "params_file=$params_file"
    echo "settings: maxdist=$maxdist normshift=$normshift normsizediff=$normsizediff normdist=$normdist"
    echo "pipelines=${pipelines[*]}"
    image_checksum "$svanalyzer_sif"
    image_checksum "$truvari_sif"
    image_checksum "$analysis_sif"
} > "$out/SVANALYZER_MANIFEST.txt"

prefilter() {  # vcf bed output id-prefix
    "$engine" exec -B "$repo_root" -B "$run_root" "$truvari_sif" python3 "$repo_root/bin/python/prefilter_vcf.py" \
        --vcf "$1" --bed "$2" --overlap 1 --stable-ids --id-prefix "$4" --output "$3"
}
prefilter "$truth_vcf" "$hci_bed" "$out/inputs/truth.hci.vcf.gz" truth_
prefilter "$truth_vcf" "$target_bed" "$out/inputs/truth.target.vcf.gz" truth_

export engine repo_root run_root out results truvari_sif analysis_sif svanalyzer_sif assembly
export hci_bed target_bed maxdist normshift normsizediff normdist
export -f prefilter
run_pipeline() {
    local spec=$1 pattern=$2 name=${1/:/_}
    local comp
    comp=$(compgen -G "$results/sv_calls/$pattern" | head -1)
    prefilter "$comp" "$hci_bed" "$out/inputs/$name.hci.vcf.gz" cand_
    prefilter "$comp" "$target_bed" "$out/inputs/$name.target.vcf.gz" cand_
    for region in hci target; do
        mkdir -p "$out/runs/$name.$region"
        ( cd "$out/runs/$name.$region" && "$engine" exec -B "$run_root" "$svanalyzer_sif" \
            svanalyzer benchmark --ref "$out/reference/reference.fasta" \
                --test "$out/inputs/$name.$region.vcf.gz" --truth "$out/inputs/truth.$region.vcf.gz" \
                --maxdist "$maxdist" --normshift "$normshift" --normsizediff "$normsizediff" \
                --normdist "$normdist" --prefix "$out/runs/$name.$region/benchmark" )
    done
    "$engine" exec -B "$repo_root" -B "$run_root" "$analysis_sif" python3 "$repo_root/bin/python/svanalyzer_eligibility.py" \
        --assembly "$assembly" --pipeline "${spec/:/ }" \
        --hci-prefix "$out/runs/$name.hci/benchmark" --target-prefix "$out/runs/$name.target/benchmark" \
        --target-truth "$out/inputs/truth.target.vcf.gz" \
        --target-test "$out/inputs/$name.target.vcf.gz" \
        --normshift "$normshift" --normsizediff "$normsizediff" --normdist "$normdist" \
        --summary "$out/runs/$name.summary.tsv" --records "$out/runs/$name.transitions.tsv"
}
export -f run_pipeline
for spec in "${pipelines[@]}"; do
    printf '%s\t%s\n' "$spec" "${raw_vcf[$spec]}"
done | xargs -P "${SV_POSTHOC_JOBS:-$(getconf _NPROCESSORS_ONLN)}" -L1 bash -c \
    'run_pipeline "$0" "$1" > "$out/logs/${0/:/_}.log" 2>&1'

# One summary table for the assembly.
first=1
for file in "$out"/runs/*.summary.tsv; do
    if [[ $first == 1 ]]; then cat "$file"; first=0; else tail -n +2 "$file"; fi
done > "$out/svanalyzer_summary.tsv"
echo "finished_at=$(date --iso-8601=seconds)" >> "$out/SVANALYZER_MANIFEST.txt"
