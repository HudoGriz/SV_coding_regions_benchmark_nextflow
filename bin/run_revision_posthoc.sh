#!/usr/bin/env bash
#SBATCH --job-name=sv_posthoc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#
# Post-hoc analyses of one finished revision run, for one assembly:
#
#   bin/run_revision_posthoc.sh <run_root> <assembly>
#
# Reads <run_root>/results-<assembly> and the params file named in the run
# manifest; writes <run_root>/posthoc/<assembly>/ and refuses to overwrite it.
# Every script runs from this repository inside a pinned image: Truvari's for
# the steps that use Truvari's own code, the analysis image for the rest.
#
# Environment (optional):
#   ANALYSIS_SIF, TRUVARI_SIF   local images instead of the published URIs
#   SV_IMAGE_CACHE              where published images are pulled to
#   SV_POSTHOC_JOBS             parallel per-pipeline jobs (default: cpus)
#   SV_COMPARE_RESULTS          results dir of an earlier run to compare against
#                               (GRCh37: reproducibility; GRCh38: old vs new BAMs)
set -euo pipefail

run_root=${1:?usage: run_revision_posthoc.sh <run_root> <assembly>}
assembly=${2:?usage: run_revision_posthoc.sh <run_root> <assembly>}
run_root=$(cd "$run_root" && pwd)

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

results="$run_root/results-$assembly"
out="$run_root/posthoc/$assembly"
[[ -d "$results" ]] || { echo "ERROR: no results at $results" >&2; exit 1; }
[[ ! -e "$out" ]] || { echo "Refusing to overwrite $out" >&2; exit 1; }

manifest="$run_root/RUN_MANIFEST.$assembly.txt"
params_file=$(awk -F= '/^params_file=/{print $2; exit}' "$manifest")
[[ -f "$params_file" ]] || { echo "ERROR: params file from $manifest not found: $params_file" >&2; exit 1; }
param() { awk -F"'" -v key="$1" '$0 ~ "^"key":" {print $2; exit}' "$params_file"; }
hci_bed=$(param high_confidence_targets)
gene_panel_bed=$(param gene_panel_targets)
target_bed=$(param wes_utr_targets)
truth_vcf=$(param benchmark_vcf)
reference=$(param fasta)
tandem_repeats=$(param tandem_repeats)
for required in "$hci_bed" "$gene_panel_bed" "$target_bed" "$truth_vcf" "$reference"; do
    [[ -f "$required" ]] || { echo "ERROR: missing $required" >&2; exit 1; }
done

mkdir -p "$out"/{decomposition,stratified,uncertainty,sensitivity,annotations,logs}
image_cache=${SV_IMAGE_CACHE:-$HOME/.cache/sv-benchmark-images}
analysis_sif=$(resolve_image "${ANALYSIS_SIF:-$ANALYSIS_IMAGE_DEFAULT}" "$image_cache")
truvari_sif=$(resolve_image "${TRUVARI_SIF:-$TRUVARI_IMAGE_DEFAULT}" "$image_cache")
engine=$(container_engine)
# One bind per distinct directory; the engine warns about repeats.
binds=()
declare -A bound=()
for path in "$repo_root" "$run_root" "$(dirname "$params_file")" \
            "$(dirname "$(readlink -f "$hci_bed")")" "$(dirname "$(readlink -f "$gene_panel_bed")")" \
            "$(dirname "$(readlink -f "$target_bed")")" "$(dirname "$(readlink -f "$truth_vcf")")" \
            "$(dirname "$(readlink -f "$reference")")" ${SV_COMPARE_RESULTS:+"$SV_COMPARE_RESULTS"}; do
    [[ -n "${bound[$path]:-}" ]] && continue
    bound[$path]=1
    binds+=(-B "$path")
done
analysis() { "$engine" exec "${binds[@]}" "$analysis_sif" python3 "$@"; }
truvari_py() { "$engine" exec "${binds[@]}" "$truvari_sif" python3 "$@"; }
py="$repo_root/bin/python"
jobs=${SV_POSTHOC_JOBS:-${SLURM_CPUS_PER_TASK:-4}}

# Annotation BEDs for the simulation-fidelity table. Pinned GIAB v3.3
# stratifications; the checksums go into the manifest.
giab=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3
curl -sSfL -o "$out/annotations/segdups.bed.gz" "$giab/${assembly}@all/SegmentalDuplications/${assembly}_segdups.bed.gz"
curl -sSfL -o "$out/annotations/lowmappability.bed.gz" "$giab/${assembly}@all/Mappability/${assembly}_lowmappabilityall.bed.gz"

pipeline_files="$out/logs/pipeline_files.sha256"
( cd "$repo_root" && find main.nf nextflow.config conf modules workflows bin preparation -type f \
      ! -name '*.pyc' ! -path '*/__pycache__/*' -print0 | sort -z | xargs -0 sha256sum ) > "$pipeline_files"
{
    echo "run_root=$run_root"
    echo "assembly=$assembly"
    echo "started_at=$(date --iso-8601=seconds)"
    echo "slurm_job_id=${SLURM_JOB_ID:-}"
    echo "params_file=$params_file"
    echo "code_tree_sha256=$(sha256sum < "$pipeline_files" | cut -d' ' -f1)"
    command -v git >/dev/null 2>&1 && git -C "$repo_root" rev-parse HEAD | sed 's/^/git_head=/'
    image_checksum "$analysis_sif"
    image_checksum "$truvari_sif"
    sha256sum "$out"/annotations/*.bed.gz
} > "$out/POSTHOC_MANIFEST.txt"

log() { echo "[$(date --iso-8601=seconds)] $*" | tee -a "$out/POSTHOC_MANIFEST.txt"; }

# Pipelines scored against the simulated sets (WES is not).
mapfile -t pipelines < <(
    find "$results/simulations/benchmarks" -mindepth 2 -maxdepth 2 -type d -printf '%P\n' | sort | tr '/' ':'
)
log "pipelines: ${pipelines[*]}"

log "metrics and Truvari parameters"
analysis "$py/collect_benchmark_metrics.py" --results "$results" --assembly "$assembly" \
    --include-simulations --output "$out/metrics.tsv" --params-output "$out/truvari_parameters.tsv"

if [[ -n "${SV_COMPARE_RESULTS:-}" ]]; then
    log "comparison with $SV_COMPARE_RESULTS"
    analysis "$py/collect_benchmark_metrics.py" --results "$SV_COMPARE_RESULTS" --assembly "$assembly" \
        --include-simulations --output "$out/compare_baseline_metrics.tsv"
    analysis "$py/compare_runs.py" --old "$out/compare_baseline_metrics.tsv" --new "$out/metrics.tsv" \
        --prefix "$out/compare" | tee -a "$out/POSTHOC_MANIFEST.txt"
fi

log "SV-type record accounting"
truvari_py "$py/svtype_accounting.py" --results "$results" --assembly "$assembly" \
    --hci-bed "$hci_bed" --truth-vcf "$truth_vcf" --output "$out/svtype_accounting.tsv" \
    | tee -a "$out/POSTHOC_MANIFEST.txt"

export -f analysis truvari_py log
export engine analysis_sif truvari_sif py results assembly out target_bed gene_panel_bed
export binds_string="${binds[*]}"

log "decomposition, post-matching stratification and bootstrap per pipeline"
per_pipeline() {
    local spec=$1 name=${1/:/_}
    read -r -a binds <<< "$binds_string"
    "$engine" exec "${binds[@]}" "$analysis_sif" python3 "$py/metric_decomposition.py" \
        --results "$results" --assembly "$assembly" --pipeline "$spec" \
        --target "wes_utr=$target_bed" --target "gene_panel=$gene_panel_bed" --simulations \
        --prefix "$out/decomposition/$name"
    local tech=${spec%%:*} caller=${spec##*:}
    "$engine" exec "${binds[@]}" "$truvari_sif" python3 "$py/stratify_post_matching.py" \
        --hci-bench "$results/real_intervals/$tech/truvari/$caller/high_confidence" \
        --target "wes_utr=$target_bed" --target "gene_panel=$gene_panel_bed" \
        --assembly "$assembly" --pipeline "$tech $caller" --output "$out/stratified/$name.tsv"
    "$engine" exec "${binds[@]}" "$analysis_sif" python3 "$py/bootstrap_metrics.py" \
        --results "$results" --assembly "$assembly" --pipeline "$spec" \
        --target wes_utr --target-bed "$target_bed" --output "$out/uncertainty/$name.tsv"
}
export -f per_pipeline
printf '%s\n' "${pipelines[@]}" | xargs -P "$jobs" -I{} bash -c 'per_pipeline "$@" > "$out/logs/$(echo "$1" | tr : _).log" 2>&1' _ {}

# Post-matching stratification and the decomposition's composition-only terms
# score the same truth records and the same HCI false positives, so those must
# agree exactly. They may differ in one way only: a 30-49 bp candidate that
# matched in HCI stays a TP after stratification, while an independently
# restricted benchmark drops it (neither TP nor FP) once its truth partner is
# outside the target. Stratified TP-comp may therefore exceed the decomposition's
# HCI-TP candidates, and never fall below them.
log "stratification against the composition-only values"
analysis - "$out" <<'EOF' | tee -a "$out/POSTHOC_MANIFEST.txt"
import csv, sys
from pathlib import Path
out = Path(sys.argv[1])
bad = 0
for strat in sorted((out / "stratified").glob("*.tsv")):
    dec = out / "decomposition" / f"{strat.stem}.decomposition.tsv"
    comp = {r["target"]: r for r in csv.DictReader(dec.open(), delimiter="\t")}
    for r in csv.DictReader(strat.open(), delimiter="\t"):
        d = comp[r["target"]]
        extra_tp = int(r["TP-comp"]) - int(d["candidate_hci_tp"])
        ok = (int(r["TP-base"]) == int(d["truth_hci_tp"])
              and int(r["truth_denominator"]) == int(d["n_truth"])
              and int(r["FP"]) == int(d["candidate_hci_fp"])
              and extra_tp >= 0)
        bad += not ok
        print(f"{r['pipeline']} {r['target']}: {'ok' if ok else 'MISMATCH'}"
              f" (match-only candidates kept as TP by stratification: {extra_tp})")
sys.exit(1 if bad else 0)
EOF

if [[ -d "$results/sensitivity" ]]; then
    log "sensitivity audits"
    analysis "$py/sensitivity_transitions.py" --results "$results" --assembly "$assembly" \
        --target-bed "$target_bed" --prefix "$out/sensitivity/$assembly"
else
    log "no sensitivity benchmarks in this run; skipped"
fi

log "composition standardisation"
strata_args=()
for file in "$out"/decomposition/*.strata.tsv; do strata_args+=(--strata "$file"); done
analysis "$py/composition_standardisation.py" "${strata_args[@]}" --output "$out/composition_standardisation.tsv"

log "simulation fidelity"
annotation_args=(--annotation "segdups=$out/annotations/segdups.bed.gz"
                 --annotation "lowmappability=$out/annotations/lowmappability.bed.gz")
if [[ -f "$tandem_repeats" ]]; then
    trf_dir=$(dirname "$(readlink -f "$tandem_repeats")")
    [[ -z "${bound[$trf_dir]:-}" ]] && binds+=(-B "$trf_dir")
    annotation_args+=(--annotation "tandem_repeats=$tandem_repeats")
fi
analysis "$py/simulation_fidelity.py" --target-bed "$target_bed" \
    --simulation-dir "$results/simulations/target_regions" --reference "$reference" \
    "${annotation_args[@]}" --prefix "$out/fidelity"

log "done"
echo "finished_at=$(date --iso-8601=seconds)" >> "$out/POSTHOC_MANIFEST.txt"
