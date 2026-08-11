#!/usr/bin/env bash
# Regenerate the containment-versus-overlap interval-membership sensitivity.
#
# The primary benchmark scores a truth or candidate record whenever it
# intersects the target (--bench-overlaps). Stock Truvari instead requires the
# record to be fully contained in an --includebed interval. Supplementary
# Table 9 reports what changes between those two rules, and this script is the
# evidence behind it.
#
# Every threshold is read back from the params.json that the primary run wrote
# next to each benchmark, so the containment rerun differs from the published
# overlap run in exactly one flag: --bench-overlaps is omitted. Nothing is
# re-specified by hand.
#
# Submit, do not run on the login node:
#   sbatch bin/run_containment_overlap_sensitivity.sh <run-dir> <assembly> <outdir>
#SBATCH --job-name=truvari-containment
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=%x-%j.log
set -euo pipefail

run_dir=${1:?usage: $0 <clean-rerun dir> <assembly> <outdir>}
assembly=${2:?}
outdir=${3:?}
sif=${TRUVARI_SIF:-/home/45483vrhovsek/sif/truvari-5.4.0-bench-overlaps-numeric.sif}

results="$run_dir/results-$assembly"
work="$outdir/containment_runs"
mkdir -p "$work"

echo "run_dir=$run_dir assembly=$assembly"
echo "container=$sif"
sha256sum "$sif"

# One containment rerun per published overlap benchmark.
find "$results/real_intervals" -mindepth 6 -maxdepth 6 -name params.json -print0 |
while IFS= read -r -d '' pj; do
    bench_dir=$(dirname "$pj")            # .../<target>/<prefix>
    prefix=$(basename "$bench_dir")
    out="$work/$prefix"
    if [[ -s "$out/summary.json" ]]; then
        echo "skip $prefix (already done)"
        continue
    fi
    rm -rf "$out"
    echo "[$(date --iso-8601=seconds)] containment: $prefix"
    cmd=$(python3 - "$pj" "$out" <<'PY'
import json, shlex, sys, os
p = json.load(open(sys.argv[1]))
sif = os.environ.get("TRUVARI_SIF", "/home/45483vrhovsek/sif/truvari-5.4.0-bench-overlaps-numeric.sif")
a = ["singularity", "exec", "-B", "/home", "--pwd", "/tmp", sif, "truvari", "bench",
     "-b", p["base"], "-c", p["comp"], "-f", p["reference"],
     "--includebed", p["includebed"], "-o", sys.argv[2],
     "-r", str(p["refdist"]), "-p", str(p["pctseq"]), "-P", str(p["pctsize"]),
     "-O", str(p["pctovl"]), "-s", str(p["sizemin"]), "-S", str(p["sizefilt"]),
     "--sizemax", str(p["sizemax"]), "--bnddist", str(p["bnddist"]),
     "--pick", p["pick"], "--bSample", p["bSample"], "--cSample", p["cSample"],
     "--chunksize", str(p["chunksize"])]
if p.get("passonly"):
    a.append("--passonly")
if p.get("dup_to_ins"):
    a.append("--dup-to-ins")
# --bench-overlaps deliberately omitted: that is the whole experiment.
print(" ".join(shlex.quote(x) for x in a))
PY
)
    TRUVARI_SIF="$sif" bash -c "$cmd" > "$work/$prefix.log" 2>&1 || {
        echo "FAILED $prefix"; tail -20 "$work/$prefix.log"; exit 1; }
done

# Pair each containment rerun with its published overlap counterpart.
python3 - "$results" "$work" "$outdir/containment_vs_overlap_metrics.tsv" <<'PY'
import json, sys
from pathlib import Path

results, work, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
cols = ["assembly", "technology", "caller", "target", "mode",
        "TP_base", "FP", "FN", "precision", "recall", "f1", "truth_n", "called_n"]
asm = results.name.replace("results-", "")


def row(summary, tech, caller, target, mode):
    d = json.load(open(summary))
    return dict(assembly=asm, technology=tech, caller=caller, target=target, mode=mode,
                TP_base=d["TP-base"], FP=d["FP"], FN=d["FN"],
                precision=d["precision"], recall=d["recall"], f1=d["f1"],
                truth_n=d["base cnt"], called_n=d["comp cnt"])


rows = []
for pj in sorted(results.glob("real_intervals/*/truvari/*/*/*/params.json")):
    bench = pj.parent
    prefix = bench.name
    tech, caller, target = prefix.split("-", 2)
    ov = bench.parent / f"{prefix}.summary.json"
    if not ov.exists():
        ov = bench / "summary.json"
    rows.append(row(ov, tech, caller, target, "overlap"))
    cont = work / prefix / "summary.json"
    if cont.exists():
        rows.append(row(cont, tech, caller, target, "containment"))

with open(out, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")
print(f"wrote {out} ({len(rows)} rows)")
PY

echo "[$(date --iso-8601=seconds)] done"
