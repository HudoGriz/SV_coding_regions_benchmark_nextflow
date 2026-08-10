#!/usr/bin/env python3
"""Type and size composition of the truth variants each target scores.

The manuscript claims that a target built from short dispersed intervals does
not sample the truth set uniformly: because a deletion occupies its full
reference span while an insertion occupies one position, short intervals
intercept deletions and long variants preferentially. That claim needs a
composition table for HCI, the two real diagnostic targets, and the 500
simulated interval sets, on both assemblies.

The scored truth set for a benchmark is exactly Truvari's tp-base plus fn, so
this reads those files rather than re-deriving membership from the BEDs. No
filtering logic is duplicated here: whatever Truvari scored is what is counted.

The base count is invariant across pipelines within a target (the truth set
does not depend on which caller it is compared against), so one pipeline is
read per target and the invariance is asserted against the metrics table.

Usage:
  target_composition.py --run-dir <clean-rerun dir> --outdir <evidence dir>
"""
from __future__ import annotations

import argparse
import gzip
import statistics as st
from pathlib import Path

TARGETS = {"high_confidence": "HCI", "gene_panel": "GP", "wes_utr": "EX+UTR"}
LONG_BP = 1000


def parse_info(info: str) -> dict[str, str]:
    out = {}
    for field in info.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            out[k] = v
        else:
            out[field] = ""
    return out


def read_truth(path: Path) -> list[tuple[str, int]]:
    """Return (svtype, abs(svlen)) for every record in a Truvari VCF."""
    recs = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = parse_info(f[7])
            svtype = info.get("SVTYPE", "")
            svlen = info.get("SVLEN", "")
            if svlen in ("", "."):
                # fall back to REF/ALT lengths when SVLEN is absent
                length = abs(len(f[4]) - len(f[3]))
            else:
                length = abs(int(float(svlen.split(",")[0])))
            if not svtype:
                svtype = "DEL" if len(f[3]) > len(f[4]) else "INS"
            recs.append((svtype, length))
    return recs


def scored(bench_dir: Path, prefix: str) -> list[tuple[str, int]]:
    tp = bench_dir / f"{prefix}.tp-base.vcf.gz"
    fn = bench_dir / f"{prefix}.fn.vcf.gz"
    return read_truth(tp) + read_truth(fn)


def summarize(recs: list[tuple[str, int]]) -> dict[str, float]:
    n = len(recs)
    if n == 0:
        return {"n": 0, "pct_del": float("nan"), "pct_ins": float("nan"),
                "pct_ge_1kb": float("nan"), "median_len": float("nan")}
    dels = sum(1 for t, _ in recs if t == "DEL")
    ins = sum(1 for t, _ in recs if t == "INS")
    long_ = sum(1 for _, l in recs if l >= LONG_BP)
    return {
        "n": n,
        "pct_del": 100.0 * dels / n,
        "pct_ins": 100.0 * ins / n,
        "pct_ge_1kb": 100.0 * long_ / n,
        "median_len": st.median(l for _, l in recs),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True, type=Path,
                    help="clean-rerun root holding results-GRCh37 / results-GRCh38")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--pipeline", default="ONT/Sniffles",
                    help="tech/tool whose Truvari output supplies the truth set")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    tech, tool = args.pipeline.split("/")

    rows = []
    for asm in ("GRCh37", "GRCh38"):
        root = args.run_dir / f"results-{asm}"
        if not root.is_dir():
            print(f"skip {asm}: {root} absent")
            continue

        for target, label in TARGETS.items():
            d = root / "real_intervals" / tech / "truvari" / tool / target
            if not d.is_dir():
                print(f"skip {asm} {label}: {d} absent")
                continue
            s = summarize(scored(d, f"{tech}-{tool}-{target}"))
            rows.append(dict(assembly=asm, target_set="real", target=label,
                             replicate="observed", **s))
            print(f"{asm} {label:7s} n={s['n']:6d} DEL={s['pct_del']:5.1f}% "
                  f">=1kb={s['pct_ge_1kb']:5.1f}%")

        simroot = root / "simulations" / "benchmarks" / tech / tool
        sims = sorted(p for p in simroot.glob(f"{tech}-{tool}-simulation*")
                      if p.is_dir()) if simroot.is_dir() else []
        per_sim = []
        for p in sims:
            s = summarize(scored(p.parent, p.name))
            per_sim.append(s)
            rows.append(dict(assembly=asm, target_set="simulated", target="simulated",
                             replicate=p.name.split("-")[-1], **s))
        if per_sim:
            def med(k):
                return st.median(x[k] for x in per_sim)
            print(f"{asm} simulated ({len(per_sim)} sets) median n={med('n'):.1f} "
                  f"DEL={med('pct_del'):.1f}% >=1kb={med('pct_ge_1kb'):.1f}%")
            rows.append(dict(assembly=asm, target_set="simulated", target="simulated",
                             replicate="median_of_sets", n=med("n"),
                             pct_del=med("pct_del"), pct_ins=med("pct_ins"),
                             pct_ge_1kb=med("pct_ge_1kb"), median_len=med("median_len")))
            pooled = dict(assembly=asm, target_set="simulated", target="simulated",
                          replicate="pooled_over_sets")
            tot = sum(x["n"] for x in per_sim)
            pooled.update(n=tot,
                          pct_del=sum(x["pct_del"] * x["n"] for x in per_sim) / tot,
                          pct_ins=sum(x["pct_ins"] * x["n"] for x in per_sim) / tot,
                          pct_ge_1kb=sum(x["pct_ge_1kb"] * x["n"] for x in per_sim) / tot,
                          median_len=med("median_len"))
            rows.append(pooled)

    out = args.outdir / "target_composition.tsv"
    cols = ["assembly", "target_set", "target", "replicate", "n",
            "pct_del", "pct_ins", "pct_ge_1kb", "median_len"]
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(
                f"{r[c]:.6f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    print(f"\nwrote {out} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
