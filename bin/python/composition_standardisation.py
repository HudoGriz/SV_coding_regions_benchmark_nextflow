#!/usr/bin/env python3
"""Does the SV mix explain the gap between the simulated sets and the target?

The simulated interval sets retain more deletions and longer events than the
real target. If that composition difference is what separates their metrics
from the target's, reweighting each simulated set to the target's own mix
should close the gap. This script tests that directly, from the per-stratum
counts that metric_decomposition.py writes.

For each simulated set, the metric is recomputed as a weighted mean of its
per-stratum values, with weights taken from the real target:

  recall     truth strata (scored type x size bin), weights = target truth mix
  precision  candidate strata, weights = target candidate mix

Two stratifications are reported: type only (DEL vs INS, the deletion-
composition argument) and type x size. A stratum a simulated set does not
populate is dropped from that set and the remaining weights renormalised; the
share of target weight that could be used is reported.

Output: one row per assembly x pipeline x metric x stratification, with the
target value, the simulated median raw and standardised, both gaps, and the
target's percentile in each distribution.
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np

METRICS = {"recall": "truth", "precision": "candidate"}


def percentile(values: np.ndarray, value: float) -> float:
    return float((values < value).mean() * 100 + (values == value).mean() * 50)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--strata", action="append", required=True, type=Path,
                        help="<prefix>.strata.tsv from metric_decomposition.py; repeatable")
    parser.add_argument("--target", default="wes_utr")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    # counts[(assembly, pipeline, target_set, target, side)][(svtype, size_bin)] = [n, tp]
    counts = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    for path in args.strata:
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                key = (row["assembly"], row["pipeline"], row["target_set"], row["target"], row["side"])
                cell = counts[key][(row["svtype"], row["size_bin"])]
                cell[0] += int(row["n"])
                cell[1] += int(row["target_tp"])

    def collapse(strata: dict, scheme: str) -> dict:
        if scheme == "type_size":
            return strata
        out = defaultdict(lambda: [0, 0])
        for (svtype, _), (n, tp) in strata.items():
            out[(svtype,)][0] += n
            out[(svtype,)][1] += tp
        return out

    rows = []
    pipelines = sorted({(key[0], key[1]) for key in counts})
    for assembly, pipeline in pipelines:
        for metric, side in METRICS.items():
            target_strata = counts.get((assembly, pipeline, "real", args.target, side))
            if not target_strata:
                continue
            simulations = [key for key in counts
                           if key[:3] == (assembly, pipeline, "simulated") and key[4] == side]
            if not simulations:
                continue
            for scheme in ("type", "type_size"):
                target = collapse(target_strata, scheme)
                total = sum(n for n, _ in target.values())
                weights = {stratum: n / total for stratum, (n, _) in target.items() if n}
                target_value = sum(tp for _, tp in target.values()) / total
                raw, standardised, weight_used = [], [], []
                for key in simulations:
                    simulated = collapse(counts[key], scheme)
                    n_all = sum(n for n, _ in simulated.values())
                    raw.append(sum(tp for _, tp in simulated.values()) / n_all if n_all else np.nan)
                    usable = {s: w for s, w in weights.items() if simulated.get(s, [0, 0])[0] > 0}
                    used = sum(usable.values())
                    weight_used.append(used)
                    standardised.append(
                        sum(w * simulated[s][1] / simulated[s][0] for s, w in usable.items()) / used
                        if used else np.nan
                    )
                raw, standardised = np.array(raw), np.array(standardised)
                raw, standardised = raw[~np.isnan(raw)], standardised[~np.isnan(standardised)]
                rows.append({
                    "assembly": assembly, "pipeline": pipeline, "metric": metric,
                    "stratification": scheme, "target": args.target,
                    "target_value": target_value,
                    "simulated_median_raw": float(np.median(raw)),
                    "simulated_median_standardised": float(np.median(standardised)),
                    "gap_raw": float(np.median(raw)) - target_value,
                    "gap_standardised": float(np.median(standardised)) - target_value,
                    "target_percentile_raw": percentile(raw, target_value),
                    "target_percentile_standardised": percentile(standardised, target_value),
                    "median_target_weight_used": float(np.median(weight_used)),
                    "n_simulated": len(raw),
                })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows")


if __name__ == "__main__":
    main()
