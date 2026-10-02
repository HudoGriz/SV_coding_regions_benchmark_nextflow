#!/usr/bin/env python3
"""Compare two runs' Truvari metrics, row for row.

Both inputs are long tables from collect_benchmark_metrics.py. Rows are joined
on setting, target set, target and pipeline. Counts must agree exactly for a
row to be identical; metrics are compared to 1e-12.

Two uses in the revision:
  reproducibility  fresh GRCh37 run vs the published run: every shared row is
                   expected to be identical
  old vs new       GRCh38 from the old strictly filtered BAMs vs the rebuilt
                   BAMs: the differences are the result

Outputs <prefix>.rows.tsv (every joined row with old, new and delta) and
<prefix>.summary.tsv (per setting x target set x pipeline: rows compared,
identical, different, only in one run).
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

KEY = ("setting", "target_set", "target", "pipeline")
COUNTS = ("TP-base", "FN", "TP-comp", "FP")
METRICS = ("precision", "recall", "f1")


def read(path: Path) -> dict:
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    table = {tuple(row[k] for k in KEY): row for row in rows}
    if len(table) != len(rows):
        raise SystemExit(f"duplicate keys in {path}")
    return table


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--old", required=True, type=Path)
    parser.add_argument("--new", required=True, type=Path)
    parser.add_argument("--prefix", required=True, type=Path)
    args = parser.parse_args()

    old, new = read(args.old), read(args.new)
    rows, summary = [], Counter()
    for key in sorted(set(old) | set(new)):
        before, after = old.get(key), new.get(key)
        group = key[:2] + (key[3],)
        if before is None or after is None:
            summary[group + ("only_new" if before is None else "only_old",)] += 1
            continue
        row = dict(zip(KEY, key))
        identical = True
        for count in COUNTS:
            a, b = int(before[count]), int(after[count])
            row[f"{count}_old"], row[f"{count}_new"], row[f"{count}_delta"] = a, b, b - a
            identical &= a == b
        for metric in METRICS:
            a, b = float(before[metric]), float(after[metric])
            row[f"{metric}_old"], row[f"{metric}_new"] = a, b
            row[f"{metric}_delta_pp"] = 100 * (b - a)
            identical &= abs(b - a) <= 1e-12
        row["identical"] = identical
        rows.append(row)
        summary[group + ("identical" if identical else "different",)] += 1

    args.prefix.parent.mkdir(parents=True, exist_ok=True)
    if rows:
        with open(f"{args.prefix}.rows.tsv", "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    groups = sorted({key[:3] for key in summary})
    with open(f"{args.prefix}.summary.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["setting", "target_set", "pipeline", "identical", "different", "only_old", "only_new"])
        for group in groups:
            writer.writerow([*group] + [summary[group + (s,)] for s in ("identical", "different", "only_old", "only_new")])
    total = Counter()
    for key, n in summary.items():
        total[key[3]] += n
    print(dict(total))


if __name__ == "__main__":
    main()
