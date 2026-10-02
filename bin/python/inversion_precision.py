#!/usr/bin/env python3
"""Precision with and without candidate inversions, every pipeline and real target.

Neither truth set contains inversions, so every candidate inversion a benchmark
scores is a false positive. Removing them gives the precision on the scored
deletions and insertions alone.

    inversion_precision.py --results RESULTS_DIR --metrics metrics.tsv --output inversion_precision.tsv

RESULTS_DIR holds the primary real-target benchmarks in the published layout,
real_intervals/<technology>/truvari/<caller>/<target>/<technology>-<caller>-<target>.fp.vcf.gz.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path

INVERSION = re.compile(r"(^|;)SVTYPE=INV(;|$)")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    with open(args.metrics) as handle:
        primary = {(r["pipeline"], r["target"]): r for r in csv.DictReader(handle, delimiter="\t")
                   if r["setting"] == "primary" and r["target_set"] == "real"}

    rows = []
    for (pipeline, target), r in sorted(primary.items()):
        technology, caller = pipeline.split(" ")
        fp_vcf = (args.results / "real_intervals" / technology / "truvari" / caller / target
                  / f"{technology}-{caller}-{target}.fp.vcf.gz")
        inversions = 0
        with gzip.open(fp_vcf, "rt") as handle:
            for line in handle:
                if line[0] != "#" and INVERSION.search(line.split("\t", 8)[7]):
                    inversions += 1
        tp, fp = int(r["TP-comp"]), int(r["FP"])
        precision = tp / (tp + fp) if tp + fp else float("nan")
        without = tp / (tp + fp - inversions) if tp + fp - inversions else float("nan")
        rows.append([pipeline, target, fp, inversions, round(precision, 4), round(without, 4),
                     round(100 * (without - precision), 2)])

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["pipeline", "target", "fp", "fp_inversions", "precision", "precision_without_inversions",
                         "difference_pp"])
        writer.writerows(rows)


if __name__ == "__main__":
    main()
