#!/usr/bin/env python3
"""Stratify a broad HCI benchmark by target after matching (GA4GH style).

Matching happens once, in HCI. Truth records (TP-base, FN) are then assigned to
a target by their own location and candidate records (TP-comp, FP) by theirs,
and precision, recall and F1 are recomputed from the records inside. Nothing is
re-matched, so no candidate can lose its HCI partner to the target boundary.

Membership uses Truvari's own region code, exactly as `truvari bench
--includebed BED --bench-overlaps N` selects records: the BED is read into a
region tree with build_region_tree, overlapping regions are merged, and records
are taken with VariantFile.fetch_regions(overlap=N); N = 0 means containment.
The result can therefore be compared record for record with an independently
restricted benchmark of the same target.

Runs inside the Truvari image, which provides the truvari package.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import truvari


def output_vcf(bench_dir: Path, name: str) -> Path:
    matches = sorted(bench_dir.glob(f"*.{name}.vcf.gz")) or sorted(bench_dir.glob(f"{name}.vcf.gz"))
    if len(matches) != 1:
        raise RuntimeError(f"expected one {name}.vcf.gz in {bench_dir}, found {len(matches)}")
    return matches[0]


def count_inside(vcf_path: Path, bed: Path, overlap: int) -> int:
    vcf = truvari.VariantFile(str(vcf_path))
    tree = truvari.build_region_tree(vcf, None, str(bed))
    truvari.merge_region_tree_overlaps(tree)
    return sum(1 for _ in vcf.fetch_regions(tree, overlap=overlap))


def metrics(tp_base: int, fn: int, tp_comp: int, fp: int) -> tuple[float, float, float]:
    recall = tp_base / (tp_base + fn) if tp_base + fn else float("nan")
    precision = tp_comp / (tp_comp + fp) if tp_comp + fp else float("nan")
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else float("nan")
    return precision, recall, f1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hci-bench", required=True, type=Path, help="HCI benchmark directory")
    parser.add_argument("--target", action="append", required=True, help="NAME=BED; repeatable")
    parser.add_argument("--overlap", type=int, default=1, help="--bench-overlaps value (default 1)")
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--pipeline", required=True, help="label, e.g. 'ONT Sniffles'")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    rows = []
    for spec in args.target:
        name, bed = spec.split("=", 1)
        counts = {
            key: count_inside(output_vcf(args.hci_bench, key), Path(bed), args.overlap)
            for key in ("tp-base", "fn", "tp-comp", "fp")
        }
        precision, recall, f1 = metrics(counts["tp-base"], counts["fn"], counts["tp-comp"], counts["fp"])
        rows.append({
            "assembly": args.assembly, "pipeline": args.pipeline, "target": name,
            "method": "post_matching_stratification", "overlap": args.overlap,
            "TP-base": counts["tp-base"], "FN": counts["fn"],
            "TP-comp": counts["tp-comp"], "FP": counts["fp"],
            "truth_denominator": counts["tp-base"] + counts["fn"],
            "precision": precision, "recall": recall, "f1": f1,
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
