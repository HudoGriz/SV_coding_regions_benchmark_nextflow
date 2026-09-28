#!/usr/bin/env python3
"""Collect every Truvari summary of one pipeline run into a single long table.

Covers the primary real-target benchmarks, the breakend sensitivity, the
sensitivity benchmarks (thresholds, containment, --extend, padding) and,
optionally, the simulated interval sets. Each row carries the setting, target
and pipeline, the four counts and the three metrics exactly as Truvari wrote
them, so no metric is recomputed here.

Layouts read (as published by conf/modules.config):
  real_intervals/<tech>/truvari/<caller>/<target>/<tech>-<caller>-<target>.summary.json
  bnd_sensitivity/<mode>/real_intervals/<tech>/truvari/<caller>/<target>/...
  sensitivity/<setting>/<tech>/<caller>/<target>/<tech>-<caller>-<target>.summary.json
  simulations/benchmarks/<tech>/<caller>/<tech>-<caller>-simulation<N>.summary.json
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

FIELDS = [
    "assembly", "setting", "target_set", "target", "pipeline", "technology", "caller",
    "TP-base", "FN", "TP-comp", "FP", "truth_denominator", "precision", "recall", "f1",
    "summary_path",
]


def summary_row(path: Path, assembly: str, setting: str, target_set: str,
                technology: str, caller: str, target: str) -> dict:
    data = json.loads(path.read_text())
    return {
        "assembly": assembly, "setting": setting, "target_set": target_set,
        "target": target, "pipeline": f"{technology} {caller}",
        "technology": technology, "caller": caller,
        "TP-base": data["TP-base"], "FN": data["FN"],
        "TP-comp": data["TP-comp"], "FP": data["FP"],
        "truth_denominator": data["TP-base"] + data["FN"],
        "precision": data["precision"], "recall": data["recall"], "f1": data["f1"],
        "summary_path": str(path),
    }


def real_rows(root: Path, assembly: str, setting: str, target_set: str) -> list[dict]:
    rows = []
    for path in sorted(root.glob("*/truvari/*/*/*.summary.json")):
        target_dir = path.parent
        caller_dir = target_dir.parent
        technology = caller_dir.parent.parent.name
        rows.append(summary_row(path, assembly, setting, target_set,
                                technology, caller_dir.name, target_dir.name))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path, help="pipeline --outdir")
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--include-simulations", action="store_true")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    rows = real_rows(args.results / "real_intervals", args.assembly, "primary", "real")

    bnd_root = args.results / "bnd_sensitivity"
    if bnd_root.exists():
        for mode_dir in sorted(p for p in bnd_root.iterdir() if p.is_dir()):
            rows += real_rows(mode_dir / "real_intervals", args.assembly,
                              f"bnd_{mode_dir.name}", "real")

    sensitivity_root = args.results / "sensitivity"
    if sensitivity_root.exists():
        for path in sorted(sensitivity_root.glob("*/*/*/*/*.summary.json")):
            target_dir = path.parent
            caller_dir = target_dir.parent
            technology_dir = caller_dir.parent
            setting_dir = technology_dir.parent
            rows.append(summary_row(path, args.assembly, setting_dir.name, "real",
                                    technology_dir.name, caller_dir.name, target_dir.name))

    if args.include_simulations:
        for path in sorted((args.results / "simulations" / "benchmarks").glob("*/*/*.summary.json")):
            caller_dir = path.parent
            technology = caller_dir.parent.name
            simulation = path.name[len(f"{technology}-{caller_dir.name}-"):-len(".summary.json")]
            rows.append(summary_row(path, args.assembly, "primary", "simulated",
                                    technology, caller_dir.name, simulation))

    if not rows:
        raise SystemExit(f"no Truvari summaries found under {args.results}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows to {args.output}")


if __name__ == "__main__":
    main()
