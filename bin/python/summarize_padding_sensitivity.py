#!/usr/bin/env python3
"""Summarize and plot Truvari target-padding sensitivity benchmarks."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


NAME_RE = re.compile(
    r"^(GRCh3[78])_(Illumina_WGS|ONT|PacBio)_(Manta|CuteSV|Sniffles|Pbsv)_pad(\d+)$"
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args()

    rows = []
    for directory in sorted(args.input_dir.iterdir()):
        if not directory.is_dir():
            continue
        match = NAME_RE.match(directory.name)
        if not match:
            continue
        summaries = list(directory.glob("summary.json"))
        if len(summaries) != 1:
            continue
        stats = json.loads(summaries[0].read_text())
        assembly, technology, caller, padding = match.groups()
        rows.append({
            "assembly": assembly,
            "pipeline": f"{technology.replace('_', ' ')} {caller}",
            "padding_bp": int(padding),
            "tp": stats["TP-base"], "fp": stats["FP"], "fn": stats["FN"],
            "precision": stats["precision"], "recall": stats["recall"],
            "f1": stats["f1"], "truth_n": stats["base cnt"],
            "called_n": stats["comp cnt"],
        })
    data = pd.DataFrame(rows)
    if data.empty:
        raise SystemExit("no completed padding summaries found")
    expected = 2 * 5 * 6
    if len(data) != expected:
        raise SystemExit(f"expected {expected} summaries, found {len(data)}")

    baseline = data[data.padding_bp == 0][
        ["assembly", "pipeline", "precision", "recall", "f1"]
    ].rename(columns={metric: f"{metric}_pad0" for metric in ("precision", "recall", "f1")})
    data = data.merge(baseline, on=["assembly", "pipeline"], validate="many_to_one")
    for metric in ("precision", "recall", "f1"):
        data[f"delta_{metric}_pp"] = 100 * (data[metric] - data[f"{metric}_pad0"])

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(f"{args.output_prefix}.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8), sharex=True, constrained_layout=True)
    metrics = ("precision", "recall", "f1")
    colors = plt.cm.tab10.colors
    for row_index, assembly in enumerate(("GRCh37", "GRCh38")):
        subset = data[data.assembly == assembly]
        for column_index, metric in enumerate(metrics):
            ax = axes[row_index, column_index]
            for color, (pipeline, group) in zip(colors, subset.groupby("pipeline", sort=True)):
                group = group.sort_values("padding_bp")
                ax.plot(group.padding_bp, group[metric], marker="o", linewidth=1.8,
                        markersize=4, label=pipeline, color=color)
            ax.set_title(f"{assembly} {metric.capitalize()}", loc="left", fontweight="bold")
            ax.set_xlabel("Target padding (bp)")
            ax.set_ylabel(metric.capitalize())
            ax.grid(alpha=0.2)
            if row_index == 0 and column_index == 2:
                ax.legend(frameon=False, fontsize=7, loc="best")
    fig.suptitle("Sensitivity of target-restricted SV metrics to boundary padding",
                 fontsize=14, fontweight="bold")
    fig.savefig(f"{args.output_prefix}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
