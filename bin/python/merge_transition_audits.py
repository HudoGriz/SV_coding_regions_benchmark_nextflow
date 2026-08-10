#!/usr/bin/env python3
"""Merge batched target-transition audit tables and regenerate summaries."""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


def read_rows(paths):
    rows = []
    for path in paths:
        with path.open() as handle:
            rows.extend(csv.DictReader(handle, delimiter="\t"))
    return rows


def write_rows(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--prefix", required=True, type=Path)
    args = parser.parse_args()

    transition_paths = sorted(args.input_dir.glob("*.transitions.tsv"))
    if not transition_paths:
        raise SystemExit("no transition tables found")
    rows = read_rows(transition_paths)
    fields = list(rows[0])
    write_rows(Path(f"{args.prefix}.transitions.tsv"), rows, fields)

    counts = Counter(
        (row["assembly"], row["target"], row["pipeline"], row["direction"], row["mechanism"])
        for row in rows
    )
    totals = Counter((row["assembly"], row["target"], row["pipeline"]) for row in rows)
    summary = []
    for key, count in sorted(counts.items()):
        summary.append({
            "assembly": key[0], "target": key[1], "pipeline": key[2],
            "direction": key[3], "mechanism": key[4], "n": count,
            "percent": 100 * count / sum(
                value for other, value in counts.items() if other[:4] == key[:4]
            ),
        })
    write_rows(
        Path(f"{args.prefix}.mechanisms.tsv"), summary,
        ["assembly", "target", "pipeline", "direction", "mechanism", "n", "percent"],
    )
    metadata = {
        "batch_count": len(transition_paths),
        "transition_count": len(rows),
        "targets": len({row["target"] for row in rows}),
        "pipelines": sorted({row["pipeline"] for row in rows}),
        "mechanism_counts": Counter(row["mechanism"] for row in rows).most_common(),
    }
    with Path(f"{args.prefix}.metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2)
        handle.write("\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
