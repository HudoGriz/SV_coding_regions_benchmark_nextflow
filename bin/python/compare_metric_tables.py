#!/usr/bin/env python3
"""Compare two gather-statistics real-interval tables by pipeline and target."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


KEYS = ("tech", "caller", "range")
METRICS = ("precision", "recall", "f1")
COUNTS = ("TP.base", "FP", "FN", "base.cnt", "comp.cnt")


def read_table(path: Path) -> dict[tuple[str, ...], dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    indexed = {tuple(row[key] for key in KEYS): row for row in rows}
    if len(indexed) != len(rows):
        raise ValueError(f"duplicate technology/caller/target rows in {path}")
    return indexed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old", required=True, type=Path)
    parser.add_argument("--new", required=True, type=Path)
    parser.add_argument("--old-label", default="old")
    parser.add_argument("--new-label", default="new")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    old = read_table(args.old)
    new = read_table(args.new)
    keys = sorted(set(old) | set(new))
    fields = [*KEYS, "status"]
    for metric in METRICS:
        fields += [f"{metric}_{args.old_label}", f"{metric}_{args.new_label}", f"delta_{metric}_pp"]
    for count in COUNTS:
        safe = count.replace(".", "_")
        fields += [f"{safe}_{args.old_label}", f"{safe}_{args.new_label}", f"delta_{safe}"]

    output = []
    for key in keys:
        before, after = old.get(key), new.get(key)
        row = dict(zip(KEYS, key))
        row["status"] = "matched" if before and after else ("added" if after else "removed")
        for metric in METRICS:
            before_value = float(before[metric]) if before else None
            after_value = float(after[metric]) if after else None
            row[f"{metric}_{args.old_label}"] = before_value
            row[f"{metric}_{args.new_label}"] = after_value
            row[f"delta_{metric}_pp"] = (
                100 * (after_value - before_value)
                if before_value is not None and after_value is not None else ""
            )
        for count in COUNTS:
            safe = count.replace(".", "_")
            before_value = int(before[count]) if before else None
            after_value = int(after[count]) if after else None
            row[f"{safe}_{args.old_label}"] = before_value
            row[f"{safe}_{args.new_label}"] = after_value
            row[f"delta_{safe}"] = (
                after_value - before_value
                if before_value is not None and after_value is not None else ""
            )
        output.append(row)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output)


if __name__ == "__main__":
    main()
