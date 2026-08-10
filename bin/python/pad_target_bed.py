#!/usr/bin/env python3
"""Pad a target BED, clip it to an allowed BED, and merge the result."""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path


def read_bed(path: Path):
    result = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end, *_ = line.rstrip("\n").split("\t")
            result[chrom].append((int(start), int(end)))
    for chrom in result:
        result[chrom].sort()
    return result


def merge(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", required=True, type=Path)
    parser.add_argument("--allowed", required=True, type=Path)
    parser.add_argument("--padding", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    target = read_bed(args.target)
    allowed = read_bed(args.allowed)
    output = defaultdict(list)
    for chrom, intervals in target.items():
        bounds = allowed.get(chrom, [])
        bound_index = 0
        for start, end in intervals:
            padded_start = max(0, start - args.padding)
            padded_end = end + args.padding
            while bound_index < len(bounds) and bounds[bound_index][1] <= padded_start:
                bound_index += 1
            index = bound_index
            while index < len(bounds) and bounds[index][0] < padded_end:
                left = max(padded_start, bounds[index][0])
                right = min(padded_end, bounds[index][1])
                if right > left:
                    output[chrom].append((left, right))
                index += 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        for chrom in sorted(output):
            for start, end in merge(output[chrom]):
                handle.write(f"{chrom}\t{start}\t{end}\n")


if __name__ == "__main__":
    main()

