#!/usr/bin/env python3
"""How closely the simulated interval sets resemble the real target.

For the real target and every simulated set, reports what the simulation was
built to match and what it was not:

  matched by design   interval count, interval-length distribution (by
                      decile medians), chromosome frequency
  not matched         merged components and total bp after merging adjacent
                      intervals, spacing between intervals, GC fraction, and
                      overlap with tandem repeats, segmental duplications and
                      low-mappability sequence (or any other BED given)

It also reproduces the per-chromosome allocation of simulate_targets.R, where
each chromosome contributes round(n_chr / 10) intervals per length decile
(R's round, half to even, which Python's round matches). Chromosomes with
fewer than five real intervals therefore contribute none, and the table shows
how many intervals and bp that leaves out.

BED coordinates are 0-based half-open. GC is computed over merged intervals,
ignoring N. Outputs:
  <prefix>.per_set.tsv          one row per interval set (real + simulated)
  <prefix>.summary.tsv          real value, simulated median and 2.5-97.5%, and
                                the real value's percentile among the simulations
  <prefix>.chromosomes.tsv      per chromosome: real count, per-decile allocation,
                                expected and median simulated count, and the
                                difference rounding makes
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pysam

SIMULATION_RE = re.compile(r"^simulation(\d+)\.bed$")


def read_intervals(path: Path) -> dict[str, np.ndarray]:
    opener = gzip.open if path.suffix == ".gz" else open
    by_chrom = defaultdict(list)
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end, *_ = line.split("\t")
            by_chrom[chrom].append((int(start), int(end)))
    return {chrom: np.array(sorted(values), dtype=np.int64) for chrom, values in by_chrom.items()}


def merge(intervals: np.ndarray) -> np.ndarray:
    if len(intervals) == 0:
        return intervals
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return np.array(merged, dtype=np.int64)


class Coverage:
    """Covered bp of a merged annotation within any [start, end)."""

    def __init__(self, intervals: np.ndarray):
        self.starts = intervals[:, 0] if len(intervals) else np.zeros(0, dtype=np.int64)
        self.ends = intervals[:, 1] if len(intervals) else np.zeros(0, dtype=np.int64)
        self.cumulative = np.concatenate([[0], np.cumsum(self.ends - self.starts)])

    def upto(self, positions: np.ndarray) -> np.ndarray:
        index = np.searchsorted(self.ends, positions, side="right")
        covered = self.cumulative[index]
        inside = index < len(self.starts)
        partial = np.where(inside, np.clip(positions - self.starts[np.minimum(index, len(self.starts) - 1)], 0, None), 0)
        return covered + partial

    def within(self, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
        if len(self.starts) == 0:
            return np.zeros(len(starts), dtype=np.int64)
        return self.upto(ends) - self.upto(starts)


def quantiles(values: np.ndarray, qs=(0.1, 0.5, 0.9)) -> list[float]:
    return [float(np.quantile(values, q)) if len(values) else float("nan") for q in qs]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--target-bed", required=True, type=Path, help="the real target, as given to the simulation")
    parser.add_argument("--simulation-dir", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--annotation", action="append", default=[], help="NAME=BED[.gz]; repeatable")
    parser.add_argument("--prefix", required=True, type=Path)
    args = parser.parse_args()

    sets = {"real": read_intervals(args.target_bed)}
    for path in sorted(args.simulation_dir.iterdir(), key=lambda p: p.name):
        match = SIMULATION_RE.match(path.name)
        if match:
            sets[f"simulation{match.group(1)}"] = read_intervals(path)
    merged_sets = {name: {chrom: merge(iv) for chrom, iv in chroms.items()} for name, chroms in sets.items()}

    annotations = {}
    for spec in args.annotation:
        name, path = spec.split("=", 1)
        annotations[name] = {chrom: Coverage(merge(iv)) for chrom, iv in read_intervals(Path(path)).items()}

    # GC per set, one chromosome at a time through prefix sums.
    gc_bp = defaultdict(int)
    acgt_bp = defaultdict(int)
    annotated_bp = defaultdict(int)
    chromosomes = sorted({chrom for chroms in merged_sets.values() for chrom in chroms})
    fasta = pysam.FastaFile(str(args.reference))
    for chrom in chromosomes:
        if chrom not in fasta.references:
            raise SystemExit(f"{chrom} is not in {args.reference}")
        sequence = np.frombuffer(fasta.fetch(chrom).upper().encode(), dtype=np.uint8)
        is_gc = np.concatenate([[0], np.cumsum((sequence == ord("G")) | (sequence == ord("C")), dtype=np.int64)])
        is_acgt = np.concatenate([[0], np.cumsum(sequence != ord("N"), dtype=np.int64)])
        del sequence
        for name, chroms in merged_sets.items():
            intervals = chroms.get(chrom)
            if intervals is None or len(intervals) == 0:
                continue
            starts, ends = intervals[:, 0], np.minimum(intervals[:, 1], len(is_gc) - 1)
            gc_bp[name] += int((is_gc[ends] - is_gc[starts]).sum())
            acgt_bp[name] += int((is_acgt[ends] - is_acgt[starts]).sum())
            for annotation, by_chrom in annotations.items():
                coverage = by_chrom.get(chrom)
                if coverage is not None:
                    annotated_bp[(name, annotation)] += int(coverage.within(intervals[:, 0], intervals[:, 1]).sum())

    rows = []
    for name, chroms in sets.items():
        lengths = np.concatenate([iv[:, 1] - iv[:, 0] for iv in chroms.values()]) if chroms else np.zeros(0)
        merged = merged_sets[name]
        merged_lengths = np.concatenate([iv[:, 1] - iv[:, 0] for iv in merged.values()])
        gaps = np.concatenate([iv[1:, 0] - iv[:-1, 1] for iv in merged.values() if len(iv) > 1])
        total_bp = int(merged_lengths.sum())
        p10, p50, p90 = quantiles(lengths)
        row = {
            "set": name,
            "intervals": int(len(lengths)),
            "merged_components": int(len(merged_lengths)),
            "total_bp": total_bp,
            "length_p10": p10, "length_median": p50, "length_p90": p90,
            "length_mean": float(lengths.mean()) if len(lengths) else float("nan"),
            "merged_length_median": float(np.median(merged_lengths)),
            "spacing_median": float(np.median(gaps)) if len(gaps) else float("nan"),
            "spacing_p10": quantiles(gaps)[0] if len(gaps) else float("nan"),
            "chromosomes": len(chroms),
            "gc_fraction": gc_bp[name] / acgt_bp[name] if acgt_bp[name] else float("nan"),
        }
        for annotation in annotations:
            row[f"{annotation}_fraction_bp"] = annotated_bp[(name, annotation)] / total_bp if total_bp else float("nan")
        rows.append(row)

    args.prefix.parent.mkdir(parents=True, exist_ok=True)
    with open(f"{args.prefix}.per_set.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    real = rows[0]
    simulated = rows[1:]
    summary = []
    for key in rows[0]:
        if key == "set":
            continue
        values = np.array([row[key] for row in simulated], dtype=float)
        values = values[~np.isnan(values)]
        summary.append({
            "feature": key, "real": real[key],
            "simulated_median": float(np.median(values)) if len(values) else float("nan"),
            "simulated_p2.5": float(np.quantile(values, 0.025)) if len(values) else float("nan"),
            "simulated_p97.5": float(np.quantile(values, 0.975)) if len(values) else float("nan"),
            "real_percentile": float((values < real[key]).mean() * 100 + (values == real[key]).mean() * 50) if len(values) else float("nan"),
            "n_simulated": len(values),
        })
    with open(f"{args.prefix}.summary.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)

    real_counts = {chrom: len(iv) for chrom, iv in sets["real"].items()}
    chromosome_rows = []
    for chrom in sorted(set(real_counts) | {c for name in sets for c in sets[name]}):
        n_real = real_counts.get(chrom, 0)
        per_decile = round(n_real / 10)
        simulated_counts = [len(sets[name].get(chrom, [])) for name in sets if name != "real"]
        chromosome_rows.append({
            "chrom": chrom, "real_intervals": n_real,
            "real_bp": int(sum(e - s for s, e in sets["real"].get(chrom, []))),
            "allocated_per_decile": per_decile, "expected_simulated": 10 * per_decile,
            "simulated_median": float(np.median(simulated_counts)) if simulated_counts else 0.0,
            # Negative where rounding drops intervals (all of them below 5), positive where it adds.
            "rounding_difference": 10 * per_decile - n_real,
        })
    with open(f"{args.prefix}.chromosomes.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(chromosome_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(chromosome_rows)
    print(f"{len(rows) - 1} simulated sets compared with the real target")


if __name__ == "__main__":
    main()
