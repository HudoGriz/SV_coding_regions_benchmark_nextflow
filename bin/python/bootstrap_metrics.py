#!/usr/bin/env python3
"""Sampling uncertainty of the target metrics, and a like-for-like null.

Two problems with the raw percentile ranks: the observed target metric is a
single estimate with no interval, and each simulated set scores roughly twice
as many truth variants as the real target, so the simulated distribution is
narrower than the observed value's own sampling spread.

1. Block bootstrap of the observed target. Records are grouped by the merged
   target component they overlap, because SVs in the same gene are not
   independent. Components (with every truth and candidate record in them) are
   resampled with replacement; precision, recall and F1 are recomputed each
   time. Reports the 2.5 and 97.5 percentiles.

2. Rarefied simulated sets. Each simulated set is subsampled, by whole
   components in a random order, until it holds at least as many truth records
   as the real target. That gives a simulated distribution with the target's
   own sample size, and the target's percentile rank is recomputed against it.

Records outside every component (possible at an insertion's end boundary,
where the membership conventions differ by a base) join their nearest component.
Every draw uses a seeded generator; the seeds are written to the output.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from target_transition_audit import find_bench, read_bed, read_benchmark, read_simulation_benchmark  # noqa: E402

SIMULATION_RE = re.compile(r"^simulation(\d+)\.bed$")
COLUMNS = ("tp_base", "fn", "tp_comp", "fp")


def component_counts(bench: dict, bed) -> np.ndarray:
    """Rows = components, columns = TP-base, FN, TP-comp, FP."""
    index = {}
    starts_by_chrom = {chrom: starts for chrom, (starts, _) in bed.items()}
    ends_by_chrom = {chrom: ends for chrom, (_, ends) in bed.items()}

    def component(record) -> tuple:
        starts, ends = starts_by_chrom.get(record.chrom), ends_by_chrom.get(record.chrom)
        if not starts:
            return (record.chrom, -1)
        i = bisect.bisect_right(ends, record.start)
        if i < len(starts) and starts[i] < record.end:
            return (record.chrom, i)
        # Outside every component: take the nearer neighbour.
        left, right = i - 1, i
        if right >= len(starts):
            return (record.chrom, left)
        if left < 0:
            return (record.chrom, right)
        return (record.chrom, left if record.start - ends[left] <= starts[right] - record.end else right)

    rows = []
    for column, name in enumerate(COLUMNS):
        for record in bench[name].values():
            key = component(record)
            if key not in index:
                index[key] = len(rows)
                rows.append([0, 0, 0, 0])
            rows[index[key]][column] += 1
    return np.array(rows, dtype=np.int64) if rows else np.zeros((0, 4), dtype=np.int64)


def metrics(totals: np.ndarray) -> np.ndarray:
    """totals: (..., 4) -> (..., 3) precision, recall, F1."""
    tp_base, fn, tp_comp, fp = (totals[..., i].astype(float) for i in range(4))
    with np.errstate(invalid="ignore", divide="ignore"):
        recall = tp_base / (tp_base + fn)
        precision = tp_comp / (tp_comp + fp)
        f1 = 2 * precision * recall / (precision + recall)
    return np.stack([precision, recall, f1], axis=-1)


def rarefy(components: np.ndarray, truth_target: int, rng: np.random.Generator) -> np.ndarray:
    order = rng.permutation(len(components))
    truth = np.cumsum(components[order, 0] + components[order, 1])
    stop = int(np.searchsorted(truth, truth_target)) + 1
    return components[order[:stop]].sum(axis=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--pipeline", required=True, help="TECH:CALLER")
    parser.add_argument("--target", default="wes_utr")
    parser.add_argument("--target-bed", required=True, type=Path)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--rarefy-draws", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260928)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    technology, caller = args.pipeline.split(":")
    observed = component_counts(read_benchmark(find_bench(args.results, technology, caller, args.target)),
                                read_bed(args.target_bed))
    observed_totals = observed.sum(axis=0)
    observed_metrics = metrics(observed_totals)
    truth_target = int(observed_totals[0] + observed_totals[1])

    rng = np.random.default_rng(args.seed)
    draws = rng.integers(0, len(observed), size=(args.bootstrap, len(observed)))
    boot = metrics(observed[draws].sum(axis=1))

    bed_dir = args.results / "simulations" / "target_regions"
    simulations = sorted(
        (int(m.group(1)), path) for path in bed_dir.iterdir() if (m := SIMULATION_RE.match(path.name))
    )
    rarefied, full = [], []
    for number, bed_path in simulations:
        bench = read_simulation_benchmark(args.results, technology, caller, f"simulation{number}")
        components = component_counts(bench, read_bed(bed_path))
        full.append(metrics(components.sum(axis=0)))
        sim_rng = np.random.default_rng([args.seed, number])
        for _ in range(args.rarefy_draws):
            rarefied.append(metrics(rarefy(components, truth_target, sim_rng)))
    rarefied, full = np.array(rarefied), np.array(full)

    rows = []
    for i, metric in enumerate(("precision", "recall", "f1")):
        value = observed_metrics[i]
        def rank(values):
            values = values[~np.isnan(values)]
            return float((values < value).mean() * 100 + (values == value).mean() * 50)
        rows.append({
            "assembly": args.assembly, "pipeline": f"{technology} {caller}", "target": args.target,
            "metric": metric, "observed": float(value),
            "bootstrap_ci_low": float(np.nanpercentile(boot[:, i], 2.5)),
            "bootstrap_ci_high": float(np.nanpercentile(boot[:, i], 97.5)),
            "observed_truth_n": truth_target, "components": len(observed),
            "simulated_median_full": float(np.nanmedian(full[:, i])),
            "percentile_full": rank(full[:, i]),
            "simulated_median_rarefied": float(np.nanmedian(rarefied[:, i])),
            "simulated_p2.5_rarefied": float(np.nanpercentile(rarefied[:, i], 2.5)),
            "simulated_p97.5_rarefied": float(np.nanpercentile(rarefied[:, i], 97.5)),
            "percentile_rarefied": rank(rarefied[:, i]),
            "bootstrap_draws": args.bootstrap, "rarefy_draws_per_set": args.rarefy_draws,
            "n_simulated": len(simulations), "seed": args.seed,
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"{technology} {caller}: bootstrap and rarefaction written")


if __name__ == "__main__":
    main()
