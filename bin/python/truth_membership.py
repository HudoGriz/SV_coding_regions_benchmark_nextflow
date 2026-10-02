#!/usr/bin/env python3
"""Truth records a target admits by containment and by any overlap, counted directly.

The real targets were benchmarked under both membership rules (the containment
sensitivity setting), so Truvari gives their counts. This script adds two counts
that need no benchmark:

* simulated sets: the truth records each simulated interval set admits under
  either rule, with Truvari's own boundary conventions, so the simulated sets can
  be set beside the real targets without 500 more benchmarks per rule.
    - any overlap (the --bench-overlaps fork): a record spans [POS-1, POS-1+len(REF));
      an insertion therefore occupies one base. Admitted when it shares a base
      with any interval.
    - containment (stock --includebed): an insertion's span ends one base later,
      [POS-1, POS+1). Admitted when one interval contains the whole span.
  The conventions are checked against Truvari's counts on the real target (from
  metrics.tsv) before any simulated set is counted; the script stops if they
  disagree.
* independent intersection: the real target counted with plain VCF-BED
  conventions instead ([POS-1, POS-1+len(REF)) for every record except an
  insertion, which occupies the one base [POS-1, POS)). It lists the records that
  any overlap admits and containment does not, with their type and length.

The truth records are those an HCI benchmark scored (tp-base plus fn); the truth
side is the same for every pipeline. The target BED is the real target clipped
to HCI and merged, as the benchmarks use it (pad_target_bed.py --padding 0).

    truth_membership.py --tp-base HCI.tp-base.vcf.gz --fn HCI.fn.vcf.gz \\
        --target-bed TARGET.bed --simulation-dir DIR --metrics metrics.tsv \\
        --target wes_utr --prefix membership

writes <prefix>.simulated.tsv (one row per simulated set), <prefix>.summary.tsv
and <prefix>.boundary_records.tsv.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path

SVTYPE = re.compile(r"(?:^|;)SVTYPE=([^;]+)")
SVLEN = re.compile(r"(?:^|;)SVLEN=(-?\d+)")


def truth_records(paths):
    records = []
    for path in paths:
        with gzip.open(path, "rt") as handle:
            for line in handle:
                if line[0] == "#":
                    continue
                chrom, pos, _, ref, alt, _, _, info = line.split("\t", 8)[:8]
                svtype = SVTYPE.search(info)
                svlen = SVLEN.search(info)
                records.append({
                    "chrom": chrom, "pos": int(pos), "ref_len": len(ref), "alt_len": len(alt),
                    "svtype": svtype.group(1) if svtype else "NA",
                    "svlen": abs(int(svlen.group(1))) if svlen else None,
                })
    return records


def read_bed(path):
    by_chrom = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end = line.split("\t")[:3]
            by_chrom[chrom].append((int(start), int(end)))
    index = {}
    for chrom, intervals in by_chrom.items():
        intervals.sort()
        starts = [s for s, _ in intervals]
        reach, top = [], -1
        for _, end in intervals:
            top = max(top, end)
            reach.append(top)          # furthest end among intervals starting at or before this one
        index[chrom] = (starts, reach)
    return index


def overlaps(index, chrom, start, end):
    """Any interval shares a base with [start, end)."""
    starts, reach = index.get(chrom, ((), ()))
    j = bisect.bisect_left(starts, end) - 1          # last interval starting before end
    return j >= 0 and reach[j] > start


def contained(index, chrom, start, end):
    """One interval holds all of [start, end). Merged intervals make reach exact."""
    starts, reach = index.get(chrom, ((), ()))
    i = bisect.bisect_right(starts, start) - 1       # intervals starting at or before start
    return i >= 0 and reach[i] >= end


def is_insertion(record):
    return record["ref_len"] == 1 and record["alt_len"] > 1


def truvari_counts(records, index):
    """(any overlap, containment) with Truvari's conventions."""
    n_overlap = n_contained = 0
    for r in records:
        start = r["pos"] - 1
        n_overlap += overlaps(index, r["chrom"], start, start + r["ref_len"])
        n_contained += contained(index, r["chrom"], start, r["pos"] + 1 if is_insertion(r) else start + r["ref_len"])
    return n_overlap, n_contained


def independent_counts(records, index):
    """Plain VCF-BED intersection; returns (overlap, contained, records gained by overlap)."""
    n_overlap = n_contained = 0
    gained = []
    for r in records:
        start = r["pos"] - 1
        end = r["pos"] if r["svtype"] == "INS" else start + r["ref_len"]
        o = overlaps(index, r["chrom"], start, end)
        c = contained(index, r["chrom"], start, end)
        n_overlap += o
        n_contained += c
        if o and not c:
            gained.append(r)
    return n_overlap, n_contained, gained


def truvari_truth_counts(metrics, target):
    """(any overlap, containment) truth denominators of the real target in metrics.tsv."""
    found = {}
    with open(metrics) as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            if r["target_set"] == "real" and r["target"] == target and r["setting"] in ("primary", "containment"):
                found.setdefault(r["setting"], set()).add(int(r["truth_denominator"]))
    for setting, values in found.items():
        if len(values) != 1:
            raise SystemExit(f"{target} {setting}: truth denominators differ between pipelines: {sorted(values)}")
    if "primary" not in found or "containment" not in found:
        raise SystemExit(f"metrics.tsv lacks the primary or containment benchmark of {target}")
    return found["primary"].pop(), found["containment"].pop()


def percentile(values, q):
    """Linear interpolation between order statistics (numpy's default)."""
    v = sorted(values)
    k = (len(v) - 1) * q / 100
    lo = int(k)
    hi = min(lo + 1, len(v) - 1)
    return v[lo] + (v[hi] - v[lo]) * (k - lo)


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tp-base", required=True, type=Path)
    parser.add_argument("--fn", required=True, type=Path)
    parser.add_argument("--target-bed", required=True, type=Path)
    parser.add_argument("--simulation-dir", required=True, type=Path)
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--target", default="wes_utr")
    parser.add_argument("--prefix", default="membership")
    args = parser.parse_args()

    records = truth_records([args.tp_base, args.fn])
    target = read_bed(args.target_bed)

    want = truvari_truth_counts(args.metrics, args.target)
    got = truvari_counts(records, target)
    if got != want:
        raise SystemExit(f"{args.target}: direct count (overlap {got[0]}, containment {got[1]}) differs from "
                         f"Truvari (overlap {want[0]}, containment {want[1]})")

    sims = sorted(args.simulation_dir.glob("simulation*.bed"), key=lambda p: int(p.stem.replace("simulation", "")))
    rows, increases = [], []
    for bed in sims:
        n_overlap, n_contained = truvari_counts(records, read_bed(bed))
        increases.append(100 * (n_overlap / n_contained - 1))
        rows.append([bed.stem, n_overlap, n_contained, f"{increases[-1]:.4f}"])
    write_tsv(f"{args.prefix}.simulated.tsv", ["set", "overlap", "containment", "increase_pct"], rows)

    ind_overlap, ind_contained, gained = independent_counts(records, target)
    lengths = sorted(r["svlen"] for r in gained if r["svlen"] is not None)
    types = Counter(r["svtype"] for r in gained)
    write_tsv(f"{args.prefix}.boundary_records.tsv", ["chrom", "pos", "svtype", "svlen"],
              [[r["chrom"], r["pos"], r["svtype"], "" if r["svlen"] is None else r["svlen"]] for r in gained])

    summary = [
        ("truvari_conventions", "target_overlap", want[0]),
        ("truvari_conventions", "target_containment", want[1]),
        ("truvari_conventions", "target_increase_pct", f"{100 * (want[0] / want[1] - 1):.4f}"),
        ("independent_intersection", "target_overlap", ind_overlap),
        ("independent_intersection", "target_containment", ind_contained),
        ("independent_intersection", "gained_by_overlap", len(gained)),
        ("independent_intersection", "gained_increase_pct", f"{100 * (ind_overlap / ind_contained - 1):.4f}"),
        ("independent_intersection", "gained_types", ",".join(f"{t}:{n}" for t, n in sorted(types.items()))),
        ("independent_intersection", "gained_svlen_median", statistics.median(lengths) if lengths else ""),
        ("independent_intersection", "gained_svlen_min", lengths[0] if lengths else ""),
        ("independent_intersection", "gained_svlen_max", lengths[-1] if lengths else ""),
        ("simulated", "sets", len(rows)),
    ]
    if rows:
        summary += [
            ("simulated", "overlap_median", statistics.median(r[1] for r in rows)),
            ("simulated", "containment_median", statistics.median(r[2] for r in rows)),
            ("simulated", "increase_pct_median", f"{statistics.median(increases):.4f}"),
            ("simulated", "increase_pct_p2.5", f"{percentile(increases, 2.5):.4f}"),
            ("simulated", "increase_pct_p97.5", f"{percentile(increases, 97.5):.4f}"),
        ]
    write_tsv(f"{args.prefix}.summary.tsv", ["count", "quantity", "value"], summary)
    print(f"{args.target}: overlap {want[0]}, containment {want[1]} (checked against Truvari); "
          f"independent {ind_overlap}/{ind_contained}, {len(gained)} gained; {len(rows)} simulated sets")


if __name__ == "__main__":
    main()
