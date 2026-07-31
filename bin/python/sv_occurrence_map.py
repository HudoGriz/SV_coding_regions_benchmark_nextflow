#!/usr/bin/env python3
"""Map every benchmark variant onto the interval sets it occupies.

The benchmark scores three nested targets separately, which answers "how well
does each pipeline do in each target" but not "which variants actually live
where, and are they the same variants". This builds the joint picture: one row
per truth variant, carrying its interval-set membership, its class and size, and
whether each pipeline recovered it.

Interval semantics match the benchmark exactly. A variant belongs to a target if
it intersects it by at least one base, which is what `truvari bench
--bench-overlaps` scores; the containment flag is emitted alongside so the two
can be compared. Regions are read the way `truvari.read_bed_tree` stores them -
inflated by one base, so the real half-open region is [begin, end - 1).

Outputs
-------
<prefix>.variants.tsv     one row per truth variant
<prefix>.crosstab.tsv     counts by interval-set combination x SV class
<prefix>.summary.txt      human-readable summary

Usage
-----
    sv_occurrence_map.py --results results-GRCh37 \
                         --truth  HG002_SVs_Tier1_v0.6.vcf.gz \
                         --hci    hci.bed --gp panel.bed --exutr exons.bed \
                         --prefix occurrence_GRCh37
"""
from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
from collections import Counter, defaultdict

import pandas as pd
from intervaltree import IntervalTree

SVTYPE_RE = re.compile(r'SVTYPE=(\w+)')
SVLEN_RE = re.compile(r'SVLEN=(-?\d+)')
END_RE = re.compile(r'[;\t]END=(\d+)')

# (technology, caller) as laid out under <results>/real_intervals/
PIPELINES = [
    ("Illumina_WES", "Manta"),
    ("Illumina_WGS", "Manta"),
    ("ONT", "CuteSV"),
    ("ONT", "Sniffles"),
    ("PacBio", "CuteSV"),
    ("PacBio", "Pbsv"),
]

SIZE_BINS = [0, 100, 300, 1000, 10 ** 9]
SIZE_NAMES = ["50-99bp", "100-299bp", "300-999bp", ">=1kb"]


def read_bed_tree(path: str) -> dict[str, IntervalTree]:
    """Read a BED the way truvari.read_bed_tree does: addi(start, end + 1)."""
    tree = defaultdict(IntervalTree)
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        for line in fh:
            if line.startswith(('#', 'track', 'browser')):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) >= 3:
                tree[f[0]].addi(int(f[1]), int(f[2]) + 1)
    return tree


def parse_vcf(path: str, passonly: bool = True):
    """Yield (chrom, start, end, svtype, svlen, key) with 0-based half-open coords."""
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split('\t', 8)
            if passonly and f[6] != 'PASS':
                continue
            chrom, pos, info = f[0], int(f[1]), f[7]
            m = SVLEN_RE.search(info)
            svlen = abs(int(m.group(1))) if m else 0
            m = SVTYPE_RE.search(info)
            svtype = m.group(1) if m else 'NA'
            m = END_RE.search(';' + info)
            end = int(m.group(1)) if m else pos + max(svlen, 1)
            start = pos - 1
            if svtype == 'INS':
                end = start + 1
            yield chrom, start, max(end, start + 1), svtype, svlen, (chrom, pos, svtype, svlen)


def membership(tree: dict, chrom: str, start: int, end: int) -> tuple[bool, bool]:
    """Return (intersects, contained) against a read_bed_tree-style tree."""
    if chrom not in tree:
        return False, False
    hits = tree[chrom].overlap(start, end)
    if not hits:
        return False, False
    intersects = any(max(start, i.begin) < min(end, i.end - 1) for i in hits)
    contained = any(start >= i.begin and end <= i.end - 1 for i in hits)
    return intersects, contained


def recall_labels(results: str, tech: str, caller: str, target: str = "high_confidence") -> dict:
    """Map truth-variant key -> True (TP) / False (FN) from a bench result."""
    base = os.path.join(results, "real_intervals", tech, "truvari", caller, target)
    stem = f"{tech}-{caller}-{target}"
    out = {}
    for fname, label in ((f"{stem}.tp-base.vcf.gz", True), (f"{stem}.fn.vcf.gz", False)):
        p = os.path.join(base, fname)
        if not os.path.exists(p):
            return {}
        for *_, key in parse_vcf(p, passonly=False):
            out[key] = label
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results', required=True, help='pipeline results directory')
    ap.add_argument('--truth', required=True)
    ap.add_argument('--hci', required=True)
    ap.add_argument('--gp', required=True)
    ap.add_argument('--exutr', required=True)
    ap.add_argument('--prefix', required=True)
    ap.add_argument('--sizemin', type=int, default=50)
    ap.add_argument('--sizemax', type=int, default=50000)
    args = ap.parse_args()

    trees = {'HCI': read_bed_tree(args.hci),
             'GP': read_bed_tree(args.gp),
             'EXUTR': read_bed_tree(args.exutr)}

    labels = {}
    for tech, caller in PIPELINES:
        lab = recall_labels(args.results, tech, caller)
        if lab:
            labels[f"{tech}|{caller}"] = lab
    if not labels:
        print(f"warning: no bench results found under {args.results}/real_intervals",
              file=sys.stderr)

    rows = []
    for chrom, start, end, svtype, svlen, key in parse_vcf(args.truth):
        if svlen < args.sizemin or (args.sizemax > 0 and svlen > args.sizemax):
            continue
        row = {'chrom': chrom, 'start': start, 'end': end,
               'svtype': svtype, 'svlen': svlen}
        for name, tree in trees.items():
            inter, cont = membership(tree, chrom, start, end)
            row[f'in_{name}'] = inter
            row[f'contained_{name}'] = cont
        # a variant that intersects but is not contained spans a boundary and is
        # only scored because of --bench-overlaps
        for name in trees:
            row[f'spans_{name}'] = row[f'in_{name}'] and not row[f'contained_{name}']
        for pipe, lab in labels.items():
            row[f'recalled::{pipe}'] = lab.get(key)
        rows.append(row)

    df = pd.DataFrame(rows)
    if df.empty:
        print("no variants parsed", file=sys.stderr)
        return 1
    df['size_class'] = pd.cut(df['svlen'], bins=SIZE_BINS, labels=SIZE_NAMES, right=False)

    var_path = f"{args.prefix}.variants.tsv"
    df.to_csv(var_path, sep='\t', index=False)

    # combination of interval sets each variant occupies
    def combo(r):
        parts = [n for n in ('HCI', 'GP', 'EXUTR') if r[f'in_{n}']]
        return '+'.join(parts) if parts else 'none'
    df['combination'] = df.apply(combo, axis=1)

    ct = (df.groupby(['combination', 'svtype', 'size_class'], observed=True)
            .size().rename('n').reset_index())
    ct_path = f"{args.prefix}.crosstab.tsv"
    ct.to_csv(ct_path, sep='\t', index=False)

    lines = []
    add = lines.append
    add(f"SV occurrence map — {args.results}")
    add(f"truth: {os.path.basename(args.truth)}  "
        f"({args.sizemin} <= SVLEN <= {args.sizemax}, PASS only)")
    add(f"total truth variants: {len(df)}")
    add("")
    add("Interval-set occupancy (>=1bp intersection, as scored):")
    for name in ('HCI', 'GP', 'EXUTR'):
        n = int(df[f'in_{name}'].sum())
        c = int(df[f'contained_{name}'].sum())
        s = int(df[f'spans_{name}'].sum())
        add(f"  {name:6s} intersecting {n:6d}   contained {c:6d}   boundary-spanning {s:5d}")
    add("")
    add("Combinations:")
    for k, v in df['combination'].value_counts().items():
        add(f"  {k:20s} {v:6d}")
    add("")
    add("Deletion share by target:")
    for name in ('HCI', 'GP', 'EXUTR'):
        sub = df[df[f'in_{name}']]
        if len(sub):
            add(f"  {name:6s} {(sub['svtype'] == 'DEL').mean() * 100:5.1f}%  (n={len(sub)})")
    add("")
    add("Size composition by target (% of variants in that target):")
    add(f"  {'target':8s}" + "".join(f"{s:>12s}" for s in SIZE_NAMES))
    for name in ('HCI', 'GP', 'EXUTR'):
        sub = df[df[f'in_{name}']]
        if not len(sub):
            continue
        frac = sub['size_class'].value_counts(normalize=True).reindex(SIZE_NAMES).fillna(0) * 100
        add(f"  {name:8s}" + "".join(f"{frac[s]:11.1f}%" for s in SIZE_NAMES))
    if labels:
        add("")
        add("Recall by target and pipeline (from the HCI bench labels, so this is")
        add("per-variant recallability re-tabulated by where the variant lives):")
        add(f"  {'pipeline':22s}" + "".join(f"{n:>10s}" for n in ('HCI', 'GP', 'EXUTR')))
        for pipe in labels:
            col = f'recalled::{pipe}'
            cells = []
            for name in ('HCI', 'GP', 'EXUTR'):
                sub = df[df[f'in_{name}'] & df[col].notna()]
                cells.append(f"{sub[col].mean() * 100:9.1f}%" if len(sub) else "        -")
            add(f"  {pipe:22s}" + "".join(cells))

    summary = "\n".join(lines)
    with open(f"{args.prefix}.summary.txt", 'w') as fh:
        fh.write(summary + "\n")
    print(summary)
    print(f"\nwrote {var_path}, {ct_path}, {args.prefix}.summary.txt")
    return 0


if __name__ == '__main__':
    sys.exit(main())
