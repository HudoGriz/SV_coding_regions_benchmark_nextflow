#!/usr/bin/env python3
"""Attribute the EX+UTR performance deficit to measurable factors.

The benchmark shows lower metrics in exons+UTRs than genome-wide. Three things
could produce that, and they are separable:

  1. Composition - the exonic target holds a different mix of variants
     (more deletions, more >=1kb) than the genome as a whole.
  2. Scoring geometry - the same variant, scored inside a small fragmented
     target rather than genome-wide, can flip its outcome: breakpoints fall
     outside the interval, matching candidates are truncated away.
  3. Genomic context - exonic sequence is genuinely harder to call.

Step 1 decomposes the observed gap into (1) and "everything else", by taking
the per-variant true-positive / false-negative labels Truvari assigned in the
HCI benchmark and reweighting them onto the exonic variant set. That asks what
the metric would be if only the mix of scored variants changed.

Step 2 attacks "everything else". For variants present in both benchmarks it
identifies discordant calls - the same variant labelled TP genome-wide and FN in
the exonic target, or the reverse - and relates discordance to covariates that
can be computed from the reference and the annotation shipped with the run:
SV class and size, whether the variant spans a target boundary, the size of the
containing interval, distance to the nearest interval edge, GC content, tandem
repeat overlap, and local SV density.

Mappability and segmental-duplication tracks are not available in this run
directory, so those two covariates are absent rather than estimated. Tandem
repeat overlap is the closest available proxy for hard-to-align sequence.

Model: stratified discordance rates are the primary output because they are
directly interpretable. A multivariable logistic regression, fitted with
iteratively reweighted least squares in numpy, is reported alongside to show
which associations survive adjustment for the others. No scipy or sklearn
dependency, so this runs in the same container as the benchmark.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam
from intervaltree import IntervalTree

SVTYPE_RE = re.compile(r'SVTYPE=(\w+)')
SVLEN_RE = re.compile(r'SVLEN=(-?\d+)')
END_RE = re.compile(r'[;\t]END=(\d+)')

PIPELINES = [
    ("Illumina_WGS", "Manta"),
    ("ONT", "CuteSV"),
    ("ONT", "Sniffles"),
    ("PacBio", "CuteSV"),
    ("PacBio", "Pbsv"),
]

GC_FLANK = 500  # bp either side, for local GC


# ---------------------------------------------------------------- inputs

def read_bed_tree(path):
    """Read a BED as truvari.read_bed_tree does: addi(start, end + 1)."""
    tree = defaultdict(IntervalTree)
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as fh:
        for line in fh:
            if line.startswith(('#', 'track', 'browser')):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) >= 3:
                tree[f[0]].addi(int(f[1]), int(f[2]) + 1)
    return tree


def parse_vcf(path, passonly=True):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as fh:
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


def labels_for(results, tech, caller, target):
    """key -> True (TP) / False (FN) from one bench result."""
    base = os.path.join(results, "real_intervals", tech, "truvari", caller, target)
    stem = f"{tech}-{caller}-{target}"
    out = {}
    for fn, lab in ((f"{stem}.tp-base.vcf.gz", True), (f"{stem}.fn.vcf.gz", False)):
        p = os.path.join(base, fn)
        if not os.path.exists(p):
            return {}
        for *_, key in parse_vcf(p, passonly=False):
            out[key] = lab
    return out


# ---------------------------------------------------------------- covariates

def interval_context(tree, chrom, start, end):
    """(intersects, contained, containing interval length, distance to edge)."""
    if chrom not in tree:
        return False, False, np.nan, np.nan
    hits = list(tree[chrom].overlap(start, end))
    if not hits:
        return False, False, np.nan, np.nan
    inter = any(max(start, i.begin) < min(end, i.end - 1) for i in hits)
    if not inter:
        return False, False, np.nan, np.nan
    cont = any(start >= i.begin and end <= i.end - 1 for i in hits)
    # the interval the variant starts in, else the first it touches
    best = min(hits, key=lambda i: abs(i.begin - start))
    ilen = (best.end - 1) - best.begin
    edge = min(abs(start - best.begin), abs(end - (best.end - 1)))
    return inter, cont, ilen, edge


def gc_fraction(fasta, chrom, start, end, flank=GC_FLANK):
    lo = max(0, start - flank)
    try:
        seq = fasta.fetch(chrom, lo, end + flank).upper()
    except (KeyError, ValueError):
        return np.nan
    if not seq:
        return np.nan
    gc = seq.count('G') + seq.count('C')
    acgt = gc + seq.count('A') + seq.count('T')
    return gc / acgt if acgt else np.nan


def tr_overlap_fraction(tree, chrom, start, end):
    if chrom not in tree:
        return 0.0
    span = max(end - start, 1)
    covered = 0
    for i in tree[chrom].overlap(start, end):
        covered += max(0, min(end, i.end - 1) - max(start, i.begin))
    return min(covered / span, 1.0)


def local_density(positions_by_chrom, chrom, start, window=100_000):
    arr = positions_by_chrom.get(chrom)
    if arr is None or not len(arr):
        return np.nan
    lo = np.searchsorted(arr, start - window // 2)
    hi = np.searchsorted(arr, start + window // 2)
    return float(hi - lo)


# ---------------------------------------------------------------- model

def logistic_irls(X, y, max_iter=50, tol=1e-8, ridge=1e-6):
    """Plain IRLS logistic regression. Returns (beta, se)."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    n, p = X.shape
    beta = np.zeros(p)
    for _ in range(max_iter):
        eta = np.clip(X @ beta, -30, 30)
        mu = 1.0 / (1.0 + np.exp(-eta))
        w = np.clip(mu * (1 - mu), 1e-9, None)
        z = eta + (y - mu) / w
        XtW = X.T * w
        H = XtW @ X + ridge * np.eye(p)
        try:
            new = np.linalg.solve(H, XtW @ z)
        except np.linalg.LinAlgError:
            return beta, np.full(p, np.nan)
        if np.max(np.abs(new - beta)) < tol:
            beta = new
            break
        beta = new
    eta = np.clip(X @ beta, -30, 30)
    mu = 1.0 / (1.0 + np.exp(-eta))
    w = np.clip(mu * (1 - mu), 1e-9, None)
    H = (X.T * w) @ X + ridge * np.eye(p)
    try:
        se = np.sqrt(np.diag(np.linalg.inv(H)))
    except np.linalg.LinAlgError:
        se = np.full(p, np.nan)
    return beta, se


# ---------------------------------------------------------------- main

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results', required=True)
    ap.add_argument('--truth', required=True)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--hci', required=True)
    ap.add_argument('--exutr', required=True)
    ap.add_argument('--trf', default=None, help='tandem repeat BED')
    ap.add_argument('--prefix', required=True)
    ap.add_argument('--sizemin', type=int, default=50)
    ap.add_argument('--sizemax', type=int, default=50000)
    args = ap.parse_args()

    hci_tree = read_bed_tree(args.hci)
    ex_tree = read_bed_tree(args.exutr)
    trf_tree = read_bed_tree(args.trf) if args.trf and os.path.exists(args.trf) else None
    fasta = pysam.FastaFile(args.fasta)

    variants = [v for v in parse_vcf(args.truth)
                if args.sizemin <= v[4] and (args.sizemax <= 0 or v[4] <= args.sizemax)]
    pos_by_chrom = {}
    for c in {v[0] for v in variants}:
        pos_by_chrom[c] = np.sort(np.array([v[1] for v in variants if v[0] == c]))

    rows = []
    for chrom, start, end, svtype, svlen, key in variants:
        in_hci, _, _, _ = interval_context(hci_tree, chrom, start, end)
        if not in_hci:
            continue
        in_ex, cont_ex, ilen, edge = interval_context(ex_tree, chrom, start, end)
        rows.append({
            'key': key, 'chrom': chrom, 'start': start, 'end': end,
            'svtype': svtype, 'svlen': svlen,
            'in_exutr': in_ex,
            'spans_boundary': bool(in_ex and not cont_ex),
            'interval_len': ilen, 'edge_dist': edge,
            'gc': gc_fraction(fasta, chrom, start, end),
            'tr_frac': tr_overlap_fraction(trf_tree, chrom, start, end) if trf_tree else np.nan,
            'sv_density_100kb': local_density(pos_by_chrom, chrom, start),
        })
    df = pd.DataFrame(rows)
    if df.empty:
        print("no variants within HCI", file=sys.stderr)
        return 1

    out = []
    add = out.append
    add(f"EX+UTR factor attribution — {args.results}")
    add(f"HCI truth variants: {len(df)};  also in EX+UTR: {int(df['in_exutr'].sum())}")
    add("")

    # ---------------- step 1: composition vs everything else
    add("=" * 78)
    add("STEP 1  Decomposition of the recall gap")
    add("=" * 78)
    add("")
    add("  observed HCI      recall over all HCI variants, HCI benchmark labels")
    add("  composition-only  the same labels, restricted to the EX+UTR variants;")
    add("                    i.e. what recall would be if only the variant mix changed")
    add("  observed EX+UTR   the EX+UTR benchmark's own recall")
    add("")
    add(f"  {'pipeline':22s} {'HCI':>8s} {'comp-only':>10s} {'EX+UTR':>8s} "
        f"{'composition':>12s} {'residual':>10s}")
    decomp = []
    for tech, caller in PIPELINES:
        hl = labels_for(args.results, tech, caller, 'high_confidence')
        el = labels_for(args.results, tech, caller, 'wes_utr')
        if not hl or not el:
            continue
        name = f"{tech} {caller}"
        hci_lab = df['key'].map(hl)
        ex_mask = df['in_exutr'] & hci_lab.notna()
        r_hci = hci_lab.mean()
        r_comp = hci_lab[ex_mask].mean()
        ex_lab = df.loc[df['in_exutr'], 'key'].map(el)
        r_ex = ex_lab.mean()
        comp_eff = (r_comp - r_hci) * 100
        resid = (r_ex - r_comp) * 100
        decomp.append({'pipeline': name, 'hci': r_hci, 'comp_only': r_comp,
                       'exutr': r_ex, 'composition_pp': comp_eff, 'residual_pp': resid})
        add(f"  {name:22s} {r_hci:8.3f} {r_comp:10.3f} {r_ex:8.3f} "
            f"{comp_eff:+11.1f} {resid:+9.1f}")
    add("")
    add("  composition = how much of the gap is explained by which variants the")
    add("  target contains; residual = everything else (scoring geometry plus any")
    add("  genuine context effect).")

    # ---------------- step 2: discordance
    add("")
    add("=" * 78)
    add("STEP 2  What predicts a variant flipping outcome between the two targets")
    add("=" * 78)
    sub = df[df['in_exutr']].copy()
    disc_frames = []
    for tech, caller in PIPELINES:
        hl = labels_for(args.results, tech, caller, 'high_confidence')
        el = labels_for(args.results, tech, caller, 'wes_utr')
        if not hl or not el:
            continue
        s = sub.copy()
        s['hci_tp'] = s['key'].map(hl)
        s['ex_tp'] = s['key'].map(el)
        s = s[s['hci_tp'].notna() & s['ex_tp'].notna()]
        s['pipeline'] = f"{tech} {caller}"
        # lost = recovered genome-wide but missed in the exonic benchmark
        s['lost'] = s['hci_tp'] & ~s['ex_tp']
        disc_frames.append(s)
    if not disc_frames:
        add("  no paired labels available")
        print("\n".join(out))
        return 0
    D = pd.concat(disc_frames, ignore_index=True)
    add("")
    add(f"  paired variant-by-pipeline observations: {len(D)}")
    add(f"  lost between HCI and EX+UTR:             {int(D['lost'].sum())} "
        f"({D['lost'].mean() * 100:.1f}%)")
    add("")

    def strat(col, bins, names, label):
        add(f"  {label}")
        D['_b'] = pd.cut(D[col], bins=bins, labels=names, right=False) \
            if not pd.api.types.is_bool_dtype(D[col]) else D[col]
        g = D.groupby('_b', observed=True)['lost'].agg(['mean', 'size'])
        for k, r in g.iterrows():
            add(f"    {str(k):>14s}  lost {r['mean'] * 100:5.1f}%   n={int(r['size'])}")
        add("")

    strat('svtype', None, None, 'by SV class') if False else None
    add("  by SV class")
    for k, r in D.groupby('svtype')['lost'].agg(['mean', 'size']).iterrows():
        add(f"    {k:>14s}  lost {r['mean'] * 100:5.1f}%   n={int(r['size'])}")
    add("")
    add("  by boundary spanning")
    for k, r in D.groupby('spans_boundary')['lost'].agg(['mean', 'size']).iterrows():
        add(f"    {str(k):>14s}  lost {r['mean'] * 100:5.1f}%   n={int(r['size'])}")
    add("")
    strat('svlen', [0, 100, 300, 1000, 10 ** 9],
          ['50-99bp', '100-299bp', '300-999bp', '>=1kb'], 'by SV length')
    strat('interval_len', [0, 100, 200, 500, 10 ** 9],
          ['<100bp', '100-199bp', '200-499bp', '>=500bp'], 'by containing interval length')
    strat('gc', [0, .35, .45, .55, 1.01],
          ['<35%', '35-45%', '45-55%', '>=55%'], 'by local GC content')
    if trf_tree is not None:
        strat('tr_frac', [0, .01, .5, 1.01],
              ['none', '<50%', '>=50%'], 'by tandem repeat overlap')
    strat('sv_density_100kb', [0, 5, 15, 40, 10 ** 9],
          ['<5', '5-14', '15-39', '>=40'], 'by local SV density per 100kb')

    # ---------------- multivariable
    add("=" * 78)
    add("  Multivariable logistic regression on P(lost), IRLS")
    add("=" * 78)
    # SV class and boundary spanning cannot both enter as independent terms:
    # an insertion occupies a single reference base and so can never span a
    # boundary, which makes spans_boundary a strict subset of is_del and
    # separates the fit (|beta| > 6 with standard errors in the hundreds).
    # Encode the three states against insertions as the reference level.
    feats = ['del_contained', 'del_spanning', 'log_svlen', 'log_interval',
             'gc', 'tr_frac', 'density']
    M = D.copy()
    is_del = (M['svtype'] == 'DEL')
    M['del_contained'] = (is_del & ~M['spans_boundary'].astype(bool)).astype(float)
    M['del_spanning'] = (is_del & M['spans_boundary'].astype(bool)).astype(float)
    M['log_svlen'] = np.log10(M['svlen'].clip(lower=1))
    M['log_interval'] = np.log10(M['interval_len'].fillna(M['interval_len'].median()).clip(lower=1))
    M['density'] = M['sv_density_100kb'].fillna(M['sv_density_100kb'].median())
    M['gc'] = M['gc'].fillna(M['gc'].median())
    M['tr_frac'] = M['tr_frac'].fillna(0.0)
    use = [f for f in feats if M[f].notna().all() and M[f].std() > 0]
    # drop any term with no outcome variation within a level, which would
    # separate the fit the same way
    keep = []
    for f in use:
        if M[f].nunique() == 2:
            g = M.groupby(f)['lost'].mean()
            if (g == 0).any() or (g == 1).any():
                add(f"  note: dropping '{f}' — no outcome variation in one level "
                    f"(would separate the fit)")
                continue
        keep.append(f)
    use = keep
    Xr = M[use].to_numpy(float)
    Xr = (Xr - Xr.mean(0)) / np.where(Xr.std(0) == 0, 1, Xr.std(0))
    X = np.column_stack([np.ones(len(Xr)), Xr])
    beta, se = logistic_irls(X, M['lost'].to_numpy(float))
    add("")
    add("  coefficients are per standard deviation, so they are comparable")
    add(f"  {'term':16s} {'beta':>9s} {'se':>8s} {'z':>7s}   {'odds/SD':>8s}")
    names = ['intercept'] + use
    for i, nm in enumerate(names):
        z = beta[i] / se[i] if se[i] and np.isfinite(se[i]) and se[i] > 0 else np.nan
        add(f"  {nm:16s} {beta[i]:9.3f} {se[i]:8.3f} {z:7.2f}   {np.exp(beta[i]):8.3f}")
    add("")
    add("  Mappability and segmental-duplication tracks are not present in this run")
    add("  directory, so they are not in the model. Tandem repeat overlap is the")
    add("  available proxy for hard-to-align sequence.")

    text = "\n".join(out)
    print(text)
    with open(f"{args.prefix}.attribution.txt", 'w') as fh:
        fh.write(text + "\n")
    D.drop(columns=['_b'], errors='ignore').to_csv(
        f"{args.prefix}.discordance.tsv", sep='\t', index=False)
    pd.DataFrame(decomp).to_csv(f"{args.prefix}.decomposition.tsv", sep='\t', index=False)
    print(f"\nwrote {args.prefix}.attribution.txt, "
          f"{args.prefix}.decomposition.tsv, {args.prefix}.discordance.tsv")
    return 0


if __name__ == '__main__':
    sys.exit(main())
