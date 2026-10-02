#!/usr/bin/env python3
"""Every value the manuscript and its supplementary tables report, for one assembly.

    manuscript_values.py --results RESULTS_DIR --assembly A --outdir DIR

RESULTS_DIR has the pipeline's published layout: a run's --outdir, or the copy a
task stages of the parts read here.

    posthoc/                    metrics.tsv, svtype_accounting.{tsv,log},
                                composition_standardisation.tsv, fidelity.*.tsv,
                                inversion_precision.tsv, decomposition/, stratified/,
                                sensitivity/, membership/, svanalyzer/
    statistics/tables/          truvari_metrics_{real_intervals,simulated_intervals,
                                simulated_intervals_raw}.tsv
    target_transition_evidence/ tables/target_transition_evidence.transitions.tsv,
                                simulations/tables/simulation_transition_evidence.mechanisms.tsv
    benchmarked_calls/          */*/*.homref_excluded.tsv

Writes
    DIR/<A>.manuscript_numbers.md      every number the text and the figure captions
                                       quote, with how it is derived
    DIR/supplementary_tables/<A>.<table>.tsv
                                       each supplementary table, with the rounding it
                                       is printed with; supplementary_tables.index.tsv
                                       maps the files to the table numbers

A section whose inputs are missing (for example without --sensitivity_benchmarks)
is skipped and says so. Percentile ranks of the observed EX+UTR value within the
simulated sets are the share of simulated values <= observed, as in the figures
and Supplementary Tables 2 and 5; rarefied ranks come from bootstrap_metrics.py
(ties count half). Differences are in percentage points.
"""
from __future__ import annotations

import argparse
import csv
import statistics
from collections import Counter, defaultdict
from pathlib import Path

WGS = ["Illumina_WGS Manta", "Illumina_WGS Delly", "ONT CuteSV", "ONT Sniffles", "PacBio CuteSV", "PacBio Pbsv"]
WES = ["Illumina_WES Manta", "Illumina_WES Delly"]
LONG = ["ONT CuteSV", "ONT Sniffles", "PacBio CuteSV", "PacBio Pbsv"]
LABEL = {"Illumina_WES Manta": "Illumina WES Manta", "Illumina_WES Delly": "Illumina WES Delly",
         "Illumina_WGS Manta": "Illumina WGS Manta", "Illumina_WGS Delly": "Illumina WGS Delly",
         "ONT CuteSV": "ONT cuteSV", "ONT Sniffles": "ONT Sniffles", "PacBio CuteSV": "PacBio cuteSV",
         "PacBio Pbsv": "PacBio pbsv"}
TECH = {"Illumina_WES": "Illumina WES", "Illumina_WGS": "Illumina WGS", "ONT": "ONT", "PacBio": "PacBio"}
CALLER = {"Manta": "Manta", "Delly": "Delly", "CuteSV": "cuteSV", "Sniffles": "Sniffles", "Pbsv": "pbsv"}
TARGET = {"high_confidence": "HCI", "gene_panel": "GP", "wes_utr": "EX+UTR"}
TARGETS = list(TARGET.items())
METRICS = ("precision", "recall", "f1")
SIZE_BINS = ("50-99", "100-299", "300-999", "1000-4999", "5000-")
MECH = {
    "hci_candidate_excluded_by_target": "Candidate excluded by target",
    "hci_candidate_retained_as_fp": "Candidate retained as FP",
    "hci_candidate_reassigned_to_other_truth": "Candidate paired with another truth record",
    "hci_candidate_absent_despite_target_overlap": "Candidate absent despite target overlap",
    "hci_candidate_not_contained_in_target": "Candidate overlaps but is not contained",
    "target_candidate_reassigned_from_other_truth": "Candidate taken from another truth record",
}
DIRECTION = {"HCI_TP_to_target_FN": "HCI-TP to target-FN", "HCI_FN_to_target_TP": "HCI-FN to target-TP"}
GRID = [("primary", "Primary (--refdist 500)"), ("refdist100", "--refdist 100"), ("refdist200", "--refdist 200"),
        ("refdist1000", "--refdist 1000"), ("pctsize0.5", "--pctsize 0.5"), ("pctsize0.9", "--pctsize 0.9"),
        ("pctseq0.7", "--pctseq 0.7"), ("containment", "Containment")]
FIDELITY_LABEL = {
    "intervals": ("Intervals", 0), "merged_components": ("Merged components", 0), "total_bp": ("Total bp (merged)", 0),
    "length_p10": ("Interval length, 10th percentile (bp)", 1), "length_median": ("Interval length, median (bp)", 1),
    "length_p90": ("Interval length, 90th percentile (bp)", 1), "length_mean": ("Interval length, mean (bp)", 1),
    "merged_length_median": ("Merged length, median (bp)", 1), "spacing_median": ("Spacing, median (bp)", 0),
    "spacing_p10": ("Spacing, 10th percentile (bp)", 0), "chromosomes": ("Chromosomes", 0), "gc_fraction": ("GC fraction", 3),
    "segdups_fraction_bp": ("Segmental duplications, fraction of bp", 4),
    "lowmappability_fraction_bp": ("Low mappability, fraction of bp", 4),
    "tandem_repeats_fraction_bp": ("Tandem repeats, fraction of bp", 4),
}
# Supplementary table each file belongs to, by assembly (GRCh37, GRCh38).
TABLE_INDEX = [
    ("observed_metrics", "1", "4", "Truvari metrics for the observed targets"),
    ("reference_percentiles", "2", "5", "Percentile ranks and KDE within the simulated sets"),
    ("reference_summary", "3", "6", "Simulated medians, differences and SD"),
    ("target_transition_mechanisms", "7", "", "HCI-TP to EX+UTR-FN transitions by mechanism"),
    ("simulation_transition_mechanisms", "8", "8", "Transition mechanisms in the simulated sets"),
    ("membership", "9", "9", "Containment against any overlap: truth records"),
    ("membership_metric_changes", "9", "9", "Containment against any overlap: metric changes"),
    ("breakend_scope", "10", "10", "Breakend-scope sensitivity"),
    ("target_transition_counts", "", "11", "HCI-to-EX+UTR losses and gains"),
    ("padding_recovery", "12", "12", "Recovery of losses after target padding"),
    ("svtype_accounting", "13", "13", "SV-type accounting from caller VCF to HCI benchmark"),
    ("fidelity", "14", "14", "Simulated sets against EX+UTR"),
    ("fidelity_chromosomes", "14", "14", "Per-chromosome allocation"),
    ("per_type_metrics", "15", "15", "Precision and recall by SV type"),
    ("inversion_precision", "15", "15", "Precision with and without inversions"),
    ("uncertainty", "16", "16", "Bootstrap intervals and rarefied ranks"),
    ("composition_standardisation", "17", "17", "Simulated reference standardised to EX+UTR composition"),
    ("decomposition_recall", "18", "18", "Decomposition of the recall change"),
    ("decomposition_precision", "18", "18", "Decomposition of the precision change"),
    ("svanalyzer", "19", "19", "Replication with SVanalyzer"),
    ("stratification", "20", "20", "Stratification after matching"),
    ("extension_recovery", "20", "20", "Recovery after candidate-side extension"),
    ("extension_recall", "20", "20", "EX+UTR recall under extension"),
    ("threshold_transitions", "21", "21", "Threshold grid: transitions"),
    ("threshold_metrics", "21", "21", "Threshold grid: metrics"),
]


# ----------------------------------------------------------------------------- helpers

def rows(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def stem(pipe):
    return pipe.replace(" ", "_")


def n(x):
    return f"{int(x):,}"


def f3(x):
    return f"{float(x):.3f}"


def pp(x):
    return f"{100 * float(x):+.2f}"


def term(x):
    """A decomposition component of metric_decomposition.py (HCI minus target, positive
    = decrease) as the manuscript reports it: target minus HCI. 0.0 - x keeps an exact
    zero positive, so it prints as +0.00."""
    return 0.0 - float(x)


def pctl(values, obs):
    # Simulated values come from an R table printed to 15 significant digits;
    # compare with a tolerance so an exact tie counts as <=, as in the figures.
    return 100 * sum(v <= obs + 1e-12 for v in values) / len(values)


def quantile(values, q):
    v = sorted(values)
    k = (len(v) - 1) * q
    lo = int(k)
    hi = min(lo + 1, len(v) - 1)
    return v[lo] + (v[hi] - v[lo]) * (k - lo)


def fmt(v, digits):
    v = float(v)
    return f"{v:,.0f}" if digits == 0 else f"{v:,.{digits}f}"


class Results:
    def __init__(self, root: Path, asm: str):
        self.root, self.asm = root, asm
        self.post = root / "posthoc"
        self.metrics = defaultdict(dict)
        for r in rows(self.post / "metrics.tsv"):
            self.metrics[(r["setting"], r["target_set"])][(r["pipeline"], r["target"])] = r
        self.prim = self.metrics[("primary", "real")]

    def setting(self, s):
        return self.metrics[(s, "real")]

    def settings(self):
        return sorted(s for s, ts in self.metrics if ts == "real" and s != "primary")

    def pipes(self, wgs_only=True):
        return [p for p in (WGS if wgs_only else WES + WGS) if (p, "high_confidence") in self.prim]

    def path(self, *parts):
        return self.root.joinpath(*parts)

    def has(self, *parts):
        return self.path(*parts).exists()

    def decomposition(self, pipe, kind):
        return rows(self.post / "decomposition" / f"{stem(pipe)}.{kind}.tsv")

    def real_transitions(self):
        tr = rows(self.path("target_transition_evidence", "tables", "target_transition_evidence.transitions.tsv"))
        return [r for r in tr if r["target"] == "EX+UTR"]

    def sim_summary(self):
        # Written by R with row names: the header has one field fewer than each row.
        out = {}
        with open(self.path("statistics", "tables", "truvari_metrics_simulated_intervals.tsv")) as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
            for row in reader:
                assert len(row) == len(header) + 1, (len(row), len(header))
                out[row[0]] = dict(zip(header, row[1:]))
        return out

    def recovery(self):
        out = defaultdict(lambda: [0, 0, 0])
        for r in rows(self.post / "sensitivity" / f"{self.asm}.recovery_summary.tsv"):
            x = out[r["setting"]]
            x[0] += int(r["primary_losses"])
            x[1] += int(r["truth_restored_tp"])
            x[2] += int(r["restored_with_same_candidate"])
        return out

    def truth_composition(self):
        # The truth side of the strata is the same for every pipeline; the first
        # WGS pipeline is used, as in Figure 5A-B.
        return [r for r in self.decomposition(self.pipes()[0], "strata") if r["side"] == "truth"]


# ----------------------------------------------------------------------------- the report

class Report:
    def __init__(self):
        self.lines = []

    def section(self, title):
        self.lines.append(f"\n## {title}\n")

    def __call__(self, text=""):
        self.lines.append(text)

    def missing(self, what):
        self.lines.append(f"- not available: {what}")

    def text(self):
        return "\n".join(self.lines) + "\n"


def report_metrics(d, out):
    out.section("Primary metrics (P / R / F1; TP-base, FN, TP-comp, FP; candidate records = TP-comp + FP)")
    prim, pipes = d.prim, d.pipes()
    for pipe in d.pipes(wgs_only=False):
        cells = []
        for t in TARGET:
            r = prim[(pipe, t)]
            cells.append(f"{TARGET[t]} {f3(r['precision'])}/{f3(r['recall'])}/{f3(r['f1'])} "
                         f"(TPb {r['TP-base']}, FN {r['FN']}, TPc {r['TP-comp']}, FP {r['FP']}, "
                         f"cand {int(r['TP-comp']) + int(r['FP'])})")
        out(f"- {LABEL[pipe]}: " + "; ".join(cells))
    out()
    for target in ("gene_panel", "wes_utr"):
        for m in METRICS:
            dif = {p: float(prim[(p, target)][m]) - float(prim[(p, "high_confidence")][m]) for p in pipes}
            out(f"- {TARGET[target]} minus HCI, {m}: mean {100 * statistics.mean(dif.values()):+.2f} pp; "
                + ", ".join(f"{LABEL[p]} {pp(v)}" for p, v in dif.items()))
        out(f"  (means over {len(pipes)} WGS pipelines)")
    submitted = [p for p in LONG + ["Illumina_WGS Manta"] if p in pipes]
    for target in ("gene_panel", "wes_utr"):
        dif = [float(prim[(p, target)]["f1"]) - float(prim[(p, "high_confidence")]["f1"]) for p in submitted]
        out(f"- {TARGET[target]} minus HCI F1, the {len(submitted)} WGS pipelines without Delly, mean "
            f"{100 * statistics.mean(dif):+.2f} pp")
    wes = [p for p in WES if (p, "high_confidence") in prim]
    if wes:
        out(f"- Illumina WES, highest HCI recall: {max(float(prim[(p, 'high_confidence')]['recall']) for p in wes):.4f} ("
            + ", ".join(f"{LABEL[p]} {float(prim[(p, 'high_confidence')]['recall']):.4f}" for p in wes) + ")")


def report_homref(d, out):
    out.section("Calls genotyped 0/0 removed before benchmarking")
    files = sorted(d.root.glob("benchmarked_calls/**/*.homref_excluded.tsv"))
    if not files:
        out.missing("benchmarked_calls/*.homref_excluded.tsv")
        return
    for f in files:
        for r in rows(f):
            out(f"- {TECH[r['technology']]} {CALLER[r['caller']]}: {n(r['homref_removed'])} records removed, "
                f"{n(r['homref_removed_pass'])} of them PASS; {n(r['records'])} records, {n(r['kept'])} kept")


def report_reference(d, out):
    out.section("EX+UTR within the simulated sets (percentile <= observed)")
    raw = d.path("statistics", "tables", "truvari_metrics_simulated_intervals_raw.tsv")
    if not raw.exists():
        out.missing(str(raw.relative_to(d.root)))
        return
    prim, pipes = d.prim, d.pipes()
    sims = defaultdict(lambda: defaultdict(list))
    base_cnt = defaultdict(list)
    for r in rows(raw):
        key = f"{r['tech']} {r['caller']}"
        for m in METRICS:
            sims[key][m].append(float(r[m]))
        base_cnt[key].append(int(r["base.cnt"]))
    out(f"- truth records per simulated set ({LABEL[pipes[0]]}): median {statistics.median(base_cnt[pipes[0]]):.0f}")
    low = high = mid = 0
    absdiff = defaultdict(list)
    highest_low = None
    for pipe in pipes:
        for m in METRICS:
            obs = float(prim[(pipe, "wes_utr")][m])
            v = sims[pipe][m]
            rank = pctl(v, obs)
            med = statistics.median(v)
            absdiff[m].append(abs(med - obs))
            hci = float(prim[(pipe, "high_confidence")][m])
            tag = "LOW" if rank < 5 else "HIGH" if rank > 95 else "mid"
            low += tag == "LOW"
            high += tag == "HIGH"
            mid += tag == "mid"
            if tag == "LOW" and (highest_low is None or rank > highest_low[0]):
                highest_low = (rank, pipe, m)
            out(f"- {LABEL[pipe]} {m}: EX+UTR {f3(obs)}, sim median {f3(med)} (n={len(v)}), "
                f"rank {rank:.1f} [{tag}]; HCI {f3(hci)} at rank {pctl(v, hci):.1f} of sims")
    out(f"\n- combinations: {low} below 5th, {high} above 95th, {mid} central, of {low + high + mid}")
    if highest_low:
        out(f"- highest rank among the low tail: {highest_low[0]:.1f} ({LABEL[highest_low[1]]} {highest_low[2]})")
    for m in METRICS:
        out(f"- mean |sim median - EX+UTR| {m}: {100 * statistics.mean(absdiff[m]):.2f} pp")


def report_uncertainty(d, out):
    out.section("Bootstrap CI and rarefied ranks (bootstrap_metrics.py)")
    for pipe in d.pipes():
        for r in d.decomposition(pipe, "uncertainty"):
            out(f"- {LABEL[pipe]} {r['metric']}: {f3(r['observed'])} "
                f"[{f3(r['bootstrap_ci_low'])}, {f3(r['bootstrap_ci_high'])}], "
                f"components {r['components']}; full rank {float(r['percentile_full']):.1f} (mid-rank), "
                f"rarefied median {f3(r['simulated_median_rarefied'])} "
                f"[{f3(r['simulated_p2.5_rarefied'])}, {f3(r['simulated_p97.5_rarefied'])}], "
                f"rarefied rank {float(r['percentile_rarefied']):.1f}")


def truth_shares(strata):
    """(records, DEL share, >=1 kb share) of a list of truth strata rows."""
    total = sum(int(r["n"]) for r in strata)
    dels = sum(int(r["n"]) for r in strata if r["svtype"] == "DEL")
    big = sum(int(r["n"]) for r in strata if r["size_bin"] in ("1000-4999", "5000-"))
    return total, dels / total, big / total


def report_composition(d, out):
    out.section("Composition of the truth records each target scores (truth strata; Figure 5A-B)")
    t = d.truth_composition()
    shares = {}
    for key, label in TARGETS:
        sub = [r for r in t if r["target_set"] == "real" and r["target"] == key]
        shares[key] = truth_shares(sub)
        total, dels, big = shares[key]
        out(f"- {label}: {n(total)} truth records, {100 * dels:.1f}% deletions, {100 * big:.1f}% at least 1 kb")
    sim = [r for r in t if r["target_set"] == "simulated"]
    if sim:
        per_set = defaultdict(list)
        for r in sim:
            per_set[r["target"]].append(r)
        total, dels, big = truth_shares(sim)
        set_shares = [truth_shares(v) for v in per_set.values()]
        out(f"- Simulated (pooled over {len(per_set)} sets): {n(total)} truth records, {100 * dels:.1f}% deletions "
            f"(2.5-97.5% across sets {100 * quantile([s[1] for s in set_shares], 0.025):.1f}-"
            f"{100 * quantile([s[1] for s in set_shares], 0.975):.1f}), {100 * big:.1f}% at least 1 kb "
            f"({100 * quantile([s[2] for s in set_shares], 0.025):.1f}-{100 * quantile([s[2] for s in set_shares], 0.975):.1f})")
        out(f"- deletion share, simulated minus EX+UTR: {100 * (dels - shares['wes_utr'][1]):.1f} pp; "
            f"simulated minus HCI: {100 * (dels - shares['high_confidence'][1]):.1f} pp")
    out()
    out("Recall by type in HCI (truth strata; Figure 5D):")
    for pipe in d.pipes():
        g = defaultdict(lambda: [0, 0])
        for r in d.decomposition(pipe, "strata"):
            if r["target_set"] == "real" and r["target"] == "high_confidence" and r["side"] == "truth":
                g[r["svtype"]][0] += int(r["n"])
                g[r["svtype"]][1] += int(r["hci_tp"])
        out(f"- {LABEL[pipe]}: DEL {g['DEL'][1] / g['DEL'][0]:.3f}, INS {g['INS'][1] / g['INS'][0]:.3f}")


def report_standardisation(d, out):
    out.section("Composition standardisation (gap = simulated median - target)")
    for r in rows(d.post / "composition_standardisation.tsv"):
        out(f"- {LABEL[r['pipeline']]} {r['metric']} by {r['stratification']}: target {f3(r['target_value'])}, "
            f"gap raw {pp(r['gap_raw'])} -> standardised {pp(r['gap_standardised'])}; "
            f"rank raw {float(r['target_percentile_raw']):.1f} -> {float(r['target_percentile_standardised']):.1f}")


def report_decomposition(d, out):
    out.section("Decomposition, real targets (terms in pp, target minus HCI; transition counts)")
    pipes = d.pipes()
    for pipe in pipes:
        for r in d.decomposition(pipe, "decomposition"):
            if r["target_set"] != "real" or r["target"] == "high_confidence":
                continue
            out(f"- {LABEL[pipe]} {TARGET[r['target']]}: n_truth {r['n_truth']}, n_cand {r['n_candidate']}; "
                f"recall {f3(r['recall_hci'])} -> comp {f3(r['recall_composition_only'])} -> {f3(r['recall_target'])} "
                f"(comp {pp(term(r['recall_composition_component']))}, trans {pp(term(r['recall_transition_component']))}; "
                f"losses {r['truth_losses']}, gains {r['truth_gains']}); "
                f"precision {f3(r['precision_hci'])} -> comp {f3(r['precision_composition_only'])} -> {f3(r['precision_target'])} "
                f"(comp {pp(term(r['precision_composition_component']))}, trans {pp(term(r['precision_transition_component']))}; "
                f"cand losses {r['candidate_losses']}, gains {r['candidate_gains']}); "
                f"cand loss truth excluded {r['candidate_loss_hci_truth_excluded_by_target']}, "
                f"FP in target {r.get('candidate_hci_fp', '')}")
    out.section("Decomposition, simulated sets (medians; terms in pp, target minus HCI)")
    for pipe in pipes:
        sim = [r for r in d.decomposition(pipe, "decomposition") if r["target_set"] == "simulated"]
        if not sim:
            continue
        med = lambda k: statistics.median(term(r[k]) for r in sim)
        out(f"- {LABEL[pipe]}: recall comp {pp(med('recall_composition_component'))}, "
            f"trans {pp(med('recall_transition_component'))}; precision comp "
            f"{pp(med('precision_composition_component'))}, trans {pp(med('precision_transition_component'))} "
            f"(n={len(sim)})")
    out.section("False-positive origin audit (EX+UTR): FPs that were HCI TPs whose truth was excluded")
    for pipe in pipes:
        r = next(r for r in d.decomposition(pipe, "decomposition") if r["target_set"] == "real" and r["target"] == "wes_utr")
        fp = int(d.prim[(pipe, "wes_utr")]["FP"])
        excl = int(r["candidate_loss_hci_truth_excluded_by_target"])
        other = {k: r[k] for k in r if k.startswith("candidate_loss_") and r[k] not in ("0", "")}
        out(f"- {LABEL[pipe]}: {excl} of {fp} EX+UTR FPs ({100 * excl / fp if fp else 0:.1f}%); "
            f"all candidate-loss mechanisms {other}")
    out.section("Per-type metrics from strata (truth side recall, candidate side precision)")
    for pipe in pipes:
        strata = d.decomposition(pipe, "strata")
        for t in ("high_confidence", "wes_utr"):
            agg = defaultdict(lambda: [0, 0])
            for r in strata:
                if r["target_set"] == "real" and r["target"] == t and r["svtype"] in ("DEL", "INS"):
                    a = agg[(r["side"], r["svtype"])]
                    a[0] += int(r["n"])
                    a[1] += int(r["target_tp"])
            cells = []
            for side, name in (("truth", "recall"), ("candidate", "precision")):
                for st in ("DEL", "INS"):
                    k, tp = agg[(side, st)]
                    cells.append(f"{st} {name} {tp / k:.3f} (n={k})" if k else f"{st} {name} NA")
            out(f"- {LABEL[pipe]} {TARGET[t]}: " + "; ".join(cells))
        other = Counter()
        for r in strata:
            if (r["target_set"] == "real" and r["target"] == "high_confidence" and r["side"] == "candidate"
                    and r["svtype"] not in ("DEL", "INS")):
                other[r["svtype"]] += int(r["n"])
        if other:
            out(f"  other scored candidate types in HCI: {dict(other)}")


def report_membership(d, out):
    out.section("Membership rule: containment against any overlap")
    pipes = d.pipes()
    cont = d.setting("containment")
    if not cont:
        out.missing("containment benchmarks (--sensitivity_containment)")
    else:
        for t in TARGET:
            k = (pipes[0], t)
            if k in cont:
                n_c, n_o = int(cont[k]["truth_denominator"]), int(d.prim[k]["truth_denominator"])
                out(f"- {TARGET[t]} truth: containment {n_c}, overlap {n_o} (+{n_o - n_c}, +{100 * (n_o - n_c) / n_c:.2f}%)")
        for pipe in pipes:
            cells = []
            for t in TARGET:
                c, o = cont.get((pipe, t)), d.prim[(pipe, t)]
                if c:
                    cells.append(f"{TARGET[t]} dP {pp(float(o['precision']) - float(c['precision']))}, "
                                 f"dR {pp(float(o['recall']) - float(c['recall']))}, dF1 {pp(float(o['f1']) - float(c['f1']))}")
            out(f"- {LABEL[pipe]} (overlap minus containment): " + "; ".join(cells))
        long_f1 = [abs(float(d.prim[(p, "wes_utr")]["f1"]) - float(cont[(p, "wes_utr")]["f1"])) for p in LONG
                   if (p, "wes_utr") in cont]
        if long_f1:
            out(f"- long-read pipelines, |EX+UTR F1 change|: {100 * min(long_f1):.2f} to {100 * max(long_f1):.2f} pp")
        hg = [abs(float(d.prim[(p, t)]["f1"]) - float(cont[(p, t)]["f1"])) for p in pipes
              for t in ("high_confidence", "gene_panel") if (p, t) in cont]
        if hg:
            out(f"- largest |F1 change| in HCI and GP: {100 * max(hg):.2f} pp")
    summary = d.post / "membership" / "membership.summary.tsv"
    if not summary.exists():
        out.missing("posthoc/membership/membership.summary.tsv")
        return
    s = {(r["count"], r["quantity"]): r["value"] for r in rows(summary)}
    out(f"- independent VCF-BED intersection, EX+UTR: containment {s[('independent_intersection', 'target_containment')]}, "
        f"overlap {s[('independent_intersection', 'target_overlap')]}; gained by overlap "
        f"{s[('independent_intersection', 'gained_by_overlap')]} ({s[('independent_intersection', 'gained_types')]}), "
        f"median length {float(s[('independent_intersection', 'gained_svlen_median')]):g} bp "
        f"(range {s[('independent_intersection', 'gained_svlen_min')]}-{s[('independent_intersection', 'gained_svlen_max')]})")
    if s.get(("simulated", "sets"), "0") != "0":
        out(f"- simulated sets ({s[('simulated', 'sets')]}), Truvari conventions: median containment "
            f"{float(s[('simulated', 'containment_median')]):g}, median overlap {float(s[('simulated', 'overlap_median')]):g}, "
            f"median increase {float(s[('simulated', 'increase_pct_median')]):.1f}% "
            f"(95% of sets {float(s[('simulated', 'increase_pct_p2.5')]):.1f}-{float(s[('simulated', 'increase_pct_p97.5')]):.1f}%); "
            f"EX+UTR {float(s[('truvari_conventions', 'target_increase_pct')]):.1f}%")


def report_transitions(d, out):
    out.section("Record tracing, real EX+UTR")
    if not d.has("target_transition_evidence", "tables", "target_transition_evidence.transitions.tsv"):
        out.missing("target_transition_evidence (--generate_transition_evidence)")
        return
    tr = d.real_transitions()
    by = Counter((r["pipeline"], r["direction"], r["mechanism"]) for r in tr)
    for k, v in sorted(by.items()):
        out(f"- {LABEL.get(k[0], k[0])} {k[1]} {k[2]}: {v}")
    losses = [r for r in tr if r["direction"] == "HCI_TP_to_target_FN"]
    gains = [r for r in tr if r["direction"] == "HCI_FN_to_target_TP"]
    dist = sorted(int(float(r["candidate_nearest_edge_distance"])) for r in losses
                  if r["mechanism"] == "hci_candidate_excluded_by_target")
    types = Counter(r["truth_svtype"] for r in losses)
    out(f"- losses {len(losses)} ({', '.join(f'{t} {k}' for t, k in sorted(types.items()))}), gains {len(gains)}")
    if dist:
        out(f"- excluded-candidate distance to nearest edge: min {dist[0]}, max {dist[-1]}, "
            f"median {statistics.median(dist):g}, p90 {quantile(dist, 0.9):.0f}, <=100 bp {sum(x <= 100 for x in dist)}, "
            f"<=500 bp {sum(x <= 500 for x in dist)}, of {len(dist)}")
    out(f"- losses without Delly: {sum(r['pipeline'] != 'Illumina_WGS Delly' for r in losses)}")

    out.section("Simulation audit")
    mech_file = d.path("target_transition_evidence", "simulations", "tables", "simulation_transition_evidence.mechanisms.tsv")
    if not mech_file.exists():
        out.missing(str(mech_file.relative_to(d.root)))
        return
    tot = Counter()
    mech = defaultdict(Counter)
    for r in rows(mech_file):
        tot[r["direction"]] += int(r["n"])
        mech[(r["pipeline"], r["direction"])][r["mechanism"]] += int(r["n"])
    out(f"- totals by direction: {dict(tot)}")
    all_loss = sum(c for (p, dr), m in mech.items() if dr == "HCI_TP_to_target_FN" for c in m.values())
    all_excl = sum(m["hci_candidate_excluded_by_target"] for (p, dr), m in mech.items() if dr == "HCI_TP_to_target_FN")
    if all_loss:
        out(f"- HCI-TP->target-FN {all_loss}, exact candidate excluded {all_excl} ({100 * all_excl / all_loss:.2f}%)")
    for (p, dr), m in sorted(mech.items()):
        k = sum(m.values())
        out(f"  - {LABEL.get(p, p)} {dr}: {k}; " + ", ".join(f"{a} {b} ({100 * b / k:.2f}%)" for a, b in m.most_common()))


def report_sensitivity(d, out):
    out.section("Padding and extension: fixed-record recovery")
    rec_file = d.post / "sensitivity" / f"{d.asm}.recovery_summary.tsv"
    if not rec_file.exists():
        out.missing(f"posthoc/sensitivity/{d.asm}.recovery_summary.tsv (--sensitivity_benchmarks)")
        return
    for r in rows(rec_file):
        out(f"- {r['setting']} {LABEL.get(r['pipeline'], r['pipeline'])}: losses {r['primary_losses']}, "
            f"restored {r['truth_restored_tp']}, same candidate {r['restored_with_same_candidate']}")
    for s, (lost, restored, same) in sorted(d.recovery().items()):
        out(f"- TOTAL {s}: restored {restored} of {lost} (same candidate {same})")
    pipes = d.pipes()
    out.section("Padding and extension: aggregate metrics (EX+UTR)")
    for s in [k for k in d.settings() if k.startswith(("pad", "extend"))]:
        for pipe in pipes:
            r = d.setting(s).get((pipe, "wes_utr"))
            if r:
                out(f"- {s} {LABEL[pipe]}: truth {r['truth_denominator']}, {f3(r['precision'])}/{f3(r['recall'])}/{f3(r['f1'])}")
    out.section("Threshold grid: EX+UTR transitions per setting")
    grid = defaultdict(Counter)
    for r in rows(d.post / "sensitivity" / f"{d.asm}.threshold_mechanisms.tsv"):
        grid[(r["setting"], r["direction"])][(r["pipeline"], r["mechanism"])] += int(r["n"])
    for (s, dr), c in sorted(grid.items()):
        mechs = Counter()
        for (p, m), k in c.items():
            mechs[m] += k
        out(f"- {s} {dr}: total {sum(c.values())}; " + ", ".join(f"{m} {k}" for m, k in mechs.most_common()))
    out.section("Threshold grid: metrics")
    for s in [k for k in d.settings() if k.startswith(("refdist", "pctsize", "pctseq"))]:
        for pipe in pipes:
            cells = []
            for t in TARGET:
                r = d.setting(s).get((pipe, t))
                if r:
                    cells.append(f"{TARGET[t]} {f3(r['precision'])}/{f3(r['recall'])}/{f3(r['f1'])}")
            out(f"- {s} {LABEL[pipe]}: " + "; ".join(cells))


def report_stratification(d, out):
    out.section("Post-matching stratification vs independently restricted")
    for pipe in d.pipes():
        h = d.prim[(pipe, "high_confidence")]
        for r in rows(d.post / "stratified" / f"{stem(pipe)}.stratified.tsv"):
            o = d.prim[(pipe, r["target"])]
            out(f"- {LABEL[pipe]} {TARGET[r['target']]}: stratified {f3(r['precision'])}/{f3(r['recall'])}/{f3(r['f1'])} "
                f"(truth {r['truth_denominator']}); restricted {f3(o['precision'])}/{f3(o['recall'])}/{f3(o['f1'])}; "
                f"HCI {f3(h['precision'])}/{f3(h['recall'])}/{f3(h['f1'])}")
    gaps = defaultdict(list)
    for pipe in d.pipes():
        dec = {r["target"]: r for r in d.decomposition(pipe, "decomposition") if r["target_set"] == "real"}
        for r in rows(d.post / "stratified" / f"{stem(pipe)}.stratified.tsv"):
            if r["target"] == "wes_utr":
                for m in ("precision", "recall"):
                    gaps[m].append((abs(float(r[m]) - float(dec["wes_utr"][f"{m}_composition_only"])), pipe))
    for m, v in gaps.items():
        g, pipe = max(v)
        out(f"- largest |stratified - composition-only| EX+UTR {m}: {100 * g:.3f} pp ({LABEL[pipe]})")


def report_fidelity(d, out):
    out.section("Simulation fidelity")
    summ = d.post / "fidelity.summary.tsv"
    if not summ.exists():
        out.missing("posthoc/fidelity.summary.tsv")
        return
    f = {r["feature"]: r for r in rows(summ)}
    for r in f.values():
        out(f"- {r['feature']}: real {r['real']}, sim median {r['simulated_median']} "
            f"[{r['simulated_p2.5']}, {r['simulated_p97.5']}], real rank {r['real_percentile']}")
    out()
    if "total_bp" in f:
        real, sim = float(f["total_bp"]["real"]), float(f["total_bp"]["simulated_median"])
        out(f"- simulated sets cover {100 * (1 - sim / real):.1f}% less sequence than EX+UTR (median)")
    if "spacing_median" in f:
        out(f"- median spacing: real {float(f['spacing_median']['real']) / 1000:.1f} kb, "
            f"simulated {float(f['spacing_median']['simulated_median']) / 1000:.1f} kb")
    for key, name in (("gc_fraction", "GC"), ("segdups_fraction_bp", "segmental duplications"),
                      ("tandem_repeats_fraction_bp", "tandem repeats"), ("lowmappability_fraction_bp", "low mappability")):
        if key in f:
            out(f"- {name}: real {100 * float(f[key]['real']):.1f}%, simulated {100 * float(f[key]['simulated_median']):.1f}%")
    chrom = rows(d.post / "fidelity.chromosomes.tsv")
    dropped = [r for r in chrom if int(r["rounding_difference"]) != 0]
    out(f"- rounding: {len(dropped)} chromosomes differ; net {sum(int(r['rounding_difference']) for r in chrom)} intervals; "
        f"chromosomes with <5 real intervals: {[r['chrom'] for r in chrom if 0 < int(r['real_intervals']) < 5]}; "
        f"bp in dropped chromosomes {sum(int(r['real_bp']) for r in chrom if int(r['allocated_per_decile']) == 0)}")
    for r in dropped:
        out(f"  - {r['chrom']}: real {r['real_intervals']}, expected {r['expected_simulated']}, sim median {r['simulated_median']}")


def report_accounting(d, out):
    out.section("SV-type accounting log")
    log = d.post / "svtype_accounting.log"
    out(log.read_text() if log.exists() else "- not available: posthoc/svtype_accounting.log")
    out.section("Breakend-retained sensitivity (F1 change)")
    bnd = d.setting("bnd_withbnd")
    if not bnd:
        out.missing("breakend-retained benchmarks")
        return
    for pipe in d.pipes() + [p for p in WES if (p, "high_confidence") in d.prim]:
        cells = []
        for t in TARGET:
            b, o = bnd.get((pipe, t)), d.prim.get((pipe, t))
            if b and o:
                cells.append(f"{TARGET[t]} dF1 {pp(float(b['f1']) - float(o['f1']))}")
        out(f"- {LABEL[pipe]}: " + "; ".join(cells))


def report_svanalyzer(d, out):
    out.section("SVanalyzer")
    sva = d.post / "svanalyzer" / "svanalyzer_summary.tsv"
    if not sva.exists():
        out.missing("posthoc/svanalyzer/svanalyzer_summary.tsv")
        return
    tot = Counter()
    for r in rows(sva):
        out(f"- {LABEL[r['pipeline']]}: HCI {f3(r['hci_precision'])}/{f3(r['hci_recall'])}, EX+UTR "
            f"{f3(r['target_precision'])}/{f3(r['target_recall'])}; TP->FN {r['hci_tp_to_target_fn']}, "
            f"all partners excluded {r['all_partners_excluded']}, FN->TP {r['hci_fn_to_target_tp']}")
        for k in ("hci_tp_to_target_fn", "all_partners_excluded", "hci_fn_to_target_tp"):
            tot[k] += int(r[k])
    out(f"- all pipelines: TP->FN {tot['hci_tp_to_target_fn']}, all partners excluded {tot['all_partners_excluded']}, "
        f"FN->TP {tot['hci_fn_to_target_tp']}")


def report_inversions(d, out):
    out.section("Precision without inversion false positives")
    inv = d.post / "inversion_precision.tsv"
    if not inv.exists():
        out.missing("posthoc/inversion_precision.tsv")
        return
    for r in rows(inv):
        out(f"- {LABEL[r['pipeline']]} {TARGET[r['target']]}: FP {r['fp']}, inversions {r['fp_inversions']}, "
            f"precision {r['precision']} -> {r['precision_without_inversions']} ({r['difference_pp']} pp)")


# ----------------------------------------------------------------------------- supplementary tables

def t_observed(d):
    real = rows(d.path("statistics", "tables", "truvari_metrics_real_intervals.tsv"))
    by = {(r["tech"], r["caller"], r["range"]): r for r in real}
    body = []
    for t, tl in TARGETS:
        for pipe in d.pipes(wgs_only=False):
            tech, caller = pipe.split(" ")
            r = by[(tech, caller, t)]
            gt = float(r["gt_concordance"]) if r["gt_concordance"] not in ("NA", "") else 0.0
            body.append([tl, TECH[tech], CALLER[caller], n(r["TP.base"]), n(r["FP"]), n(r["FN"]), f3(r["precision"]),
                         f3(r["recall"]), f3(r["f1"]), n(r["base.cnt"]), n(r["comp.cnt"]), f"{gt:.3f}"])
    return ["Target", "Technology", "Caller", "TP", "FP", "FN", "Precision", "Recall", "F1", "Truth N", "Called N",
            "GT conc."], body


def t_percentiles(d):
    s = d.sim_summary()
    body = []
    for t, tl in TARGETS:
        for pipe in d.pipes():
            tech, caller = pipe.split(" ")
            r = s[f"{tech}-{caller}-{t}"]
            body.append([f"{LABEL[pipe]} / {tl}", f"{float(r['precision_percentile']):.2f}",
                         f"{float(r['recall_percentile']):.2f}", f"{float(r['f1_percentile']):.2f}",
                         f3(r["precision_kde"]), f3(r["recall_kde"]), f3(r["f1_kde"])])
    return ["Pipeline / target", "Prec. %ile", "Rec. %ile", "F1 %ile", "Prec. KDE", "Rec. KDE", "F1 KDE"], body


def t_sim_summary(d):
    s = d.sim_summary()
    body = []
    for t, tl in TARGETS:
        for pipe in d.pipes():
            tech, caller = pipe.split(" ")
            r = s[f"{tech}-{caller}-{t}"]
            body.append([f"{LABEL[pipe]} / {tl}"] + [f3(r[k]) for k in (
                "precision_median", "recall_median", "f1_median", "precision_diff", "recall_diff", "f1_diff",
                "precision_sd", "recall_sd", "f1_sd")])
    return ["Pipeline / target", "Prec. med", "Rec. med", "F1 med", "dPrec.", "dRec.", "dF1", "Prec. SD", "Rec. SD",
            "F1 SD"], body


def t_transition_mechanisms(d):
    tr = d.real_transitions()
    c = Counter(r["pipeline"] for r in tr if r["direction"] == "HCI_TP_to_target_FN"
                and r["mechanism"] == "hci_candidate_excluded_by_target")
    tot = Counter(r["pipeline"] for r in tr if r["direction"] == "HCI_TP_to_target_FN")
    body = [[LABEL[p], c[p], f"{100 * c[p] / tot[p]:.1f}" if tot[p] else "--"] for p in d.pipes()]
    body.append(["All WGS pipelines", sum(c.values()),
                 f"{100 * sum(c.values()) / sum(tot.values()):.1f}" if tot else "--"])
    return ["Pipeline", "Candidate excluded, n", "Percent of transitions"], body


def t_transition_counts(d):
    tr = d.real_transitions()
    loss = Counter(r["pipeline"] for r in tr if r["direction"] == "HCI_TP_to_target_FN")
    gain = Counter(r["pipeline"] for r in tr if r["direction"] == "HCI_FN_to_target_TP")
    body = [[LABEL[p], loss[p], gain[p]] for p in d.pipes()]
    body.append(["All WGS pipelines", sum(loss.values()), sum(gain.values())])
    return ["Pipeline", "Losses, n", "Gains, n"], body


def t_simulation_mechanisms(d):
    st = rows(d.path("target_transition_evidence", "simulations", "tables", "simulation_transition_evidence.mechanisms.tsv"))
    agg = defaultdict(Counter)
    for r in st:
        agg[(r["pipeline"], r["direction"])][r["mechanism"]] += int(r["n"])
    body = []
    totals = defaultdict(Counter)
    for p in d.pipes():
        for direction in ("HCI_TP_to_target_FN", "HCI_FN_to_target_TP"):
            m = agg.get((p, direction))
            if not m:
                continue
            total = sum(m.values())
            for mech, k in m.most_common():
                body.append([LABEL[p], DIRECTION[direction], MECH[mech], n(k), f"{100 * k / total:.2f}"])
                totals[direction][mech] += k
    for direction, m in totals.items():
        total = sum(m.values())
        for mech, k in m.most_common():
            body.append(["All WGS pipelines", DIRECTION[direction], MECH[mech], n(k), f"{100 * k / total:.2f}"])
    return ["Pipeline", "Direction", "Mechanism", "n", "Percent"], body


def t_membership(d):
    cont, prim = d.setting("containment"), d.prim
    p0 = d.pipes()[0]
    body = []
    for t, tl in TARGETS:
        c, o = int(cont[(p0, t)]["truth_denominator"]), int(prim[(p0, t)]["truth_denominator"])
        body.append([tl, "Truvari truth records", n(c), n(o), f"+{o - c:,} (+{100 * (o - c) / c:.2f}%)"])
    summary = d.post / "membership" / "membership.summary.tsv"
    if summary.exists():
        s = {(r["count"], r["quantity"]): r["value"] for r in rows(summary)}
        c = int(s[("independent_intersection", "target_containment")])
        o = int(s[("independent_intersection", "target_overlap")])
        body.append(["EX+UTR", "Independent intersecting truth records", str(c), str(o),
                     f"+{o - c} (+{100 * (o - c) / c:.2f}%)"])
    body2 = []
    for p in d.pipes():
        for t, tl in TARGETS:
            c, o = cont.get((p, t)), prim[(p, t)]
            if c:
                body2.append([LABEL[p], tl] + [pp(float(o[k]) - float(c[k])) for k in ("precision", "recall", "f1")])
    return (["Target", "Quantity", "Containment", "Overlap", "Difference"], body,
            ["Pipeline", "Target", "Precision (pp)", "Recall (pp)", "F1 (pp)"], body2)


def t_breakend(d):
    bnd, prim = d.setting("bnd_withbnd"), d.prim
    body = []
    for p in d.pipes(wgs_only=False):
        for t, tl in TARGETS:
            b, o = bnd.get((p, t)), prim.get((p, t))
            if b and o:
                body.append([LABEL[p], tl, n(o["FP"]), n(b["FP"]), f3(o["precision"]), f3(b["precision"]),
                             f3(o["recall"]), f3(b["recall"]), f3(o["f1"]), f3(b["f1"])])
    return ["Pipeline", "Target", "FP excl.", "FP ret.", "Prec. excl.", "Prec. ret.", "Rec. excl.", "Rec. ret.",
            "F1 excl.", "F1 ret."], body


def t_padding(d):
    rec = d.recovery()
    body = []
    for pad in (0, 20, 50, 100, 200, 300, 500):
        lost, restored, same = rec[f"pad{pad}"]
        body.append([pad, lost, restored, same, lost - restored])
    return ["Padding (bp)", "Original losses", "Truth TP", "Same candidate", "Truth FN"], body


def t_svtype(d):
    acc = rows(d.post / "svtype_accounting.tsv")
    g = defaultdict(lambda: defaultdict(int))
    for r in acc:
        key = (r["side"], r["pipeline"], r["svtype"])
        g[key][r["stage"]] += int(r["n"])
        if r["scored_type"]:
            g[key]["_scored"] = r["scored_type"]
    body = []
    for side, p in [("truth", "")] + [("candidate", p) for p in d.pipes(wgs_only=False)]:
        pname = "Truth" if side == "truth" else LABEL[p]
        for t in sorted({k[2] for k in g if k[0] == side and k[1] == p}):
            s = g[(side, p, t)]
            scored = s.get("_scored", "--") if s["scored_type"] else ("excluded" if s["bnd_excluded"] else "--")
            size, size_mo = s["size_in_range_counted"], s["size_in_range_match_only"]
            inh, inh_mo = s["in_hci_counted"], s["in_hci_match_only"]
            tp, tp_mo = s["hci_tp_counted"], s["hci_tp_match_only"]
            body.append([pname, t, n(s["records"]), n(s["pass"]), scored,
                         n(size) + (f" + {n(size_mo)}" if size_mo else ""),
                         n(inh) + (f" + {n(inh_mo)}" if inh_mo else ""),
                         "--" if side == "truth" else n(tp) + (f" + {n(tp_mo)}" if tp_mo else ""),
                         "--" if side == "truth" else n(s["hci_fp_counted"])])
    return ["Call set", "Type", "Records", "PASS", "Scored as", "Size range", "In HCI", "HCI TP", "HCI FP"], body


def t_fidelity(d):
    body = []
    for r in rows(d.post / "fidelity.summary.tsv"):
        lab, dig = FIDELITY_LABEL[r["feature"]]
        body.append([lab, fmt(r["real"], dig), fmt(r["simulated_median"], dig),
                     f"{fmt(r['simulated_p2.5'], dig)}--{fmt(r['simulated_p97.5'], dig)}",
                     f"{float(r['real_percentile']):.1f}"])
    chrom = rows(d.post / "fidelity.chromosomes.tsv")

    def key(r):
        c = r["chrom"].replace("chr", "")
        return (0, int(c)) if c.isdigit() else (1, c)
    body2 = [[r["chrom"], n(r["real_intervals"]), r["allocated_per_decile"], n(r["expected_simulated"]),
              f"{int(r['rounding_difference']):+d}"] for r in sorted(chrom, key=key)]
    body2.append(["All", n(sum(int(r["real_intervals"]) for r in chrom)), "",
                  n(sum(int(r["expected_simulated"]) for r in chrom)),
                  f"{sum(int(r['rounding_difference']) for r in chrom):+d}"])
    return (["Property", "EX+UTR", "Simulated median", "Simulated 2.5--97.5%", "Real rank"], body,
            ["Chromosome", "Real", "Allocated", "Expected", "Difference"], body2)


def t_uncertainty(d):
    body = []
    for p in d.pipes():
        for r in d.decomposition(p, "uncertainty"):
            body.append([LABEL[p], r["metric"].replace("f1", "F1"), f3(r["observed"]),
                         f"{f3(r['bootstrap_ci_low'])}--{f3(r['bootstrap_ci_high'])}", n(r["components"]),
                         f"{float(r['percentile_full']):.1f}", f3(r["simulated_median_rarefied"]),
                         f"{f3(r['simulated_p2.5_rarefied'])}--{f3(r['simulated_p97.5_rarefied'])}",
                         f"{float(r['percentile_rarefied']):.1f}"])
    return ["Pipeline", "Metric", "Observed", "Bootstrap 95%", "Components", "Full rank", "Rarefied median",
            "Rarefied 2.5--97.5%", "Rarefied rank"], body


def t_standardisation(d):
    body = []
    for r in rows(d.post / "composition_standardisation.tsv"):
        body.append([LABEL[r["pipeline"]], r["metric"], r["stratification"].replace("type_size", "type and size"),
                     f3(r["target_value"]), f3(r["simulated_median_raw"]), f3(r["simulated_median_standardised"]),
                     pp(r["gap_raw"]), pp(r["gap_standardised"]), f"{float(r['target_percentile_raw']):.1f}",
                     f"{float(r['target_percentile_standardised']):.1f}"])
    return ["Pipeline", "Metric", "Strata", "EX+UTR", "Sim. median", "Standardised", "Gap", "Gap std.", "Rank",
            "Rank std."], body


def t_per_type(d):
    body = []
    for p in d.pipes():
        strata = d.decomposition(p, "strata")
        for t, tl in TARGETS:
            agg = defaultdict(lambda: [0, 0])
            for r in strata:
                if r["target_set"] == "real" and r["target"] == t and r["svtype"] in ("DEL", "INS"):
                    x = agg[(r["side"], r["svtype"])]
                    x[0] += int(r["n"])
                    x[1] += int(r["target_tp"])
            cells = []
            for side in ("truth", "candidate"):
                for st in ("DEL", "INS"):
                    k, tp = agg[(side, st)]
                    cells.append(f"{tp / k:.3f} ({k:,})" if k else "--")
            body.append([LABEL[p], tl] + cells)
    return ["Pipeline", "Target", "DEL recall", "INS recall", "DEL precision", "INS precision"], body


def t_inversions(d):
    inv = {(r["pipeline"], r["target"]): int(r["fp_inversions"]) for r in rows(d.post / "inversion_precision.tsv")}
    body = []
    for p in d.pipes(wgs_only=False):
        for t, tl in TARGETS:
            r = d.prim[(p, t)]
            tp, fp = int(r["TP-comp"]), int(r["FP"])
            k = inv[(p, t)]
            body.append([LABEL[p], tl, n(fp), n(k), f3(tp / (tp + fp)) if tp + fp else "--",
                         f3(tp / (tp + fp - k)) if tp + fp - k else "--"])
    return ["Pipeline", "Target", "FP", "Inversion FP", "Precision", "Without inversions"], body


def t_decomposition(d):
    body_r, body_p = [], []
    for p in d.pipes():
        dec = d.decomposition(p, "decomposition")
        real = {r["target"]: r for r in dec if r["target_set"] == "real"}
        sim = [r for r in dec if r["target_set"] == "simulated"]
        for t in ("gene_panel", "wes_utr"):
            r = real[t]
            body_r.append([LABEL[p], TARGET[t], n(r["n_truth"]), f3(r["recall_hci"]), f3(r["recall_composition_only"]),
                           f3(r["recall_target"]), pp(term(r["recall_composition_component"])),
                           pp(term(r["recall_transition_component"])), r["truth_losses"], r["truth_gains"]])
            body_p.append([LABEL[p], TARGET[t], n(r["n_candidate"]), f3(r["precision_hci"]),
                           f3(r["precision_composition_only"]), f3(r["precision_target"]),
                           pp(term(r["precision_composition_component"])), pp(term(r["precision_transition_component"])),
                           r["candidate_losses"], r["candidate_loss_hci_truth_excluded_by_target"]])
        if sim:
            med = lambda k: statistics.median(float(x[k]) for x in sim)
            body_r.append([LABEL[p], "Simulated (median)", f"{med('n_truth'):,.0f}", f3(med("recall_hci")),
                           f3(med("recall_composition_only")), f3(med("recall_target")),
                           pp(term(med("recall_composition_component"))), pp(term(med("recall_transition_component"))),
                           f"{med('truth_losses'):g}", f"{med('truth_gains'):g}"])
            body_p.append([LABEL[p], "Simulated (median)", f"{med('n_candidate'):,.0f}", f3(med("precision_hci")),
                           f3(med("precision_composition_only")), f3(med("precision_target")),
                           pp(term(med("precision_composition_component"))), pp(term(med("precision_transition_component"))),
                           f"{med('candidate_losses'):g}", f"{med('candidate_loss_hci_truth_excluded_by_target'):g}"])
    head = ["Pipeline", "Target", "Truth N", "HCI", "Composition only", "Target", "Composition term",
            "Transition term", "Losses", "Gains"]
    head_p = ["Pipeline", "Target", "Candidate N", "HCI", "Composition only", "Target", "Composition term",
              "Transition term", "Losses", "Truth excluded"]
    return head, body_r, head_p, body_p


def t_svanalyzer(d):
    body = [[LABEL[r["pipeline"]], f3(r["hci_precision"]), f3(r["hci_recall"]), f3(r["target_precision"]),
             f3(r["target_recall"]), r["hci_tp_to_target_fn"], r["all_partners_excluded"], r["hci_fn_to_target_tp"]]
            for r in rows(d.post / "svanalyzer" / "svanalyzer_summary.tsv")]
    return ["Pipeline", "HCI prec.", "HCI rec.", "EX+UTR prec.", "EX+UTR rec.", "TP to FN", "All partners excluded",
            "FN to TP"], body


def t_stratification(d):
    body = []
    for p in d.pipes():
        for r in rows(d.post / "stratified" / f"{stem(p)}.stratified.tsv"):
            o = d.prim[(p, r["target"])]
            body.append([LABEL[p], TARGET[r["target"]], f3(r["precision"]), f3(r["recall"]), f3(r["f1"]),
                         f3(o["precision"]), f3(o["recall"]), f3(o["f1"])])
    return ["Pipeline", "Target", "Strat. prec.", "Strat. rec.", "Strat. F1", "Restr. prec.", "Restr. rec.",
            "Restr. F1"], body


def t_extension(d):
    rec = d.recovery()
    body = []
    for dist in (20, 50, 100, 200, 500):
        lost, restored, same = rec[f"extend{dist}"]
        body.append([str(dist), str(lost), str(restored), str(same), str(rec[f"pad{dist}"][1])])
    body2 = []
    for p in d.pipes():
        dec = {r["target"]: r for r in d.decomposition(p, "decomposition") if r["target_set"] == "real"}
        cells = [f3(d.setting(f"extend{x}")[(p, "wes_utr")]["recall"]) for x in (20, 50, 100, 200, 500)]
        body2.append([LABEL[p], f3(d.prim[(p, "wes_utr")]["recall"])] + cells + [f3(dec["wes_utr"]["recall_composition_only"])])
    return (["Extension (bp)", "Original losses", "Restored", "Same candidate", "Restored by padding"], body,
            ["Pipeline", "No extension", "20 bp", "50 bp", "100 bp", "200 bp", "500 bp", "Composition only"], body2)


def t_threshold(d):
    grid = defaultdict(Counter)
    for r in rows(d.post / "sensitivity" / f"{d.asm}.threshold_mechanisms.tsv"):
        grid[r["setting"]][(r["direction"], r["mechanism"])] += int(r["n"])
    body = []
    for s, sl in GRID:
        c = grid[s]
        body.append([sl, sum(v for (dr, m), v in c.items() if dr == "HCI_TP_to_target_FN"),
                     c[("HCI_TP_to_target_FN", "hci_candidate_excluded_by_target")],
                     c[("HCI_TP_to_target_FN", "hci_candidate_not_contained_in_target")],
                     sum(v for (dr, m), v in c.items() if dr == "HCI_FN_to_target_TP")])
    body2 = []
    for s, sl in GRID:
        src = d.prim if s == "primary" else d.setting(s)
        for p in d.pipes():
            e, h = src.get((p, "wes_utr")), src.get((p, "high_confidence"))
            if e and h:
                body2.append([sl, LABEL[p], f3(e["precision"]), f3(e["recall"]), f3(e["f1"]), f3(h["precision"]),
                              f3(h["recall"]), f3(h["f1"])])
    return (["Setting", "HCI-TP to EX+UTR-FN", "Candidate excluded", "Not contained", "HCI-FN to EX+UTR-TP"], body,
            ["Setting", "Pipeline", "EX+UTR prec.", "EX+UTR rec.", "EX+UTR F1", "HCI prec.", "HCI rec.", "HCI F1"], body2)


def supplementary_tables(d):
    """(name, builder, files it needs) in TABLE_INDEX order; builders return one or two (header, body) pairs."""
    sens = d.post / "sensitivity" / f"{d.asm}.recovery_summary.tsv"
    stats = d.path("statistics", "tables", "truvari_metrics_real_intervals.tsv")
    sim_stats = d.path("statistics", "tables", "truvari_metrics_simulated_intervals.tsv")
    tr = d.path("target_transition_evidence", "tables", "target_transition_evidence.transitions.tsv")
    sim_tr = d.path("target_transition_evidence", "simulations", "tables", "simulation_transition_evidence.mechanisms.tsv")
    return [
        (("observed_metrics",), t_observed, [stats]),
        (("reference_percentiles",), t_percentiles, [sim_stats]),
        (("reference_summary",), t_sim_summary, [sim_stats]),
        (("target_transition_mechanisms",), t_transition_mechanisms, [tr]),
        (("target_transition_counts",), t_transition_counts, [tr]),
        (("simulation_transition_mechanisms",), t_simulation_mechanisms, [sim_tr]),
        (("membership", "membership_metric_changes"), t_membership, [sens]),
        (("breakend_scope",), t_breakend, []),
        (("padding_recovery",), t_padding, [sens]),
        (("svtype_accounting",), t_svtype, [d.post / "svtype_accounting.tsv"]),
        (("fidelity", "fidelity_chromosomes"), t_fidelity, [d.post / "fidelity.summary.tsv"]),
        (("per_type_metrics",), t_per_type, []),
        (("inversion_precision",), t_inversions, [d.post / "inversion_precision.tsv"]),
        (("uncertainty",), t_uncertainty, []),
        (("composition_standardisation",), t_standardisation, [d.post / "composition_standardisation.tsv"]),
        (("decomposition_recall", "decomposition_precision"), t_decomposition, []),
        (("svanalyzer",), t_svanalyzer, [d.post / "svanalyzer" / "svanalyzer_summary.tsv"]),
        (("stratification",), t_stratification, []),
        (("extension_recovery", "extension_recall"), t_extension, [sens]),
        (("threshold_transitions", "threshold_metrics"), t_threshold, [sens]),
    ]


def write_tsv(path, header, body):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(body)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--outdir", default=".", type=Path)
    args = parser.parse_args()

    d = Results(args.results, args.assembly)
    out = Report()
    out(f"# {args.assembly}: values quoted in the manuscript")
    out()
    out("Generated by manuscript_values.py from the pipeline's published outputs. WGS pipelines: "
        + ", ".join(LABEL[p] for p in d.pipes()) + ".")
    for section in (report_metrics, report_homref, report_reference, report_uncertainty, report_composition,
                    report_standardisation, report_decomposition, report_membership, report_transitions,
                    report_sensitivity, report_stratification, report_fidelity, report_accounting, report_svanalyzer,
                    report_inversions):
        section(d, out)
    (args.outdir / f"{args.assembly}.manuscript_numbers.md").write_text(out.text())

    tables = args.outdir / "supplementary_tables"
    tables.mkdir(parents=True, exist_ok=True)
    written = set()
    for names, builder, needs in supplementary_tables(d):
        lacking = [str(p.relative_to(d.root)) for p in needs if not p.exists()]
        if lacking:
            print(f"skipped {', '.join(names)}: missing {', '.join(lacking)}")
            continue
        parts = builder(d)
        for i, name in enumerate(names):
            write_tsv(tables / f"{args.assembly}.{name}.tsv", parts[2 * i], parts[2 * i + 1])
            written.add(name)
    column = 1 if args.assembly == "GRCh37" else 2
    write_tsv(tables / "supplementary_tables.index.tsv", ["file", "supplementary_table", "content"],
              [[f"{args.assembly}.{e[0]}.tsv", e[column] or "--", e[3]] for e in TABLE_INDEX if e[0] in written])
    print(f"{args.assembly}: wrote {args.assembly}.manuscript_numbers.md and {len(written)} supplementary tables")


if __name__ == "__main__":
    main()
