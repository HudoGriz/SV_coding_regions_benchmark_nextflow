#!/usr/bin/env python3
"""Trace HCI-to-target outcome changes in an SVanalyzer benchmark.

Reads two `svanalyzer benchmark` runs of the same pipeline, one on HCI and one
on the target, each given truth and candidate VCFs pre-filtered independently
to that interval set (prefilter_vcf.py, with stable IDs). SVanalyzer matches
many-to-many: a truth variant is a TP when any candidate is within the
thresholds. The pairs and their distances are in <prefix>.distances.

For every truth variant that is TP in HCI and FN in the target, the script
lists the candidates that matched it in HCI and checks whether any of them is
in the target's candidate VCF. With many-to-many matching and identical
thresholds, a retained partner would still match, so the expected answer is
that every such partner was removed by the target filter. The informative
numbers are how many truth variants change state this way and how much of the
HCI-to-target recall change they account for, compared with the same pipeline
under Truvari.

Output: one summary row per pipeline, and a record-level table.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

# The report gives recall and precision as percentages, e.g.
# "Recall (DTP/(DTP+FN)): 46.09%".
REPORT = {
    "recall": re.compile(r"^Recall \(DTP/\(DTP\+FN\)\):\s*([0-9.]+)%", re.M),
    "precision": re.compile(r"^Precision \(PTP/\(PTP\+FP\)\):\s*([0-9.]+)%", re.M),
    "fn": re.compile(r"^Number of undetected true variants \(FN\):\s*(\d+)", re.M),
    "fp": re.compile(r"^Number of false positives \(FP\):\s*(\d+)", re.M),
}


def vcf_ids(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return [line.split("\t", 3)[2] for line in handle if not line.startswith("#")]


def matches(distance_file: Path, normshift: float, normsizediff: float, normdist: float) -> dict:
    """truth id -> set of candidate ids that pass every threshold."""
    by_truth = defaultdict(set)
    with distance_file.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            # The first line is a header that also starts with DIST.
            if fields[0] != "DIST" or fields[1] == "ID1":
                continue
            test_id, truth_id = fields[1], fields[2]
            d1, d2, d3 = float(fields[10]), float(fields[11]), float(fields[12])
            if d1 <= normshift and d2 <= normsizediff and d3 <= normdist:
                by_truth[truth_id].add(test_id)
    return by_truth


def report_values(path: Path) -> dict:
    text = path.read_text()
    values = {}
    for key, pattern in REPORT.items():
        match = pattern.search(text)
        if not match:
            values[key] = ""
        elif key in ("recall", "precision"):
            values[key] = round(float(match.group(1)) / 100, 6)
        else:
            values[key] = int(match.group(1))
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--pipeline", required=True)
    parser.add_argument("--hci-prefix", required=True, type=Path, help="svanalyzer --prefix of the HCI run")
    parser.add_argument("--target-prefix", required=True, type=Path)
    parser.add_argument("--target-truth", required=True, type=Path)
    parser.add_argument("--target-test", required=True, type=Path)
    parser.add_argument("--normshift", type=float, required=True)
    parser.add_argument("--normsizediff", type=float, required=True)
    parser.add_argument("--normdist", type=float, required=True)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--records", required=True, type=Path)
    args = parser.parse_args()

    hci_matches = matches(Path(f"{args.hci_prefix}.distances"), args.normshift, args.normsizediff, args.normdist)
    target_matches = matches(Path(f"{args.target_prefix}.distances"), args.normshift, args.normsizediff, args.normdist)
    target_truth = vcf_ids(args.target_truth)
    target_test = set(vcf_ids(args.target_test))

    records = []
    lost = gained = 0
    for truth_id in target_truth:
        hci_tp = bool(hci_matches.get(truth_id))
        target_tp = bool(target_matches.get(truth_id))
        if hci_tp and not target_tp:
            lost += 1
            partners = sorted(hci_matches[truth_id])
            retained = [p for p in partners if p in target_test]
            records.append({
                "assembly": args.assembly, "pipeline": args.pipeline, "truth_id": truth_id,
                "hci_partners": ";".join(partners), "partners_in_target": ";".join(retained),
                "mechanism": "all_hci_partners_excluded_by_target" if not retained else "partner_retained_but_unmatched",
            })
        elif target_tp and not hci_tp:
            gained += 1

    hci_report = report_values(Path(f"{args.hci_prefix}.report"))
    target_report = report_values(Path(f"{args.target_prefix}.report"))
    n_target = len(target_truth)
    hci_tp_in_target = sum(1 for t in target_truth if hci_matches.get(t))
    summary = {
        "assembly": args.assembly, "pipeline": args.pipeline,
        "hci_recall": hci_report["recall"], "hci_precision": hci_report["precision"],
        "target_recall": target_report["recall"], "target_precision": target_report["precision"],
        "target_truth": n_target,
        "target_recall_composition_only": hci_tp_in_target / n_target if n_target else "",
        "hci_tp_to_target_fn": lost, "hci_fn_to_target_tp": gained,
        "all_partners_excluded": sum(r["mechanism"] == "all_hci_partners_excluded_by_target" for r in records),
        "recall_transition_component": (lost - gained) / n_target if n_target else "",
        "normshift": args.normshift, "normsizediff": args.normsizediff, "normdist": args.normdist,
    }

    for path, rows in ((args.summary, [summary]), (args.records, records)):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as handle:
            fields = list(rows[0]) if rows else ["assembly", "pipeline", "truth_id", "hci_partners",
                                                 "partners_in_target", "mechanism"]
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    print(f"{args.pipeline}: {lost} HCI-TP to target-FN, {summary['all_partners_excluded']} with every partner excluded")


if __name__ == "__main__":
    main()
