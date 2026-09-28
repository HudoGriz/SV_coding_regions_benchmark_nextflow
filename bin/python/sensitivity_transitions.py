#!/usr/bin/env python3
"""Record-level audits of the sensitivity benchmarks.

Two questions, both answered from the published pipeline outputs:

1. Threshold settings (refdist<N>, pctsize<X>, pctseq<X>, containment).
   Each setting scored HCI and the boundary target the same way, so the
   HCI-to-target transition audit is rerun per setting, exactly as for the
   primary configuration. If target-boundary candidate exclusion drives the
   losses, their number should grow with refdist: a looser positional tolerance
   admits HCI pairs whose candidate lies further from the truth, and so further
   outside the target.

2. Boundary settings (extend<N>, pad<N>).
   Every primary HCI-TP to target-FN truth record, and the exact candidate it
   was paired with in HCI, is followed into each boundary benchmark: is the
   truth record now TP, and is it paired with that same candidate? Under
   extend<N> the truth set is unchanged, so recovery here translates directly
   into the target's recall.

Outputs:
  <prefix>.threshold_transitions.tsv  record-level audit rows with a setting column
  <prefix>.threshold_mechanisms.tsv   counts per setting x pipeline x direction x mechanism
  <prefix>.recovery.tsv               one row per original loss per boundary setting
  <prefix>.recovery_summary.tsv       counts per setting x pipeline
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from target_transition_audit import audit_one, find_bench, read_bed, read_benchmark  # noqa: E402


def split(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def pipelines_in(results: Path, setting: str) -> list[tuple[str, str]]:
    root = results / "sensitivity" / setting
    if not root.exists():
        raise FileNotFoundError(root)
    return sorted(
        (tech_dir.name, caller_dir.name)
        for tech_dir in root.iterdir() if tech_dir.is_dir()
        for caller_dir in tech_dir.iterdir() if caller_dir.is_dir()
    )


def setting_bench(results: Path, setting: str, technology: str, caller: str, target: str) -> Path:
    return results / "sensitivity" / setting / technology / caller / target


def row_key(row: dict, prefix: str) -> tuple:
    return (
        row[f"{prefix}_chrom"], int(row[f"{prefix}_pos"]), row[f"{prefix}_id"],
        row[f"{prefix}_svtype"], int(row[f"{prefix}_svlen"]), row[f"{prefix}_allele_digest"],
    )


def write(path: Path, rows: list[dict]):
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--target", default="wes_utr")
    parser.add_argument("--target-bed", required=True, type=Path)
    parser.add_argument("--target-label", default="EX+UTR")
    parser.add_argument("--hci-target", default="high_confidence")
    parser.add_argument("--threshold-settings", default="refdist100,refdist200,refdist1000,pctsize0.5,pctsize0.9,pctseq0.7,containment")
    parser.add_argument("--recovery-settings", default="extend20,extend50,extend100,extend200,extend500,pad0,pad20,pad50,pad100,pad200,pad300,pad500")
    parser.add_argument("--primary-transitions", type=Path,
                        help="default: <results>/target_transition_evidence/tables/target_transition_evidence.transitions.tsv")
    parser.add_argument("--prefix", required=True, type=Path)
    args = parser.parse_args()

    bed = read_bed(args.target_bed)
    threshold_rows = []

    # The primary configuration is the reference point of the threshold series.
    for technology, caller in pipelines_in(args.results, split(args.threshold_settings)[0]):
        hci = read_benchmark(find_bench(args.results, technology, caller, args.hci_target))
        target = read_benchmark(find_bench(args.results, technology, caller, args.target))
        for row in audit_one(args.assembly, args.target_label, f"{technology} {caller}", hci, target, bed):
            threshold_rows.append({"setting": "primary", **row})

    for setting in split(args.threshold_settings):
        for technology, caller in pipelines_in(args.results, setting):
            hci = read_benchmark(setting_bench(args.results, setting, technology, caller, args.hci_target))
            target = read_benchmark(setting_bench(args.results, setting, technology, caller, args.target))
            for row in audit_one(args.assembly, args.target_label, f"{technology} {caller}", hci, target, bed):
                # The audit tests candidate eligibility by overlap. Under full
                # containment a candidate can overlap the target and still be
                # ineligible, which the overlap test reports as "absent".
                if setting == "containment" and row["mechanism"] == "hci_candidate_absent_despite_target_overlap":
                    row["mechanism"] = "hci_candidate_not_contained_in_target"
                threshold_rows.append({"setting": setting, **row})

    mechanism_counts = Counter(
        (row["setting"], row["pipeline"], row["direction"], row["mechanism"]) for row in threshold_rows
    )
    mechanism_rows = [
        {"assembly": args.assembly, "target": args.target_label, "setting": key[0], "pipeline": key[1],
         "direction": key[2], "mechanism": key[3], "n": n}
        for key, n in sorted(mechanism_counts.items())
    ]

    primary_path = args.primary_transitions or (
        args.results / "target_transition_evidence" / "tables" / "target_transition_evidence.transitions.tsv"
    )
    with primary_path.open() as handle:
        losses = [row for row in csv.DictReader(handle, delimiter="\t")
                  if row["direction"] == "HCI_TP_to_target_FN" and row["assembly"] == args.assembly]

    recovery_rows = []
    for setting in split(args.recovery_settings):
        for technology, caller in pipelines_in(args.results, setting):
            pipeline = f"{technology} {caller}"
            pipeline_losses = [row for row in losses if row["pipeline"] == pipeline]
            if not pipeline_losses:
                continue
            bench = read_benchmark(setting_bench(args.results, setting, technology, caller, args.target))
            for row in pipeline_losses:
                truth_key, candidate_key = row_key(row, "truth"), row_key(row, "candidate")
                truth = bench["tp_base"].get(truth_key)
                candidate = bench["tp_comp"].get(candidate_key)
                truth_status = "TP" if truth else "FN" if truth_key in bench["fn"] else "absent"
                candidate_status = "TP" if candidate else "FP" if candidate_key in bench["fp"] else "absent"
                recovery_rows.append({
                    "assembly": args.assembly, "setting": setting, "pipeline": pipeline,
                    "truth_chrom": row["truth_chrom"], "truth_pos": row["truth_pos"],
                    "truth_id": row["truth_id"], "hci_match_id": row["hci_match_id"],
                    "primary_mechanism": row["mechanism"],
                    "candidate_nearest_edge_distance": row["candidate_nearest_edge_distance"],
                    "truth_status": truth_status,
                    "exact_hci_candidate_status": candidate_status,
                    "restored_with_same_candidate": bool(
                        truth and candidate and truth.match_id and truth.match_id == candidate.match_id
                    ),
                })

    summary = Counter()
    for row in recovery_rows:
        key = (row["setting"], row["pipeline"])
        summary[key + ("losses",)] += 1
        summary[key + ("truth_tp",)] += row["truth_status"] == "TP"
        summary[key + ("same_candidate",)] += row["restored_with_same_candidate"]
    summary_rows = [
        {"assembly": args.assembly, "setting": setting, "pipeline": pipeline,
         "primary_losses": summary[(setting, pipeline, "losses")],
         "truth_restored_tp": summary[(setting, pipeline, "truth_tp")],
         "restored_with_same_candidate": summary[(setting, pipeline, "same_candidate")]}
        for setting, pipeline in sorted({key[:2] for key in summary})
    ]

    args.prefix.parent.mkdir(parents=True, exist_ok=True)
    write(Path(f"{args.prefix}.threshold_transitions.tsv"), threshold_rows)
    write(Path(f"{args.prefix}.threshold_mechanisms.tsv"), mechanism_rows)
    write(Path(f"{args.prefix}.recovery.tsv"), recovery_rows)
    write(Path(f"{args.prefix}.recovery_summary.tsv"), summary_rows)
    print(f"{len(threshold_rows)} threshold transitions, {len(recovery_rows)} recovery rows")


if __name__ == "__main__":
    main()
