#!/usr/bin/env python3
"""Split the HCI-to-target change in recall and precision into its parts.

For a restricted target T, scored independently with the same truth and
candidate VCFs as the broad HCI benchmark, every record T scores was also
scored in HCI. That gives an exact, additive decomposition for each metric:

  recall:    R_HCI - R_T  =  (R_HCI - R_comp)  +  (R_comp - R_T)
  precision: P_HCI - P_T  =  (P_HCI - P_comp)  +  (P_comp - P_T)

R_comp is the recall T would have if every truth record it retains kept its
HCI label: the composition-only expectation. The first term is therefore the
effect of which truth records T retains. The second term is exactly the net
count of records whose label changed, divided by T's truth denominator:

  R_comp - R_T = (HCI-TP -> T-FN  minus  HCI-FN -> T-TP) / n_truth(T)

and each of those transitions is classified by mechanism. P_comp and the
candidate-side transitions are the same construction over T's candidate
records. The candidate side answers where a target's false positives come
from: an HCI true positive becomes a target false positive when the truth
record it matched in HCI is not in T.

Outputs (TSV):
  <prefix>.decomposition.tsv   one row per pipeline x target: denominators, the
                               three values of each metric, transition counts
                               by mechanism, and F1 at each stage
  <prefix>.strata.tsv          truth and candidate counts per type x size stratum,
                               with HCI and target labels, for per-type metrics
                               and composition standardisation
  <prefix>.candidate_transitions.tsv
                               record-level candidate-side transitions for the
                               real targets (the false-positive origin audit)

Truth-side transitions are classified exactly as target_transition_audit.py
classifies them, so the published audit and this decomposition agree.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from target_transition_audit import (  # noqa: E402
    Record, bed_context, classify_transition, find_bench, read_bed, read_benchmark,
    read_simulation_benchmark,
)

SIZE_BINS = [(50, 100), (100, 300), (300, 1000), (1000, 5000), (5000, math.inf)]
SIMULATION_RE = re.compile(r"^simulation(\d+)\.bed$")

TRUTH_LOSS_MECHANISMS = [
    "hci_candidate_excluded_by_target",
    "hci_candidate_reassigned_to_other_truth",
    "hci_candidate_retained_as_fp",
    "hci_candidate_absent_despite_target_overlap",
    "candidate_tp_pair_unresolved",
    "hci_match_pair_missing",
]
CANDIDATE_LOSS_MECHANISMS = [
    "hci_truth_excluded_by_target",
    "hci_truth_matched_to_other_candidate",
    "hci_truth_retained_as_fn",
    "hci_truth_pair_unresolved",
    "hci_match_pair_missing",
]


def size_bin(size: int) -> str:
    for low, high in SIZE_BINS:
        if low <= size < high:
            return f"{low}-{'' if math.isinf(high) else high - 1}"
    return "<50"


def scored_type(record: Record) -> str:
    """The type Truvari compares after --dup-to-ins."""
    return "INS" if record.svtype.startswith("DUP") else record.svtype


def record_size(record: Record) -> int:
    return record.svlen if record.svlen else max(record.end - record.start, 0)


def f1(precision: float, recall: float) -> float:
    return 2 * precision * recall / (precision + recall) if precision + recall else math.nan


def ratio(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else math.nan


def classify_candidate_loss(candidate: Record, hci: dict, target: dict, truth_keys: set) -> str:
    truth = hci["base_by_match"].get(candidate.match_id)
    if truth is None:
        return "hci_match_pair_missing"
    if truth.key not in truth_keys:
        return "hci_truth_excluded_by_target"
    if truth.key in target["tp_base"]:
        paired = target["comp_by_match"].get(target["tp_base"][truth.key].match_id)
        if paired is not None and paired.key != candidate.key:
            return "hci_truth_matched_to_other_candidate"
        return "hci_truth_pair_unresolved"
    if truth.key in target["fn"]:
        return "hci_truth_retained_as_fn"
    return "hci_truth_pair_unresolved"


def decompose(assembly: str, pipeline: str, target_name: str, target_set: str,
              hci: dict, target: dict, bed, candidate_rows: list | None) -> tuple[dict, list]:
    truth_keys = set(target["tp_base"]) | set(target["fn"])
    candidate_keys = set(target["tp_comp"]) | set(target["fp"])
    n_truth, n_candidate = len(truth_keys), len(candidate_keys)

    # Truth side.
    truth_hci_tp = sum(1 for key in truth_keys if key in hci["tp_base"])
    truth_absent = sum(1 for key in truth_keys if key not in hci["tp_base"] and key not in hci["fn"])
    losses = Counter()
    for key in target["fn"]:
        if key in hci["tp_base"]:
            truth = hci["tp_base"][key]
            candidate = hci["comp_by_match"].get(truth.match_id)
            mechanism, _, _ = classify_transition(truth, candidate, target, bed)
            losses[mechanism] += 1
    gains = sum(1 for key in target["tp_base"] if key in hci["fn"])

    # Candidate side.
    candidate_hci_tp = sum(1 for key in candidate_keys if key in hci["tp_comp"])
    candidate_absent = sum(1 for key in candidate_keys if key not in hci["tp_comp"] and key not in hci["fp"])
    candidate_losses = Counter()
    for key in target["fp"]:
        if key in hci["tp_comp"]:
            candidate = hci["tp_comp"][key]
            mechanism = classify_candidate_loss(candidate, hci, target, truth_keys)
            candidate_losses[mechanism] += 1
            if candidate_rows is not None:
                truth = hci["base_by_match"].get(candidate.match_id)
                truth_context = bed_context(bed, truth) if truth else {}
                candidate_context = bed_context(bed, candidate)
                candidate_rows.append({
                    "assembly": assembly, "target": target_name, "pipeline": pipeline,
                    "direction": "HCI_TP_to_target_FP", "mechanism": mechanism,
                    "candidate_chrom": candidate.chrom, "candidate_pos": candidate.pos,
                    "candidate_id": candidate.record_id, "candidate_svtype": candidate.svtype,
                    "candidate_svlen": candidate.svlen,
                    "candidate_allele_digest": candidate.allele_digest,
                    "hci_match_id": candidate.match_id,
                    "truth_chrom": truth.chrom if truth else "",
                    "truth_pos": truth.pos if truth else "",
                    "truth_id": truth.record_id if truth else "",
                    "truth_svtype": truth.svtype if truth else "",
                    "truth_svlen": truth.svlen if truth else "",
                    "truth_allele_digest": truth.allele_digest if truth else "",
                    "truth_overlaps_target": truth_context.get("overlaps", ""),
                    "truth_nearest_edge_distance": truth_context.get("nearest_edge_distance", ""),
                    "candidate_nearest_edge_distance": candidate_context.get("nearest_edge_distance", ""),
                })
    candidate_gains = sum(1 for key in target["tp_comp"] if key in hci["fp"])

    recall_hci = ratio(len(hci["tp_base"]), len(hci["tp_base"]) + len(hci["fn"]))
    precision_hci = ratio(len(hci["tp_comp"]), len(hci["tp_comp"]) + len(hci["fp"]))
    recall_comp = ratio(truth_hci_tp, n_truth)
    precision_comp = ratio(candidate_hci_tp, n_candidate)
    recall_target = ratio(len(target["tp_base"]), n_truth)
    precision_target = ratio(len(target["tp_comp"]), n_candidate)

    row = {
        "assembly": assembly, "pipeline": pipeline, "target_set": target_set, "target": target_name,
        "n_truth": n_truth, "n_candidate": n_candidate,
        "truth_absent_from_hci": truth_absent, "candidate_absent_from_hci": candidate_absent,
        "recall_hci": recall_hci, "recall_composition_only": recall_comp, "recall_target": recall_target,
        "recall_composition_component": recall_hci - recall_comp,
        "recall_transition_component": recall_comp - recall_target,
        "truth_losses": sum(losses.values()), "truth_gains": gains,
        "precision_hci": precision_hci, "precision_composition_only": precision_comp,
        "precision_target": precision_target,
        "precision_composition_component": precision_hci - precision_comp,
        "precision_transition_component": precision_comp - precision_target,
        "candidate_losses": sum(candidate_losses.values()), "candidate_gains": candidate_gains,
        "f1_hci": f1(precision_hci, recall_hci),
        "f1_composition_only": f1(precision_comp, recall_comp),
        "f1_target": f1(precision_target, recall_target),
    }
    for mechanism in TRUTH_LOSS_MECHANISMS:
        row[f"truth_loss_{mechanism}"] = losses.get(mechanism, 0)
    for mechanism in CANDIDATE_LOSS_MECHANISMS:
        row[f"candidate_loss_{mechanism}"] = candidate_losses.get(mechanism, 0)
    unexpected = (set(losses) - set(TRUTH_LOSS_MECHANISMS)) | (set(candidate_losses) - set(CANDIDATE_LOSS_MECHANISMS))
    if unexpected:
        raise RuntimeError(f"unclassified mechanisms: {sorted(unexpected)}")

    # Identities that must hold exactly.
    if n_truth and truth_absent == 0:
        assert truth_hci_tp - sum(losses.values()) + gains == len(target["tp_base"]), row
    if n_candidate and candidate_absent == 0:
        assert candidate_hci_tp - sum(candidate_losses.values()) + candidate_gains == len(target["tp_comp"]), row

    strata = Counter()
    for key in truth_keys:
        record = target["tp_base"].get(key) or target["fn"][key]
        stratum = (scored_type(record), size_bin(record_size(record)))
        strata[("truth",) + stratum + ("n",)] += 1
        strata[("truth",) + stratum + ("target_tp",)] += key in target["tp_base"]
        strata[("truth",) + stratum + ("hci_tp",)] += key in hci["tp_base"]
    for key in candidate_keys:
        record = target["tp_comp"].get(key) or target["fp"][key]
        stratum = (scored_type(record), size_bin(record_size(record)))
        strata[("candidate",) + stratum + ("n",)] += 1
        strata[("candidate",) + stratum + ("target_tp",)] += key in target["tp_comp"]
        strata[("candidate",) + stratum + ("hci_tp",)] += key in hci["tp_comp"]
    strata_rows = []
    for side, svtype, size in sorted({key[:3] for key in strata}):
        strata_rows.append({
            "assembly": assembly, "pipeline": pipeline, "target_set": target_set,
            "target": target_name, "side": side, "svtype": svtype, "size_bin": size,
            "n": strata[(side, svtype, size, "n")],
            "target_tp": strata[(side, svtype, size, "target_tp")],
            "hci_tp": strata[(side, svtype, size, "hci_tp")],
        })
    return row, strata_rows


def write(path: Path, rows: list[dict]):
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", required=True, type=Path, help="pipeline --outdir")
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--pipeline", required=True, help="TECH:CALLER")
    parser.add_argument("--target", action="append", default=[],
                        help="NAME=BED for a real target, e.g. wes_utr=/path/exutr.bed; repeatable")
    parser.add_argument("--simulations", action="store_true", help="also decompose every simulated set")
    parser.add_argument("--hci-target", default="high_confidence")
    parser.add_argument("--prefix", required=True, type=Path)
    args = parser.parse_args()

    technology, caller = args.pipeline.split(":")
    pipeline = f"{technology} {caller}"
    hci = read_benchmark(find_bench(args.results, technology, caller, args.hci_target))

    rows, strata_rows, candidate_rows = [], [], []
    # The HCI row itself anchors the per-type tables.
    hci_row, hci_strata = decompose(args.assembly, pipeline, args.hci_target, "real", hci, hci, {}, None)
    rows.append(hci_row)
    strata_rows += hci_strata

    for spec in args.target:
        name, bed_path = spec.split("=", 1)
        target = read_benchmark(find_bench(args.results, technology, caller, name))
        row, strata = decompose(args.assembly, pipeline, name, "real", hci, target,
                                read_bed(Path(bed_path)), candidate_rows)
        rows.append(row)
        strata_rows += strata

    if args.simulations:
        bed_dir = args.results / "simulations" / "target_regions"
        beds = sorted(
            (int(match.group(1)), path)
            for path in bed_dir.iterdir()
            if (match := SIMULATION_RE.match(path.name))
        )
        for number, bed_path in beds:
            simulation = f"simulation{number}"
            target = read_simulation_benchmark(args.results, technology, caller, simulation)
            row, strata = decompose(args.assembly, pipeline, simulation, "simulated", hci, target,
                                    read_bed(bed_path), None)
            rows.append(row)
            strata_rows += strata

    args.prefix.parent.mkdir(parents=True, exist_ok=True)
    write(Path(f"{args.prefix}.decomposition.tsv"), rows)
    write(Path(f"{args.prefix}.strata.tsv"), strata_rows)
    write(Path(f"{args.prefix}.candidate_transitions.tsv"), candidate_rows)
    print(f"{pipeline}: {len(rows)} targets decomposed")


if __name__ == "__main__":
    main()
