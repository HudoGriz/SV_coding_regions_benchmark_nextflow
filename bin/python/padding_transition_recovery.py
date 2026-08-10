#!/usr/bin/env python3
"""Trace original target losses through a series of padded-target benchmarks."""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from target_transition_audit import read_benchmark


DIR_RE = re.compile(
    r"^(GRCh3[78])_(Illumina_WGS|ONT|PacBio)_(Manta|CuteSV|Sniffles|Pbsv)_pad(\d+)$"
)


def row_key(row, prefix):
    return (
        row[f"{prefix}_chrom"], int(row[f"{prefix}_pos"]), row[f"{prefix}_id"],
        row[f"{prefix}_svtype"], int(row[f"{prefix}_svlen"]),
        row[f"{prefix}_allele_digest"],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--padding-dir", required=True, type=Path)
    parser.add_argument("--transitions", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    originals = {}
    for path in args.transitions:
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["direction"] != "HCI_TP_to_target_FN":
                    continue
                originals.setdefault((row["assembly"], row["pipeline"]), []).append(row)

    output = []
    for directory in sorted(args.padding_dir.iterdir()):
        match = DIR_RE.match(directory.name)
        if not match or not (directory / "summary.json").exists():
            continue
        assembly, technology, caller, padding = match.groups()
        pipeline = f"{technology} {caller}"
        rows = originals.get((assembly, pipeline), [])
        if not rows:
            continue
        bench = read_benchmark(directory)
        for row in rows:
            truth_key = row_key(row, "truth")
            candidate_key = row_key(row, "candidate")
            truth_status = (
                "TP" if truth_key in bench["tp_base"] else
                "FN" if truth_key in bench["fn"] else "absent"
            )
            candidate_status = (
                "TP" if candidate_key in bench["tp_comp"] else
                "FP" if candidate_key in bench["fp"] else "absent"
            )
            output.append({
                "assembly": assembly, "pipeline": pipeline, "padding_bp": int(padding),
                "truth_chrom": row["truth_chrom"], "truth_pos": row["truth_pos"],
                "truth_id": row["truth_id"], "hci_match_id": row["hci_match_id"],
                "truth_status": truth_status, "exact_hci_candidate_status": candidate_status,
            })

    fields = [
        "assembly", "pipeline", "padding_bp", "truth_chrom", "truth_pos", "truth_id",
        "hci_match_id", "truth_status", "exact_hci_candidate_status",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output)


if __name__ == "__main__":
    main()
