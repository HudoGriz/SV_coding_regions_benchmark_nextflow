#!/usr/bin/env python3
"""Account for every caller record between the raw VCF and the HCI benchmark.

For each pipeline, the raw call set is passed through the same filters Truvari
applied, in order, using Truvari's own code and the parameters the HCI
benchmark saved in its params.json:

  records              every record in the caller VCF
  pass                 survives filter_call (--passonly, monomorphic reference,
                       single-end breakends)
  scored_type          breakends removed (--bnddist -1); DUP relabelled INS
                       (--dup-to-ins); every other type kept as the caller wrote it
  size_in_range        survives filter_size (--sizefilt, --sizemax)
  in_hci               of those, overlapping HCI (--bench-overlaps membership)
  hci_tp / hci_fp      the HCI benchmark's TP-comp and FP records

Stages from size_in_range on are split by size class. A candidate between
--sizefilt and --sizemin (30-49 bp by default) may match a truth record, but
when it does not it is dropped rather than counted as a false positive; only
records of at least --sizemin count toward FP.

Rows are split by the caller's SVTYPE and by the type Truvari scores. Nothing
is converted or rewritten before Truvari: coordinates, lengths and alleles are
the caller's own. The script checks that every counted in-HCI record is a TP or
an FP, that no match-only record is an FP, and that every match-only TP is an
in-HCI match-only record. The truth VCF is accounted the same way with base-side filters.

Runs inside the Truvari image.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import truvari

RAW_VCF = {
    ("Illumina_WGS", "Manta"): "Illumina_WGS/Manta/*.diploid_sv.vcf.gz",
    ("Illumina_WGS", "Delly"): "Illumina_WGS/Delly/*.vcf.gz",
    ("Illumina_WES", "Manta"): "Illumina_WES/Manta/*.diploid_sv.vcf.gz",
    ("ONT", "CuteSV"): "ONT/CuteSV/*.vcf.gz",
    ("ONT", "Sniffles"): "ONT/Sniffles/*.vcf.gz",
    ("PacBio", "CuteSV"): "PacBio/CuteSV/*.vcf.gz",
    ("PacBio", "Pbsv"): "PacBio/PBSV/*.vcf.gz",
}


def single(paths: list[Path], what: str) -> Path:
    if len(paths) != 1:
        raise RuntimeError(f"expected one {what}, found {len(paths)}: {paths}")
    return paths[0]


def variant_params(params_json: Path) -> truvari.VariantParams:
    saved = json.loads(params_json.read_text())
    keep = {key: value for key, value in saved.items()
            if key in truvari.VariantParams.DEFAULTS and key != "reference"}
    return truvari.VariantParams(**keep)


def raw_type(entry) -> str:
    svtype = entry.info.get("SVTYPE")
    if isinstance(svtype, (list, tuple)):
        svtype = svtype[0]
    return str(svtype) if svtype else entry.var_type().name


def scored_type(entry, params) -> str | None:
    """The type Truvari compares, or None when the record is not scored."""
    if entry.is_bnd():
        return None if params.bnddist == -1 else "BND"
    kind = entry.var_type().name
    return "INS" if params.dup_to_ins and kind == "DUP" else kind


def account(vcf_path: Path, params, hci_bed: Path, overlap: int, base: bool, counts: Counter):
    vcf = truvari.VariantFile(str(vcf_path), params=params)
    for entry in vcf:
        kind = raw_type(entry)
        counts[(kind, "", "records")] += 1
        try:
            if entry.filter_call(base=base):
                continue
        except ValueError:
            counts[(kind, "", "multiallelic_rejected")] += 1
            continue
        counts[(kind, "", "pass")] += 1
        scored = scored_type(entry, params)
        if scored is None:
            counts[(kind, "", "bnd_excluded")] += 1
            continue
        counts[(kind, scored, "scored_type")] += 1
        if entry.filter_size(base=base):
            continue
        counts[(kind, scored, f"size_in_range_{size_class(entry, params)}")] += 1

    vcf = truvari.VariantFile(str(vcf_path), params=params)
    tree = truvari.build_region_tree(vcf, None, str(hci_bed))
    truvari.merge_region_tree_overlaps(tree)
    for entry in vcf.fetch_regions(tree, overlap=overlap):
        try:
            if entry.filter_call(base=base):
                continue
        except ValueError:
            continue
        scored = scored_type(entry, params)
        if scored is None or entry.filter_size(base=base):
            continue
        counts[(raw_type(entry), scored, f"in_hci_{size_class(entry, params)}")] += 1


def size_class(entry, params) -> str:
    return "counted" if entry.var_size() >= params.sizemin else "match_only"


def bench_outputs(bench_dir: Path, name: str, params, counts: Counter, stage: str):
    path = single(sorted(bench_dir.glob(f"*.{name}.vcf.gz")), f"{name} VCF in {bench_dir}")
    for entry in truvari.VariantFile(str(path), params=params):
        scored = scored_type(entry, params) or "BND"
        counts[(raw_type(entry), scored, f"{stage}_{size_class(entry, params)}")] += 1


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--hci-bed", required=True, type=Path)
    parser.add_argument("--truth-vcf", required=True, type=Path)
    parser.add_argument("--hci-target", default="high_confidence")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    rows = []
    checks = []
    truth_done = False
    for (technology, caller), pattern in RAW_VCF.items():
        matches = sorted((args.results / "sv_calls").glob(pattern))
        if not matches:
            continue
        raw_vcf = single(matches, f"{technology} {caller} VCF")
        bench_dir = args.results / "real_intervals" / technology / "truvari" / caller / args.hci_target
        params_json = single(sorted(bench_dir.glob("*/params.json")), f"params.json in {bench_dir}")
        params = variant_params(params_json)
        # The boolean build saved true/false, the numeric build a count.
        overlap = int(json.loads(params_json.read_text()).get("bench_overlaps", 0))

        counts = Counter()
        account(raw_vcf, params, args.hci_bed, overlap, base=False, counts=counts)
        bench_outputs(bench_dir, "tp-comp", params, counts, "hci_tp")
        bench_outputs(bench_dir, "fp", params, counts, "hci_fp")
        for (kind, scored, stage), n in sorted(counts.items()):
            rows.append({"assembly": args.assembly, "side": "candidate",
                         "pipeline": f"{technology} {caller}", "svtype": kind,
                         "scored_type": scored, "stage": stage, "n": n, "raw_vcf": str(raw_vcf)})
        total = Counter()
        for (_, _, stage), n in counts.items():
            total[stage] += n
        checks.append((
            f"{technology} {caller}",
            total["in_hci_counted"] == total["hci_tp_counted"] + total["hci_fp_counted"]
            and total["hci_fp_match_only"] == 0
            and total["hci_tp_match_only"] <= total["in_hci_match_only"],
            dict(total),
        ))

        if not truth_done:
            truth_counts = Counter()
            account(args.truth_vcf, params, args.hci_bed, overlap, base=True, counts=truth_counts)
            bench_outputs(bench_dir, "tp-base", params, truth_counts, "hci_tp")
            bench_outputs(bench_dir, "fn", params, truth_counts, "hci_fn")
            for (kind, scored, stage), n in sorted(truth_counts.items()):
                rows.append({"assembly": args.assembly, "side": "truth", "pipeline": "",
                             "svtype": kind, "scored_type": scored, "stage": stage, "n": n,
                             "raw_vcf": str(args.truth_vcf)})
            truth_done = True

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    failed = False
    for pipeline, ok, total in checks:
        print(f"{pipeline}: in HCI {total.get('in_hci_counted', 0)} counted + "
              f"{total.get('in_hci_match_only', 0)} match-only; HCI TP {total.get('hci_tp_counted', 0)}"
              f"+{total.get('hci_tp_match_only', 0)}, FP {total.get('hci_fp_counted', 0)} "
              f"[{'ok' if ok else 'MISMATCH'}]")
        failed |= not ok
    if failed:
        raise SystemExit("record accounting does not reconcile with the HCI benchmark outputs")


if __name__ == "__main__":
    main()
