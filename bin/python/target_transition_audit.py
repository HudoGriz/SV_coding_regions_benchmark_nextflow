#!/usr/bin/env python3
"""Audit why Truvari truth variants change outcome after target restriction.

The script compares a broad HCI Truvari benchmark with one or more restricted
benchmarks produced from the same truth and candidate VCFs. Truvari writes the
same ``MatchId`` on paired records in ``tp-base`` and ``tp-comp``. Those IDs are
used to recover the exact candidate paired to every HCI true positive.

For every truth record that changes from HCI TP to restricted-target FN, the
script determines whether the original HCI candidate:

* does not overlap the restricted BED and was therefore filtered before match;
* overlaps the BED but was reassigned to a different truth record;
* overlaps the BED and is an unmatched false positive; or
* has another/unresolved output state.

Both real targets and Nextflow simulation targets are supported. Output is a
record-level TSV plus mechanism and boundary-distance summary tables.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import ctypes
import gc
import gzip
import hashlib
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path



SVTYPE_RE = re.compile(r"(?:^|;)SVTYPE=([^;]+)")
SVLEN_RE = re.compile(r"(?:^|;)SVLEN=(-?\d+)")
END_RE = re.compile(r"(?:^|;)END=(\d+)")
MATCH_RE = re.compile(r"(?:^|;)MatchId=([^;]+)")


@dataclass(frozen=True)
class Record:
    chrom: str
    pos: int
    record_id: str
    svtype: str
    svlen: int
    start: int
    end: int
    match_id: str
    allele_digest: str

    @property
    def key(self) -> tuple:
        # Position/type/length are not a unique record identity.  In particular,
        # the denser GRCh38 truth VCF contains distinct resolved alleles with the
        # same values for all of those fields.  Include a compact REF/ALT digest
        # so outcome comparisons cannot collapse multiallelic records.
        return (
            self.chrom, self.pos, self.record_id, self.svtype, self.svlen,
            self.allele_digest,
        )


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def parse_record(line: str) -> Record:
    # Do not split resolved REF/ALT alleles into new Python strings: individual
    # sequence-resolved SV alleles can span tens of kilobases, and thousands of
    # such temporary copies exceed small analysis-container memory limits.
    chrom, pos_text, record_id, remainder = line.split("\t", 3)
    tab_positions = []
    offset = 0
    for _ in range(5):  # REF, ALT, QUAL, FILTER, then the end of INFO
        offset = remainder.find("\t", offset)
        if offset < 0:
            raise ValueError("malformed VCF record")
        tab_positions.append(offset)
        offset += 1
    info = remainder[tab_positions[3] + 1:tab_positions[4]]
    ref_alt = remainder[:tab_positions[1]]
    allele_digest = hashlib.sha256(ref_alt.encode()).hexdigest()[:20]
    pos = int(pos_text)
    svtype_match = SVTYPE_RE.search(info)
    svtype = svtype_match.group(1) if svtype_match else "NA"
    svlen_match = SVLEN_RE.search(info)
    if svlen_match:
        svlen = abs(int(svlen_match.group(1)))
    else:
        svlen = 0
    end_match = END_RE.search(info)
    start = pos - 1
    if svtype == "INS":
        end = start + 1
    elif end_match:
        end = max(int(end_match.group(1)), start + 1)
    else:
        end = start + max(svlen, 1)
    match = MATCH_RE.search(info)
    return Record(
        chrom=chrom,
        pos=pos,
        record_id=record_id,
        svtype=svtype,
        svlen=svlen,
        start=start,
        end=end,
        match_id=match.group(1) if match else "",
        allele_digest=allele_digest,
    )


def read_vcf(path: Path) -> list[Record]:
    if not path.exists():
        raise FileNotFoundError(path)
    records = []
    with open_text(path) as handle:
        for line in handle:
            if not line.startswith("#"):
                records.append(parse_record(line))
    return records


def read_bed(path: Path) -> dict[str, tuple[list[int], list[int]]]:
    intervals = defaultdict(list)
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            start, end = int(fields[1]), int(fields[2])
            if end > start:
                intervals[fields[0]].append((start, end))
    index = {}
    for chrom, values in intervals.items():
        values.sort()
        merged = []
        for start, end in values:
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        index[chrom] = ([value[0] for value in merged], [value[1] for value in merged])
    return index


def overlap_intervals(trees, record: Record):
    indexed = trees.get(record.chrom)
    if indexed is None:
        return []
    starts, ends = indexed
    idx = bisect.bisect_right(ends, record.start)
    hits = []
    while idx < len(starts) and starts[idx] < record.end:
        if ends[idx] > record.start:
            hits.append((starts[idx], ends[idx]))
        idx += 1
    return hits


def bed_context(trees, record: Record) -> dict:
    hits = overlap_intervals(trees, record)
    if not hits:
        indexed = trees.get(record.chrom)
        nearest_distance = math.nan
        if indexed is not None:
            starts, ends = indexed
            idx = bisect.bisect_left(starts, record.end)
            candidates = []
            for candidate_idx in (idx - 1, idx):
                if 0 <= candidate_idx < len(starts):
                    distance = max(
                        starts[candidate_idx] - record.end,
                        record.start - ends[candidate_idx],
                        0,
                    )
                    candidates.append(distance)
            if candidates:
                nearest_distance = min(candidates)
        return {
            "overlaps": False,
            "same_interval": "",
            "interval_len": math.nan,
            "left_edge_distance": math.nan,
            "right_edge_distance": math.nan,
            "nearest_edge_distance": nearest_distance,
        }
    # Choose the hit with the greatest record overlap; ties prefer the shorter
    # interval and then genomic order.
    best = max(
        hits,
        key=lambda hit: (
            min(record.end, hit[1]) - max(record.start, hit[0]),
            -(hit[1] - hit[0]),
            -hit[0],
        ),
    )
    left = record.start - best[0]
    right = best[1] - record.end
    return {
        "overlaps": True,
        "same_interval": f"{record.chrom}:{best[0]}-{best[1]}",
        "interval_len": best[1] - best[0],
        "left_edge_distance": left,
        "right_edge_distance": right,
        "nearest_edge_distance": min(abs(left), abs(right)),
    }


def output_path(bench_dir: Path, suffix: str) -> Path:
    matches = sorted(bench_dir.glob(f"*.{suffix}.vcf.gz"))
    exact = bench_dir / f"{suffix}.vcf.gz"
    if exact.exists():
        matches.append(exact)
        matches = sorted(set(matches))
    if len(matches) != 1:
        raise RuntimeError(f"expected one *.{suffix}.vcf.gz in {bench_dir}, found {len(matches)}")
    return matches[0]


def read_benchmark(bench_dir: Path) -> dict:
    outputs = {
        name: read_vcf(output_path(bench_dir, name))
        for name in ("tp-base", "tp-comp", "fn", "fp")
    }
    return {
        "tp_base": {record.key: record for record in outputs["tp-base"]},
        "fn": {record.key: record for record in outputs["fn"]},
        "tp_comp": {record.key: record for record in outputs["tp-comp"]},
        "fp": {record.key: record for record in outputs["fp"]},
        "comp_by_match": {
            record.match_id: record
            for record in outputs["tp-comp"]
            if record.match_id
        },
        "base_by_match": {
            record.match_id: record
            for record in outputs["tp-base"]
            if record.match_id
        },
    }


def find_bench(results: Path, tech: str, caller: str, target: str) -> Path:
    base = results / "real_intervals" / tech / "truvari" / caller / target
    if list(base.glob("*.tp-base.vcf.gz")):
        return base
    dirs = [path for path in base.iterdir() if path.is_dir()] if base.exists() else []
    if len(dirs) != 1:
        raise RuntimeError(f"expected one benchmark directory under {base}, found {len(dirs)}")
    return dirs[0]


def simulation_bench(results: Path, tech: str, caller: str, simulation: str) -> Path:
    base = results / "simulations" / "benchmarks" / tech / caller
    # Simulation outputs are published directly in the caller directory.
    expected = base / f"{tech}-{caller}-{simulation}.tp-base.vcf.gz"
    if not expected.exists():
        raise FileNotFoundError(expected)
    return base


def read_simulation_benchmark(results: Path, tech: str, caller: str, simulation: str) -> dict:
    base = simulation_bench(results, tech, caller, simulation)
    prefix = f"{tech}-{caller}-{simulation}"
    outputs = {
        name: read_vcf(base / f"{prefix}.{name}.vcf.gz")
        for name in ("tp-base", "tp-comp", "fn", "fp")
    }
    return {
        "tp_base": {record.key: record for record in outputs["tp-base"]},
        "fn": {record.key: record for record in outputs["fn"]},
        "tp_comp": {record.key: record for record in outputs["tp-comp"]},
        "fp": {record.key: record for record in outputs["fp"]},
        "comp_by_match": {
            record.match_id: record for record in outputs["tp-comp"] if record.match_id
        },
        "base_by_match": {
            record.match_id: record for record in outputs["tp-base"] if record.match_id
        },
    }


def classify_transition(
    truth: Record,
    hci_candidate: Record | None,
    restricted: dict,
    bed,
) -> tuple[str, dict, dict]:
    truth_context = bed_context(bed, truth)
    if hci_candidate is None:
        return "hci_match_pair_missing", truth_context, {}
    candidate_context = bed_context(bed, hci_candidate)
    if not candidate_context["overlaps"]:
        mechanism = "hci_candidate_excluded_by_target"
    elif hci_candidate.key in restricted["tp_comp"]:
        target_candidate = restricted["tp_comp"][hci_candidate.key]
        paired_truth = restricted["base_by_match"].get(target_candidate.match_id)
        if paired_truth is not None and paired_truth.key != truth.key:
            mechanism = "hci_candidate_reassigned_to_other_truth"
        else:
            mechanism = "candidate_tp_pair_unresolved"
    elif hci_candidate.key in restricted["fp"]:
        mechanism = "hci_candidate_retained_as_fp"
    else:
        mechanism = "hci_candidate_absent_despite_target_overlap"
    return mechanism, truth_context, candidate_context


def audit_one(
    assembly: str,
    target_name: str,
    pipeline: str,
    hci: dict,
    restricted: dict,
    bed,
) -> list[dict]:
    rows = []
    for key, truth in hci["tp_base"].items():
        if key not in restricted["fn"]:
            continue
        candidate = hci["comp_by_match"].get(truth.match_id)
        mechanism, truth_context, candidate_context = classify_transition(
            truth, candidate, restricted, bed
        )
        candidate = candidate or Record("", 0, "", "", 0, 0, 0, "", "")
        same_interval = bool(
            truth_context.get("same_interval")
            and truth_context.get("same_interval") == candidate_context.get("same_interval")
        )
        rows.append({
            "assembly": assembly,
            "target": target_name,
            "pipeline": pipeline,
            "direction": "HCI_TP_to_target_FN",
            "mechanism": mechanism,
            "truth_chrom": truth.chrom,
            "truth_pos": truth.pos,
            "truth_id": truth.record_id,
            "truth_svtype": truth.svtype,
            "truth_svlen": truth.svlen,
            "truth_allele_digest": truth.allele_digest,
            "truth_start": truth.start,
            "truth_end": truth.end,
            "hci_match_id": truth.match_id,
            "candidate_chrom": candidate.chrom,
            "candidate_pos": candidate.pos,
            "candidate_id": candidate.record_id,
            "candidate_svtype": candidate.svtype,
            "candidate_svlen": candidate.svlen,
            "candidate_allele_digest": candidate.allele_digest,
            "candidate_start": candidate.start,
            "candidate_end": candidate.end,
            "start_distance": candidate.start - truth.start if candidate.chrom else math.nan,
            "end_distance": candidate.end - truth.end if candidate.chrom else math.nan,
            "size_ratio": (
                min(candidate.svlen, truth.svlen) / max(candidate.svlen, truth.svlen)
                if candidate.svlen and truth.svlen else math.nan
            ),
            "truth_overlaps_target": truth_context.get("overlaps", False),
            "candidate_overlaps_target": candidate_context.get("overlaps", False),
            "truth_target_interval": truth_context.get("same_interval", ""),
            "candidate_target_interval": candidate_context.get("same_interval", ""),
            "truth_candidate_same_target_interval": same_interval,
            "target_interval_len": truth_context.get("interval_len", math.nan),
            "truth_nearest_edge_distance": truth_context.get("nearest_edge_distance", math.nan),
            "candidate_nearest_edge_distance": candidate_context.get("nearest_edge_distance", math.nan),
        })
    # Also audit the reverse transition so the net recall change is fully
    # accounted for. The restricted candidate is linked through its own
    # MatchId and then located in the broad HCI candidate outputs.
    for key, truth in restricted["tp_base"].items():
        if key not in hci["fn"]:
            continue
        candidate = restricted["comp_by_match"].get(truth.match_id)
        truth_context = bed_context(bed, truth)
        candidate_context = bed_context(bed, candidate) if candidate else {}
        if candidate is None:
            mechanism = "target_match_pair_missing"
            candidate = Record("", 0, "", "", 0, 0, 0, "", "")
        elif candidate.key in hci["fp"]:
            mechanism = "target_candidate_was_hci_fp"
        elif candidate.key in hci["tp_comp"]:
            hci_candidate = hci["tp_comp"][candidate.key]
            paired_truth = hci["base_by_match"].get(hci_candidate.match_id)
            mechanism = (
                "target_candidate_reassigned_from_other_truth"
                if paired_truth is not None and paired_truth.key != truth.key
                else "target_candidate_hci_pair_unresolved"
            )
        else:
            mechanism = "target_candidate_absent_from_hci_outputs"
        same_interval = bool(
            truth_context.get("same_interval")
            and truth_context.get("same_interval") == candidate_context.get("same_interval")
        )
        rows.append({
            "assembly": assembly, "target": target_name, "pipeline": pipeline,
            "direction": "HCI_FN_to_target_TP", "mechanism": mechanism,
            "truth_chrom": truth.chrom, "truth_pos": truth.pos, "truth_id": truth.record_id,
            "truth_svtype": truth.svtype, "truth_svlen": truth.svlen,
            "truth_allele_digest": truth.allele_digest,
            "truth_start": truth.start, "truth_end": truth.end, "hci_match_id": "",
            "candidate_chrom": candidate.chrom, "candidate_pos": candidate.pos,
            "candidate_id": candidate.record_id, "candidate_svtype": candidate.svtype,
            "candidate_svlen": candidate.svlen, "candidate_start": candidate.start,
            "candidate_allele_digest": candidate.allele_digest,
            "candidate_end": candidate.end,
            "start_distance": candidate.start - truth.start if candidate.chrom else math.nan,
            "end_distance": candidate.end - truth.end if candidate.chrom else math.nan,
            "size_ratio": (
                min(candidate.svlen, truth.svlen) / max(candidate.svlen, truth.svlen)
                if candidate.svlen and truth.svlen else math.nan
            ),
            "truth_overlaps_target": truth_context.get("overlaps", False),
            "candidate_overlaps_target": candidate_context.get("overlaps", False),
            "truth_target_interval": truth_context.get("same_interval", ""),
            "candidate_target_interval": candidate_context.get("same_interval", ""),
            "truth_candidate_same_target_interval": same_interval,
            "target_interval_len": truth_context.get("interval_len", math.nan),
            "truth_nearest_edge_distance": truth_context.get("nearest_edge_distance", math.nan),
            "candidate_nearest_edge_distance": candidate_context.get("nearest_edge_distance", math.nan),
        })
    return rows


def write_tsv(path: str, rows: list[dict], fieldnames: list[str] | None = None):
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def mechanism_summary(rows: list[dict]) -> list[dict]:
    counts = Counter(
        (row["assembly"], row["target"], row["pipeline"], row["direction"], row["mechanism"])
        for row in rows
    )
    totals = Counter(
        (row["assembly"], row["target"], row["pipeline"])
        for row in rows
    )
    return [
        {
            "assembly": key[0], "target": key[1], "pipeline": key[2],
            "direction": key[3], "mechanism": key[4], "n": count,
            "percent": 100 * count / sum(
                value for other, value in counts.items() if other[:4] == key[:4]
            ),
        }
        for key, count in sorted(counts.items())
    ]


def edge_bin(value) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "unknown"
    if value < 0:
        return "spans/outside"
    if value < 25:
        return "0-24"
    if value < 50:
        return "25-49"
    if value < 100:
        return "50-99"
    if value < 200:
        return "100-199"
    if value < 500:
        return "200-499"
    return ">=500"


def boundary_summary(rows: list[dict]) -> list[dict]:
    counts = Counter(
        (
            row["assembly"], row["target"], row["pipeline"],
            edge_bin(row["truth_nearest_edge_distance"]), row["mechanism"],
        )
        for row in rows
    )
    return [
        {
            "assembly": key[0], "target": key[1], "pipeline": key[2],
            "truth_edge_bin": key[3], "mechanism": key[4], "n": count,
        }
        for key, count in sorted(counts.items())
    ]


def parse_pipeline(value: str) -> tuple[str, str]:
    try:
        return tuple(value.split(":", 1))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("pipeline must be TECH:CALLER") from exc


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, help="published pipeline results root")
    parser.add_argument("--hci-bench", type=Path, help="direct HCI Truvari output directory")
    parser.add_argument("--target-bench", type=Path, help="direct restricted-target Truvari output directory")
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--target-bed", type=Path, help="BED for a real restricted target")
    parser.add_argument("--target-name", default="EX+UTR")
    parser.add_argument("--target-dir", default="wes_utr", help="real_intervals target directory name")
    parser.add_argument("--simulations", action="store_true", help="audit simulation target BEDs and benchmarks")
    parser.add_argument("--simulation-start", type=int, default=1, help="first numbered simulation to audit")
    parser.add_argument("--simulation-limit", type=int, default=0, help="0 means all simulations")
    parser.add_argument("--pipeline", action="append", default=[], help="TECH:CALLER; repeatable")
    parser.add_argument("--prefix", required=True, type=Path)
    args = parser.parse_args()

    direct_mode = args.hci_bench is not None or args.target_bench is not None
    if direct_mode and (args.hci_bench is None or args.target_bench is None):
        parser.error("--hci-bench and --target-bench must be supplied together")
    if not direct_mode and args.results is None:
        parser.error("--results is required unless direct benchmark directories are supplied")

    pipelines = [parse_pipeline(value) for value in args.pipeline] or [
        ("Illumina_WGS", "Manta"),
        ("ONT", "CuteSV"),
        ("ONT", "Sniffles"),
        ("PacBio", "CuteSV"),
        ("PacBio", "Pbsv"),
    ]
    if direct_mode and len(pipelines) != 1:
        parser.error("direct benchmark mode requires exactly one --pipeline TECH:CALLER")
    hci_by_pipeline = {}
    for tech, caller in pipelines:
        hci_by_pipeline[(tech, caller)] = read_benchmark(
            args.hci_bench if direct_mode else find_bench(args.results, tech, caller, "high_confidence")
        )

    rows = []
    if args.target_bed:
        bed = read_bed(args.target_bed)
        for tech, caller in pipelines:
            restricted = read_benchmark(
                args.target_bench if direct_mode else find_bench(args.results, tech, caller, args.target_dir)
            )
            rows.extend(audit_one(
                args.assembly,
                args.target_name,
                f"{tech} {caller}",
                hci_by_pipeline[(tech, caller)],
                restricted,
                bed,
            ))

    if args.simulations:
        bed_dir = args.results / "simulations" / "target_regions"
        simulations = sorted(
            (path.stem for path in bed_dir.glob("simulation*.bed")),
            key=lambda value: int(value.removeprefix("simulation")),
        )
        simulations = [
            value for value in simulations
            if int(value.removeprefix("simulation")) >= args.simulation_start
        ]
        if args.simulation_limit:
            simulations = simulations[:args.simulation_limit]
        for index, simulation in enumerate(simulations, 1):
            bed = read_bed(bed_dir / f"{simulation}.bed")
            for tech, caller in pipelines:
                restricted = read_simulation_benchmark(args.results, tech, caller, simulation)
                rows.extend(audit_one(
                    args.assembly,
                    simulation,
                    f"{tech} {caller}",
                    hci_by_pipeline[(tech, caller)],
                    restricted,
                    bed,
                ))
                del restricted
            del bed
            gc.collect()
            # Resolved SV alleles can be very large. CPython's allocator may
            # retain freed arenas across hundreds of benchmarks, so return
            # them to glibc between simulations on Linux.
            try:
                ctypes.CDLL("libc.so.6").malloc_trim(0)
            except OSError:
                pass
            if index % 25 == 0:
                print(f"audited {index}/{len(simulations)} simulations", file=sys.stderr)

    args.prefix.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(f"{args.prefix}.transitions.tsv", rows)
    write_tsv(
        f"{args.prefix}.mechanisms.tsv", mechanism_summary(rows),
        ["assembly", "target", "pipeline", "direction", "mechanism", "n", "percent"],
    )
    write_tsv(
        f"{args.prefix}.boundary.tsv", boundary_summary(rows),
        ["assembly", "target", "pipeline", "truth_edge_bin", "mechanism", "n"],
    )
    metadata = {
        "assembly": args.assembly,
        "results": str(args.results.resolve()) if args.results else None,
        "hci_bench": str(args.hci_bench.resolve()) if args.hci_bench else None,
        "target_bench": str(args.target_bench.resolve()) if args.target_bench else None,
        "pipelines": [f"{tech}:{caller}" for tech, caller in pipelines],
        "real_target": str(args.target_bed.resolve()) if args.target_bed else None,
        "simulations": args.simulations,
        "simulation_start": args.simulation_start,
        "simulation_limit": args.simulation_limit,
        "transition_count": len(rows),
        "mechanism_counts": Counter(row["mechanism"] for row in rows).most_common(),
    }
    with open(f"{args.prefix}.metadata.json", "w") as handle:
        json.dump(metadata, handle, indent=2)
        handle.write("\n")
    print(json.dumps(metadata, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
