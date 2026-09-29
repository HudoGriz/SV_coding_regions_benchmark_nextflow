#!/usr/bin/env python3
"""Remove SA-tag entries that name contigs absent from a BAM's header.

build_grch38_analysis_bams.sh restricts the GRCh38 BAMs to the contigs of the
analysis reference. A record kept on an analysis contig can still list, in its
SA tag, another part of the same read aligned to a contig that was dropped (a
decoy, HLA or alt contig). The SAM specification requires SA reference names to
be header references, and pbsv call aborts on them ("Could not find contig ...
in reference file"). This script removes exactly those entries. The record, its
flags and every other tag are kept; an SA tag left with no entries is deleted.

The BAM is processed one header contig at a time by parallel workers. Each
writes <chunk_dir>/<header index>.bam with the input header plus one @PG line;
the caller concatenates the chunks in header order (samtools cat), which keeps
the coordinate sort, and indexes the result. Every record of the input sits on
a header contig: the script stops if the index reports unplaced reads, and each
worker checks its record count against the index.

Outputs a per-contig count table and the number of removed entries per absent
contig.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import pysam

PG_ID = "strip_absent_sa_entries"


def output_header(header: dict, command_line: str) -> dict:
    header = {key: [dict(line) for line in value] if isinstance(value, list) else dict(value)
              for key, value in header.items()}
    programs = header.setdefault("PG", [])
    ids = {program["ID"] for program in programs}
    pg_id, suffix = PG_ID, 1
    while pg_id in ids:
        suffix += 1
        pg_id = f"{PG_ID}.{suffix}"
    program = {"ID": pg_id, "PN": PG_ID, "CL": command_line}
    if programs:
        program["PP"] = programs[-1]["ID"]
    programs.append(program)
    return header


def filter_contig(task):
    index, contig, expected, input_path, chunk_path, header, kept = task
    counts = Counter()
    dropped = Counter()
    with pysam.AlignmentFile(input_path, "rb") as source, \
            pysam.AlignmentFile(chunk_path, "wb", header=header) as sink:
        for record in source.fetch(contig):
            counts["records"] += 1
            if record.has_tag("SA"):
                counts["records_with_sa"] += 1
                entries = [entry for entry in record.get_tag("SA").split(";") if entry]
                retained = []
                for entry in entries:
                    reference = entry.rsplit(",", 5)[0]
                    if reference in kept:
                        retained.append(entry)
                    else:
                        dropped[reference] += 1
                if len(retained) < len(entries):
                    counts["sa_entries_removed"] += len(entries) - len(retained)
                    if retained:
                        record.set_tag("SA", ";".join(retained) + ";", value_type="Z")
                        counts["sa_tags_shortened"] += 1
                    else:
                        record.set_tag("SA", None)
                        counts["sa_tags_deleted"] += 1
            sink.write(record)
    if counts["records"] != expected:
        raise RuntimeError(f"{contig}: read {counts['records']} records, the index reports {expected}")
    return index, contig, counts, dropped


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, type=Path, help="coordinate-sorted, indexed BAM")
    parser.add_argument("--chunk-dir", required=True, type=Path, help="new directory for the per-contig BAMs")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--counts", required=True, type=Path, help="per-contig count table (TSV)")
    parser.add_argument("--dropped", required=True, type=Path, help="removed entries per absent contig (TSV)")
    args = parser.parse_args()

    if args.chunk_dir.exists():
        sys.exit(f"ERROR: refusing to write into existing {args.chunk_dir}")

    with pysam.AlignmentFile(args.input, "rb") as source:
        if source.nocoordinate:
            sys.exit(f"ERROR: {args.input} holds {source.nocoordinate} unplaced records; expected none")
        header = output_header(source.header.to_dict(), " ".join(sys.argv))
        contigs = list(source.references)
        expected = {stat.contig: stat.total for stat in source.get_index_statistics()}

    kept = frozenset(contigs)
    args.chunk_dir.mkdir(parents=True)
    width = len(str(len(contigs)))
    tasks = [
        (index, contig, expected.get(contig, 0), str(args.input),
         str(args.chunk_dir / f"{index:0{width}d}.bam"), header, kept)
        for index, contig in enumerate(contigs)
    ]
    # Largest contigs first, so the longest workers start early.
    tasks.sort(key=lambda task: -task[2])

    results = {}
    dropped = Counter()
    with Pool(args.workers) as pool:
        for index, contig, counts, contig_dropped in pool.imap_unordered(filter_contig, tasks):
            results[index] = (contig, counts)
            dropped.update(contig_dropped)

    fields = ["records", "records_with_sa", "sa_entries_removed", "sa_tags_shortened", "sa_tags_deleted"]
    totals = Counter()
    with args.counts.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", *fields])
        for index in range(len(contigs)):
            contig, counts = results[index]
            totals.update(counts)
            writer.writerow([contig, *(counts[field] for field in fields)])
        writer.writerow(["total", *(totals[field] for field in fields)])
    with args.dropped.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["absent_contig", "sa_entries_removed"])
        for contig, count in sorted(dropped.items(), key=lambda item: (-item[1], item[0])):
            writer.writerow([contig, count])

    print("\t".join(f"{field}={totals[field]}" for field in fields))


if __name__ == "__main__":
    main()
