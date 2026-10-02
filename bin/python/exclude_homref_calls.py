#!/usr/bin/env python3
"""Drop calls genotyped homozygous reference from a single-sample SV VCF.

A call whose genotype is fully called and carries no ALT allele (0/0, 0|0, or
haploid 0) is the caller stating that the sample does not have the variant.
The truth sets only count records that carry an ALT allele, so such a record
is not a call and must not be scored: unmatched it would be a false positive,
matched it would be a true positive for a variant the caller rejected.

Truvari's own --no-ref also drops calls with a missing genotype (./.), which
would remove every call of a caller that does not genotype (cuteSV here), so
this filter is applied before benchmarking instead. Calls with a missing or
partly missing genotype are kept. A VCF without a GT field is passed through.

Writes the filtered VCF (bgzipped, tabix-indexed), a one-row count table, and
the number of removed records to --removed-count for the pipeline to branch on.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import pysam


def is_homref(record) -> bool:
    if not record.samples:
        return False
    genotype = record.samples[0].get("GT")
    if not genotype:
        return False
    return all(allele == 0 for allele in genotype)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path, help="filtered .vcf.gz")
    parser.add_argument("--counts", required=True, type=Path)
    parser.add_argument("--removed-count", required=True, type=Path)
    parser.add_argument("--technology", required=True)
    parser.add_argument("--caller", required=True)
    args = parser.parse_args()

    with pysam.VariantFile(str(args.input)) as source:
        if len(source.header.samples) > 1:
            sys.exit(f"ERROR: {args.input} has {len(source.header.samples)} samples; expected one")
        header = source.header.copy()
        header.add_line(f"##exclude_homref_calls=removed records genotyped homozygous reference; "
                        f"input={args.input.name}")
        records = removed = removed_pass = 0
        with pysam.VariantFile(str(args.output), "wz", header=header) as sink:
            for record in source:
                records += 1
                if is_homref(record):
                    removed += 1
                    removed_pass += "PASS" in record.filter.keys() or not record.filter.keys()
                    continue
                sink.write(record)
    pysam.tabix_index(str(args.output), preset="vcf", force=True)

    with args.counts.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["technology", "caller", "records", "homref_removed", "homref_removed_pass", "kept"])
        writer.writerow([args.technology, args.caller, records, removed, removed_pass, records - removed])
    args.removed_count.write_text(f"{removed}\n")
    print(f"{args.technology} {args.caller}: {removed} of {records} records genotyped homozygous reference "
          f"({removed_pass} PASS) removed")


if __name__ == "__main__":
    main()
