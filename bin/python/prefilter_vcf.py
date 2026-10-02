#!/usr/bin/env python3
"""Write the records a restricted benchmark would consider, for another comparator.

Selects records from a VCF exactly as the pipeline's Truvari benchmarks do
before matching (Truvari's own code): inside a BED by --bench-overlaps
membership, surviving filter_call (--passonly and the monomorphic-reference and
single-breakend checks), breakends excluded as with --bnddist -1, and at least
--min-size and at most --max-size long. Records are written unchanged apart
from the optional ID rewrite below.

Used to hand independently restricted truth and candidate VCFs to SVanalyzer,
so the comparator is the only thing that differs from the Truvari analysis.
SVanalyzer identifies variants by the ID column and numbers "." IDs in file
order, which differs between subsets; --stable-ids therefore replaces each ID
with one derived from the record itself (chrom, pos, REF/ALT digest), identical
in every subset. Runs inside the Truvari image.
"""
from __future__ import annotations

import argparse
import hashlib

import pysam
import truvari


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--overlap", type=int, default=1, help="--bench-overlaps value; 0 = containment")
    parser.add_argument("--min-size", type=int, default=50)
    parser.add_argument("--max-size", type=int, default=50000)
    parser.add_argument("--stable-ids", action="store_true")
    parser.add_argument("--id-prefix", default="", help="prepended to stable IDs, e.g. truth_ or cand_")
    parser.add_argument("--output", required=True, help="bgzipped VCF; indexed alongside")
    args = parser.parse_args()

    params = truvari.VariantParams(passonly=True, sizemin=args.min_size, sizefilt=args.min_size,
                                   sizemax=args.max_size, bnddist=-1)
    vcf = truvari.VariantFile(args.vcf, params=params)
    tree = truvari.build_region_tree(vcf, None, args.bed)
    truvari.merge_region_tree_overlaps(tree)

    kept = 0
    with pysam.VariantFile(args.output, "wz", header=vcf.header) as out:
        for entry in vcf.fetch_regions(tree, overlap=args.overlap):
            if entry.filter_call(base=True) or entry.is_bnd() or entry.filter_size(base=True):
                continue
            record = entry.get_record()
            if args.stable_ids:
                digest = hashlib.sha256(f"{record.ref}\t{','.join(record.alts)}".encode()).hexdigest()[:16]
                record.id = f"{args.id_prefix}{record.chrom}_{record.pos}_{digest}"
            out.write(record)
            kept += 1
    pysam.tabix_index(args.output, preset="vcf", force=True)
    print(f"{args.output}: {kept} records")


if __name__ == "__main__":
    main()
