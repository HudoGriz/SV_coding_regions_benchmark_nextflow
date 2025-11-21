#!/usr/bin/env python3
"""
Simulate target regions similar to exome+UTR regions for benchmarking.

This script generates random target regions by sampling exons, CDS, and UTR features
from a GENCODE GTF file, creating simulated target sets that mimic real exome capture
regions.
"""

import argparse
import gzip
import random
import sys
from collections import defaultdict
from pathlib import Path


def parse_gtf(gtf_file):
    """
    Parse GTF file and extract relevant features (exon, CDS, UTR).
    
    Args:
        gtf_file: Path to GTF file (can be gzipped)
    
    Returns:
        dict: Dictionary of features by chromosome
    """
    features = defaultdict(list)
    
    # Determine if file is gzipped
    open_func = gzip.open if gtf_file.endswith('.gz') else open
    mode = 'rt' if gtf_file.endswith('.gz') else 'r'
    
    print(f"Parsing GTF file: {gtf_file}", file=sys.stderr)
    
    with open_func(gtf_file, mode) as f:
        for line in f:
            # Skip comments
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            chrom, source, feature_type, start, end, score, strand, frame, attributes = fields
            
            # Only keep protein-coding gene features
            if 'gene_type "protein_coding"' not in attributes:
                continue
            
            # Keep exon, CDS, and UTR features
            if feature_type in ['exon', 'CDS', 'UTR', 'five_prime_utr', 'three_prime_utr']:
                features[chrom].append({
                    'chrom': chrom,
                    'start': int(start),
                    'end': int(end),
                    'type': feature_type,
                    'strand': strand
                })
    
    print(f"Parsed {sum(len(v) for v in features.values())} features from {len(features)} chromosomes", 
          file=sys.stderr)
    
    return features


def merge_intervals(intervals):
    """
    Merge overlapping intervals.
    
    Args:
        intervals: List of (start, end) tuples
    
    Returns:
        List of merged (start, end) tuples
    """
    if not intervals:
        return []
    
    # Sort by start position
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged = [sorted_intervals[0]]
    
    for current in sorted_intervals[1:]:
        last = merged[-1]
        # Check for overlap
        if current[0] <= last[1]:
            # Merge intervals
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)
    
    return merged


def simulate_target_regions(features, num_regions=None, coverage_fraction=0.02, seed=None):
    """
    Simulate target regions by randomly sampling features.
    
    Args:
        features: Dictionary of features by chromosome
        num_regions: Target number of regions (if None, uses coverage_fraction)
        coverage_fraction: Fraction of total features to include
        seed: Random seed for reproducibility
    
    Returns:
        dict: Simulated regions by chromosome
    """
    if seed is not None:
        random.seed(seed)
    
    simulated_regions = defaultdict(list)
    
    for chrom, chrom_features in features.items():
        if not chrom_features:
            continue
        
        # Calculate how many features to sample
        if num_regions:
            n_sample = min(num_regions, len(chrom_features))
        else:
            n_sample = max(1, int(len(chrom_features) * coverage_fraction))
        
        # Randomly sample features
        sampled_features = random.sample(chrom_features, n_sample)
        
        # Extract intervals
        intervals = [(f['start'], f['end']) for f in sampled_features]
        
        # Merge overlapping intervals
        merged_intervals = merge_intervals(intervals)
        
        simulated_regions[chrom] = merged_intervals
    
    total_regions = sum(len(v) for v in simulated_regions.values())
    print(f"Generated {total_regions} simulated target regions across {len(simulated_regions)} chromosomes",
          file=sys.stderr)
    
    return simulated_regions


def write_bed(regions, output_file):
    """
    Write regions to BED file.
    
    Args:
        regions: Dictionary of regions by chromosome
        output_file: Output BED file path
    """
    with open(output_file, 'w') as f:
        for chrom in sorted(regions.keys()):
            for start, end in sorted(regions[chrom]):
                # BED format: chrom, start (0-based), end, name, score, strand
                f.write(f"{chrom}\t{start-1}\t{end}\n")
    
    print(f"Wrote simulated regions to: {output_file}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Simulate target regions from GENCODE GTF file'
    )
    parser.add_argument(
        '--gtf',
        required=True,
        help='Input GTF file (can be gzipped)'
    )
    parser.add_argument(
        '--output-prefix',
        required=True,
        help='Output file prefix for BED files'
    )
    parser.add_argument(
        '--num-simulations',
        type=int,
        default=100,
        help='Number of simulated target sets to generate (default: 100)'
    )
    parser.add_argument(
        '--coverage',
        type=float,
        default=0.02,
        help='Fraction of features to include in each simulation (default: 0.02)'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=None,
        help='Random seed for reproducibility'
    )
    
    args = parser.parse_args()
    
    # Parse GTF file
    features = parse_gtf(args.gtf)
    
    if not features:
        print("ERROR: No protein-coding features found in GTF file", file=sys.stderr)
        sys.exit(1)
    
    # Generate simulations
    print(f"Generating {args.num_simulations} simulated target sets...", file=sys.stderr)
    
    for i in range(1, args.num_simulations + 1):
        # Use different seed for each simulation
        seed = args.seed + i if args.seed is not None else None
        
        # Simulate regions
        regions = simulate_target_regions(
            features,
            coverage_fraction=args.coverage,
            seed=seed
        )
        
        # Write to BED file
        output_file = f"{args.output_prefix}_simulation_{i}.bed"
        write_bed(regions, output_file)
    
    print(f"Successfully generated {args.num_simulations} simulated target sets", file=sys.stderr)


if __name__ == '__main__':
    main()
