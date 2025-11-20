#!/bin/bash

# Example script showing different ways to run the pipeline
# This replaces the original analysis.sh

set -e

PROJECT_DIR=$(pwd)

echo "================================================"
echo "SV Calling Pipeline - Nextflow Version"
echo "================================================"
echo ""

# Option 1: Run with params file (RECOMMENDED)
echo "Option 1: Run with parameter file"
echo "-----------------------------------"
echo "nextflow run main.nf -params-file params.yaml -profile singularity"
echo ""

# Option 2: Run with command line parameters
echo "Option 2: Run with command line parameters"
echo "--------------------------------------------"
cat << 'EOF'
nextflow run main.nf \
  --fasta /path/to/reference.fasta \
  --illumina_wes_bam /path/to/wes.bam \
  --illumina_wgs_bam /path/to/wgs.bam \
  --pacbio_bam /path/to/pacbio.bam \
  --ont_bam /path/to/ont.bam \
  --benchmark_vcf /path/to/truth.vcf.gz \
  --high_confidence_targets /path/to/targets1.bed \
  --gene_panel_targets /path/to/targets2.bed \
  --wes_utr_targets /path/to/targets3.bed \
  -profile singularity
EOF
echo ""

# Option 3: Run specific technologies only
echo "Option 3: Run only specific technologies"
echo "-------------------------------------------"
cat << 'EOF'
# Edit params.yaml and comment out unwanted technologies:
# illumina_wes_bam: null  # Skip WES
# illumina_wgs_bam: '/path/to/wgs.bam'  # Run WGS

nextflow run main.nf -params-file params.yaml -profile singularity
EOF
echo ""

# Option 4: Resume a failed run
echo "Option 4: Resume a previous run"
echo "----------------------------------"
echo "nextflow run main.nf -params-file params.yaml -profile singularity -resume"
echo ""

# Option 5: Run on SLURM cluster
echo "Option 5: Run on SLURM cluster"
echo "--------------------------------"
cat << 'EOF'
nextflow run main.nf \
  -params-file params.yaml \
  -profile slurm,singularity \
  -resume
EOF
echo ""

# Option 6: Dry run to check workflow
echo "Option 6: Test workflow without execution (stub run)"
echo "-------------------------------------------------------"
echo "nextflow run main.nf -params-file params.yaml -stub-run"
echo ""

echo "================================================"
echo "Key Differences from Bash Pipeline:"
echo "================================================"
echo "1. No manual checkpoint tracking - use -resume instead"
echo "2. Automatic parallelization - all technologies run in parallel"
echo "3. Better error handling - automatic retries on failures"
echo "4. Execution reports - timeline, resource usage, DAG visualization"
echo "5. Scalability - same command works on laptop, HPC, or cloud"
echo ""

echo "================================================"
echo "Quick Start:"
echo "================================================"
echo "1. Edit params.yaml with your file paths"
echo "2. Ensure all BAM files have .bai indexes"
echo "3. Ensure reference FASTA has .fai index"
echo "4. Ensure VCF files have .tbi indexes"
echo "5. Run: nextflow run main.nf -params-file params.yaml -profile singularity"
echo ""

# Uncomment to actually run the pipeline
# nextflow run main.nf -params-file params.yaml -profile singularity -resume
