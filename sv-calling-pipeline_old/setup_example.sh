#!/bin/bash

# Example script showing how to set up the pipeline directory structure
# This mirrors the original bash script structure

# Set your project directory
PROJECT_DIR=$(pwd)

echo "Setting up SV calling pipeline directory structure..."

# Create main directories
mkdir -p ${PROJECT_DIR}/singularity_images
mkdir -p ${PROJECT_DIR}/data/references
mkdir -p ${PROJECT_DIR}/data/Illumina_wes/bam
mkdir -p ${PROJECT_DIR}/data/Illumina_wgs/bam
mkdir -p ${PROJECT_DIR}/data/Pacbio/bam
mkdir -p ${PROJECT_DIR}/data/ONT/bam

echo "Directory structure created!"
echo ""
echo "Next steps:"
echo "1. Copy your Singularity images to: ${PROJECT_DIR}/singularity_images/"
echo "   - manta_latest.sif"
echo "   - cutesv_latest.sif"
echo "   - pbsv_latest.sif"
echo "   - sniffles_latest.sif"
echo "   - samtools_latest.sif"
echo "   - truvari_modded.sif"
echo ""
echo "2. Copy your reference files to: ${PROJECT_DIR}/data/references/"
echo "   - human_hs37d5.fasta (and .fai index)"
echo "   - HG002_SVs_Tier1_v0.6.vcf.gz (and .tbi index)"
echo "   - HG002_SVs_Tier1_v0.6.bed"
echo "   - Paediatric_disorders.HG002_SVs_Tier1.bed"
echo "   - exome_utr_gtf.HG002_SVs_Tier1.bed"
echo "   - human_hs37d5.trf.bed (optional, for Sniffles)"
echo ""
echo "3. Copy your BAM files (with .bai indexes) to appropriate directories"
echo ""
echo "4. Edit params.yaml with your actual file paths"
echo ""
echo "5. Run the pipeline:"
echo "   nextflow run main.nf -params-file params.yaml -profile singularity"
