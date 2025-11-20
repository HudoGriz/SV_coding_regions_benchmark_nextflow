#!/bin/bash
set -euo pipefail

# Local testing script for SV calling pipeline
# This script creates minimal test data and runs the pipeline in stub mode

echo "=================================================="
echo "SV Pipeline Local Test"
echo "=================================================="

# Create test directories
echo "Creating test directories..."
mkdir -p test_data
mkdir -p singularity_images

# Generate synthetic test data
echo "Generating synthetic reference genome..."
echo ">chr22" > test_data/reference.fa
python3 -c "import random; random.seed(42); print(''.join(random.choices('ACGT', k=1000000)))" >> test_data/reference.fa

echo "Indexing reference..."
samtools faidx test_data/reference.fa

# Create synthetic BAM files
echo "Creating synthetic BAM file..."
cat > test_data/test.sam << EOF
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr22	LN:1000000
@RG	ID:test	SM:HG002	PL:ILLUMINA
EOF

# Add some aligned reads
for i in {1..20}; do
  pos=$((i * 50000))
  echo -e "read_${i}\t0\tchr22\t${pos}\t60\t100M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\t*" >> test_data/test.sam
done

echo "Converting SAM to sorted BAM..."
samtools view -Sb test_data/test.sam > test_data/test.bam
samtools sort test_data/test.bam -o test_data/test_sorted.bam
samtools index test_data/test_sorted.bam
rm test_data/test.sam test_data/test.bam

# Create synthetic truth VCF
echo "Creating synthetic truth VCF..."
cat > test_data/truth.vcf << 'EOF'
##fileformat=VCFv4.2
##contig=<ID=chr22,length=1000000>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG002
chr22	10000	sv1	N	<DEL>	60	PASS	SVTYPE=DEL;SVLEN=-500;END=10500	GT	0/1
chr22	50000	sv2	N	<DEL>	60	PASS	SVTYPE=DEL;SVLEN=-1000;END=51000	GT	0/1
chr22	200000	sv3	N	<INS>	60	PASS	SVTYPE=INS;SVLEN=800;END=200000	GT	0/1
chr22	500000	sv4	N	<DUP>	60	PASS	SVTYPE=DUP;SVLEN=2000;END=502000	GT	0/1
chr22	800000	sv5	N	<DEL>	60	PASS	SVTYPE=DEL;SVLEN=-300;END=800300	GT	0/1
EOF

bgzip -c test_data/truth.vcf > test_data/truth.vcf.gz
tabix -p vcf test_data/truth.vcf.gz

# Create target BED files
echo "Creating target BED files..."
echo -e "chr22\t0\t1000000" > test_data/high_confidence.bed

cat > test_data/gene_panel.bed << EOF
chr22	10000	100000
chr22	200000	300000
chr22	500000	600000
chr22	700000	900000
EOF

cat > test_data/wes_utr.bed << EOF
chr22	5000	15000
chr22	45000	55000
chr22	195000	205000
chr22	495000	505000
chr22	795000	805000
EOF

# Create test params file
echo "Creating test parameters..."
cat > test_params.yaml << EOF
# Test parameters for local testing
run_name: 'local_test'
outdir: 'test_results'

# Reference files
fasta: 'test_data/reference.fa'

# Input BAM files - using same BAM for all technologies for simplicity
illumina_wes_bam: 'test_data/test_sorted.bam'
illumina_wgs_bam: 'test_data/test_sorted.bam'

# Benchmarking files
benchmark_vcf: 'test_data/truth.vcf.gz'
high_confidence_targets: 'test_data/high_confidence.bed'
gene_panel_targets: 'test_data/gene_panel.bed'
wes_utr_targets: 'test_data/wes_utr.bed'

# Reduce resources for testing
max_cpus: 2
max_memory: '6.GB'
max_time: '10.m'
EOF

echo ""
echo "Test data created successfully!"
echo ""
echo "=================================================="
echo "Running pipeline in stub mode..."
echo "=================================================="

# Run in stub mode (doesn't execute actual processes, just tests workflow logic)
nextflow run main.nf \
  -params-file test_params.yaml \
  -profile test \
  -stub-run

echo ""
echo "=================================================="
echo "Test completed!"
echo "=================================================="
echo ""
echo "To run a full test (requires Docker/Singularity):"
echo "  nextflow run main.nf -params-file test_params.yaml -profile test"
echo ""
echo "To clean up:"
echo "  rm -rf test_data/ test_results/ test_params.yaml work/ .nextflow*"
echo ""
