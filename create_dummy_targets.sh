#!/bin/bash
# =====================================================
# Create Dummy Target BED Files
# =====================================================
# Use this if you don't have target BED files yet
# and just want to test SV calling (without benchmarking)

set -e

PROJECT_DIR="/home/45483vrhovsek/genome_in_bottle/SV_coding_regions_benchmark_nextflow"
TARGETS_DIR="${PROJECT_DIR}/data/targets"

echo "================================================="
echo "Creating dummy target BED files"
echo "================================================="

# Create targets directory if it doesn't exist
mkdir -p "${TARGETS_DIR}"

# Create dummy BED files
# These are minimal BED files with one entry covering chromosome 1
# They won't be useful for actual benchmarking, but allow the pipeline to run

echo "Creating: ${TARGETS_DIR}/dummy_high_confidence.bed"
cat > "${TARGETS_DIR}/dummy_high_confidence.bed" << EOF
chr1	1	249250621
EOF

echo "Creating: ${TARGETS_DIR}/dummy_gene_panel.bed"
cat > "${TARGETS_DIR}/dummy_gene_panel.bed" << EOF
chr1	1	249250621
EOF

echo "Creating: ${TARGETS_DIR}/dummy_wes_utr.bed"
cat > "${TARGETS_DIR}/dummy_wes_utr.bed" << EOF
chr1	1	249250621
EOF

echo ""
echo "================================================="
echo "Dummy BED files created successfully!"
echo "================================================="
echo ""
echo "Files created:"
echo "  - ${TARGETS_DIR}/dummy_high_confidence.bed"
echo "  - ${TARGETS_DIR}/dummy_gene_panel.bed"
echo "  - ${TARGETS_DIR}/dummy_wes_utr.bed"
echo ""
echo "To use these in params_minimal_example.yaml:"
echo "  high_confidence_targets: '${TARGETS_DIR}/dummy_high_confidence.bed'"
echo "  gene_panel_targets: '${TARGETS_DIR}/dummy_gene_panel.bed'"
echo "  wes_utr_targets: '${TARGETS_DIR}/dummy_wes_utr.bed'"
echo ""
echo "================================================="
echo "NOTE: These are dummy files for testing only!"
echo "For real benchmarking, use actual target BED files."
echo "================================================="
