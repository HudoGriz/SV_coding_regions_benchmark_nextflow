#!/bin/bash
# =====================================================
# Run SV Calling Pipeline - Production Version
# =====================================================
# This script runs the Nextflow pipeline on production data
# Replaces: scripts/run_all_analyses_and_bench.sh

set -e  # Exit on first error

# -----------------------------------------------------
# Configuration
# -----------------------------------------------------
PROJECT_DIR="/home/45483vrhovsek/genome_in_bottle/SV_coding_regions_benchmark_nextflow"
PARAMS_FILE="${PROJECT_DIR}/params_production.yaml"
PROFILE="singularity"  # or "docker" depending on your setup

# -----------------------------------------------------
# Validate Setup
# -----------------------------------------------------
echo "================================================="
echo "SV Calling Pipeline - Production Run"
echo "================================================="
echo "Project directory: ${PROJECT_DIR}"
echo "Parameters file:   ${PARAMS_FILE}"
echo "Profile:           ${PROFILE}"
echo "================================================="
echo ""

# Check if params file exists
if [ ! -f "$PARAMS_FILE" ]; then
    echo "ERROR: Parameters file not found: ${PARAMS_FILE}"
    echo "Please create it or adjust the path."
    exit 1
fi

# -----------------------------------------------------
# Run Nextflow Pipeline
# -----------------------------------------------------
echo "Starting Nextflow pipeline..."
echo ""

nextflow run ${PROJECT_DIR}/main.nf \
    -params-file ${PARAMS_FILE} \
    -profile ${PROFILE} \
    -resume \
    -with-report ${PROJECT_DIR}/benchmarking_run/pipeline_info/report.html \
    -with-timeline ${PROJECT_DIR}/benchmarking_run/pipeline_info/timeline.html \
    -with-trace ${PROJECT_DIR}/benchmarking_run/pipeline_info/trace.txt \
    -with-dag ${PROJECT_DIR}/benchmarking_run/pipeline_info/dag.html

# -----------------------------------------------------
# Check Exit Status
# -----------------------------------------------------
if [ $? -eq 0 ]; then
    echo ""
    echo "================================================="
    echo "Pipeline completed successfully!"
    echo "================================================="
    echo "Results directory: ${PROJECT_DIR}/benchmarking_run"
    echo ""
    echo "Key outputs:"
    echo "  - SV calls:     benchmarking_run/calls/"
    echo "  - Benchmarking: benchmarking_run/real_intervals/"
    echo "  - Reports:      benchmarking_run/pipeline_info/"
    echo "================================================="
else
    echo ""
    echo "================================================="
    echo "Pipeline failed! Check the log files."
    echo "================================================="
    exit 1
fi
