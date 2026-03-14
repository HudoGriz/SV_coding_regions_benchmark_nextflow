#!/bin/bash
set -euo pipefail

# =============================================================================
# Prepare data and generate params YAML for SV benchmarking pipeline
# =============================================================================
#
# Wrapper script that:
#   1. Downloads reference data, BAM files, and truth sets
#   2. Generates a ready-to-use Nextflow params YAML file
#
# Prerequisites:
#   - Singularity/Apptainer installed (for container images and indexing)
#   - wget installed
#   - ~500 GB free disk space per genome build
#   - Stable internet connection
#
# Usage:
#   bash prepare.sh --genome GRCh37 --outdir /path/to/data
#   bash prepare.sh --genome GRCh38 --outdir /path/to/data
#   bash prepare.sh --genome all    --outdir /path/to/data
#
# Output structure:
#   /path/to/data/
#     singularity_images/    Shared containers (used by both builds)
#     GRCh37/
#       data/                GRCh37-specific BAMs and references
#       params_GRCh37.yaml   Generated params file
#     GRCh38/
#       data/                GRCh38-specific BAMs and references
#       params_GRCh38.yaml   Generated params file
#
# After completion, run the pipeline with:
#   nextflow run main.nf -params-file /path/to/data/GRCh37/params_GRCh37.yaml -profile singularity
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<EOF
Usage: bash prepare.sh --genome <GRCh37|GRCh38|all> --outdir <path> [options]

Required:
  --genome   Reference genome build: GRCh37, GRCh38, or all (both)
  --outdir   Base output directory (will be created if it doesn't exist)

Options:
  --skip-download      Skip download phase
  --skip-preprocessing Skip preprocessing phase
  --skip-json          Skip params YAML generation
  --generate-yaml-only Backward-compatible alias for: --skip-download --skip-preprocessing
  -h, --help        Show this help message

Examples:
  # Full preparation for GRCh37
  bash prepare.sh --genome GRCh37 --outdir /scratch/sv_benchmark

  # Both genome builds (data separated into GRCh37/ and GRCh38/ subdirs)
  bash prepare.sh --genome all --outdir /scratch/sv_benchmark

  # Only generate YAML (no downloads or preprocessing)
  bash prepare.sh --genome GRCh37 --outdir /scratch/sv_benchmark --skip-download --skip-preprocessing

  # Download only (skip preprocessing and params YAML)
  bash prepare.sh --genome GRCh38 --outdir /scratch/sv_benchmark --skip-preprocessing --skip-json

Output structure:
  /scratch/sv_benchmark/
    singularity_images/     # Shared containers
    GRCh37/
      data/                 # GRCh37 BAMs and references
      params_GRCh37.yaml    # Ready-to-use params file
    GRCh38/
      data/                 # GRCh38 BAMs and references  
      params_GRCh38.yaml    # Ready-to-use params file
EOF
    exit 1
}

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
GENOME=""
OUTDIR=""
SKIP_DOWNLOAD=false
SKIP_PREPROCESSING=false
SKIP_JSON=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome)             GENOME="$2";       shift 2 ;;
        --outdir)             OUTDIR="$2";       shift 2 ;;
        --skip-download)       SKIP_DOWNLOAD=true; shift ;;
        --skip-preprocessing|--skip-preprocesing) SKIP_PREPROCESSING=true; shift ;;
        --skip-json)           SKIP_JSON=true; shift ;;
        --generate-yaml-only)
            SKIP_DOWNLOAD=true
            SKIP_PREPROCESSING=true
            shift
            ;;
        -h|--help)            usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

if [[ -z "${GENOME}" || -z "${OUTDIR}" ]]; then
    echo "ERROR: --genome and --outdir are required."
    usage
fi

GENOME_LOWER=$(echo "${GENOME}" | tr '[:upper:]' '[:lower:]')
if [[ "${GENOME_LOWER}" != "grch37" && "${GENOME_LOWER}" != "grch38" && "${GENOME_LOWER}" != "all" ]]; then
    echo "ERROR: --genome must be GRCh37, GRCh38, or all (got: ${GENOME})"
    exit 1
fi

# ---------------------------------------------------------------------------
# Create base output directory
# ---------------------------------------------------------------------------
mkdir -p "${OUTDIR}"
OUTDIR="$(cd "${OUTDIR}" && pwd)"  # Resolve to absolute path

echo "Base output directory: ${OUTDIR}"
echo ""

# ---------------------------------------------------------------------------
# Run download + YAML generation for each requested genome build
# ---------------------------------------------------------------------------

run_build() {
    local build="$1"
    local build_dir="${OUTDIR}/${build}"
    local data_dir="${build_dir}/data"
    
    # Create genome-specific subdirectory
    mkdir -p "${build_dir}"

    echo ""
    echo "============================================================"
    echo "  Preparing ${build}"
    echo "============================================================"
    echo "Build directory: ${build_dir}"
    echo ""

    # Step 1: Download + preprocessing
    local prep_flags=()
    if [[ "${SKIP_DOWNLOAD}" == "true" ]]; then
        prep_flags+=("--skip-download")
    fi
    if [[ "${SKIP_PREPROCESSING}" == "true" ]]; then
        prep_flags+=("--skip-preprocessing")
    fi

    if [[ "${build}" == "GRCh37" ]]; then
        bash "${SCRIPT_DIR}/download_and_prep_GRCh37.sh" "${build_dir}" "${OUTDIR}/singularity_images" "${prep_flags[@]}"
    else
        bash "${SCRIPT_DIR}/download_and_prep_GRCh38.sh" "${build_dir}" "${OUTDIR}/singularity_images" "${prep_flags[@]}"
    fi

    # Step 2: Generate params YAML
    if [[ "${SKIP_JSON}" == "false" ]]; then
        local params_file="${build_dir}/params_${build}.yaml"
        bash "${SCRIPT_DIR}/generate_params.sh" \
            --genome "${build}" \
            --datadir "${data_dir}" \
            --outfile "${params_file}"

        echo ""
        echo "Params file written to: ${params_file}"
        echo "Run the pipeline with:"
        echo "  nextflow run main.nf -params-file ${params_file} -profile singularity"
        echo ""
    else
        echo "Skipping params YAML generation (--skip-json)."
    fi
}

# Dispatch
if [[ "${GENOME_LOWER}" == "grch37" ]]; then
    run_build "GRCh37"
elif [[ "${GENOME_LOWER}" == "grch38" ]]; then
    run_build "GRCh38"
else
    run_build "GRCh37"
    run_build "GRCh38"
fi

echo "============================================================"
echo "  Preparation complete!"
echo "============================================================"
echo ""
echo "Output structure:"
echo "  ${OUTDIR}/"
echo "    singularity_images/     # Shared containers"
if [[ "${GENOME_LOWER}" == "grch37" || "${GENOME_LOWER}" == "all" ]]; then
    echo "    GRCh37/"
    echo "      data/                 # GRCh37 BAMs and references"
    echo "      params_GRCh37.yaml    # Ready-to-use params file"
fi
if [[ "${GENOME_LOWER}" == "grch38" || "${GENOME_LOWER}" == "all" ]]; then
    echo "    GRCh38/"
    echo "      data/                 # GRCh38 BAMs and references"
    echo "      params_GRCh38.yaml    # Ready-to-use params file"
fi
echo ""
echo "To run the pipeline:"
if [[ "${GENOME_LOWER}" == "grch37" || "${GENOME_LOWER}" == "all" ]]; then
    echo "  nextflow run main.nf -params-file ${OUTDIR}/GRCh37/params_GRCh37.yaml -profile singularity"
fi
if [[ "${GENOME_LOWER}" == "grch38" || "${GENOME_LOWER}" == "all" ]]; then
    echo "  nextflow run main.nf -params-file ${OUTDIR}/GRCh38/params_GRCh38.yaml -profile singularity"
fi
