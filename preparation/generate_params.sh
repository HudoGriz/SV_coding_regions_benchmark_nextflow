#!/bin/bash
set -euo pipefail

# =============================================================================
# Generate a Nextflow params YAML file from downloaded data
# =============================================================================
#
# Scans the data directory structure created by download_and_prep_GRCh37.sh or
# download_and_prep_GRCh38.sh and writes a params YAML with absolute paths.
#
# Usage:
#   Called by prepare.sh, or directly:
#   bash generate_params.sh --genome GRCh37|GRCh38 --datadir <data_directory> --outfile <params.yaml>
#
# The generated YAML is ready to use:
#   nextflow run main.nf -params-file params_GRCh37.yaml -profile singularity
# =============================================================================

usage() {
    echo "Usage: bash generate_params.sh --genome <GRCh37|GRCh38> --datadir <path> --outfile <path>"
    echo ""
    echo "Arguments:"
    echo "  --genome   Reference genome build (GRCh37 or GRCh38)"
    echo "  --datadir  Path to data directory (parent of Illumina_wes/, references/, etc.)"
    echo "  --outfile  Output YAML file path"
    exit 1
}

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
GENOME=""
DATA_DIR=""
OUTFILE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome)  GENOME="$2";   shift 2 ;;
        --datadir) DATA_DIR="$2"; shift 2 ;;
        --outfile) OUTFILE="$2";  shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

if [[ -z "${GENOME}" || -z "${DATA_DIR}" || -z "${OUTFILE}" ]]; then
    echo "ERROR: --genome, --datadir, and --outfile are all required."
    usage
fi

if [[ "${GENOME}" != "GRCh37" && "${GENOME}" != "GRCh38" ]]; then
    echo "ERROR: --genome must be GRCh37 or GRCh38 (got: ${GENOME})"
    exit 1
fi

# Resolve to absolute path
DATA_DIR="$(cd "${DATA_DIR}" && pwd)"
REFS="${DATA_DIR}/references"

# Calculate build directory (parent of DATA_DIR) and parent dir for absolute outdir path
BUILD_DIR="$(dirname "${DATA_DIR}")"
PARENT_DIR="$(dirname "${BUILD_DIR}")"

# Get the repository root directory (for fallback to repo data/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_DATA="${REPO_ROOT}/data"

# ---------------------------------------------------------------------------
# Helper: find a single file matching a glob, or return empty string
# ---------------------------------------------------------------------------
find_file() {
    local pattern="$1"
    # shellcheck disable=SC2086
    local matches=( ${pattern} )
    if [[ -f "${matches[0]:-}" ]]; then
        echo "${matches[0]}"
    else
        echo ""
    fi
}

# ---------------------------------------------------------------------------
# Helper: find file in REFS, fallback to repo data/ directory
# ---------------------------------------------------------------------------
find_file_with_fallback() {
    local refs_path="$1"
    local repo_filename="$2"
    
    # First try the refs directory
    local result=$(find_file "${refs_path}")
    if [[ -n "${result}" ]]; then
        echo "${result}"
        return
    fi
    
    # Fallback to repo data directory
    if [[ -n "${repo_filename}" && -f "${REPO_DATA}/${repo_filename}" ]]; then
        echo "${REPO_DATA}/${repo_filename}"
    else
        echo ""
    fi
}

# ---------------------------------------------------------------------------
# Resolve file paths per genome build
# ---------------------------------------------------------------------------

if [[ "${GENOME}" == "GRCh37" ]]; then
    FASTA=$(find_file "${REFS}/human_hs37d5.fasta")
    ILLUMINA_WES_BAM=$(find_file "${DATA_DIR}/Illumina_wes/bam/*.bam")
    WES_SEQ_TARGETS=$(find_file_with_fallback "${REFS}/agilent_sureselect_human_all_exon_v5_b37_targets.bed" "agilent_sureselect_human_all_exon_v5_b37_targets.bed")
    ILLUMINA_WGS_BAM=$(find_file "${DATA_DIR}/Illumina_wgs/bam/HG002.hs37d5.60x.1.bam")
    PACBIO_BAM=$(find_file "${DATA_DIR}/Pacbio/bam/HG002_PacBio-HiFi-Revio*GRCh37.bam")
    ONT_BAM=$(find_file "${DATA_DIR}/ONT/bam/HG002_GRCh37_ONT*.bam")
    BENCHMARK_VCF=$(find_file "${REFS}/HG002_SVs_Tier1_v0.6.vcf.gz")
    HC_TARGETS=$(find_file "${REFS}/HG002_SVs_Tier1_v0.6.bed")
    GENE_PANEL=$(find_file_with_fallback "${REFS}/Paediatric_disorders.HG002_SVs_Tier1.GRCh37.bed" "Paediatric_disorders.HG002_SVs_Tier1.GRCh37.bed")
    WES_UTR=$(find_file "${REFS}/exome_utr_gtf.HG002_SVs_Tier1.bed")
    TANDEM_REPEATS=$(find_file "${REFS}/human_hs37d5.trf.bed")
    RUN_NAME="GRCh37"
else
    # GRCh38
    FASTA=$(find_file "${REFS}/human_GRCh38_no_alt_analysis_set.fasta")
    ILLUMINA_WES_BAM=""  # No WES available for GRCh38
    WES_SEQ_TARGETS=""
    ILLUMINA_WGS_BAM=$(find_file "${DATA_DIR}/filtered_bams/HG002.Illumina.60.filtered.strict.bam")
    PACBIO_BAM=$(find_file "${DATA_DIR}/filtered_bams/HG002.PacBio.filtered.header.strict.pacbiospec.bam")
    ONT_BAM=$(find_file "${DATA_DIR}/filtered_bams/HG002.ONT.filtered.header.strict.longread.bam")
    if [[ -z "${ILLUMINA_WGS_BAM}" ]]; then
        ILLUMINA_WGS_BAM=$(find_file "${DATA_DIR}/Illumina_wgs/bam_GRCh38/HG002.GRCh38.60x.1.bam")
    fi
    if [[ -z "${PACBIO_BAM}" ]]; then
        PACBIO_BAM=$(find_file "${DATA_DIR}/Pacbio/bam_GRCh38/HG002_PacBio-HiFi-Revio*GRCh38*.bam")
    fi
    if [[ -z "${ONT_BAM}" ]]; then
        ONT_BAM=$(find_file "${DATA_DIR}/ONT/bam_GRCh38/HG002_GRCh38_ONT*.bam")
    fi
    BENCHMARK_VCF=$(find_file "${REFS}/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz")
    HC_TARGETS=$(find_file "${REFS}/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed")
    GENE_PANEL=$(find_file_with_fallback "${REFS}/Paediatric_disorders.HG002_SVs_Tier1.GRCh38.bed" "Paediatric_disorders.HG002_SVs_Tier1.GRCh38.bed")
    WES_UTR=$(find_file "${REFS}/exome_utr_gtf.GRCh38_HG002-T2TQ100-V1.0_stvar.bed")
    TANDEM_REPEATS=$(find_file "${REFS}/human_GRCh38_no_alt_analysis_set.trf.bed")
    RUN_NAME="GRCh38"
fi

# ---------------------------------------------------------------------------
# Validate required files
# ---------------------------------------------------------------------------
MISSING=0
if [[ -z "${FASTA}" ]]; then
    echo "WARNING: Reference FASTA not found for ${GENOME}"
    MISSING=1
fi
if [[ -z "${BENCHMARK_VCF}" ]]; then
    echo "WARNING: Benchmark VCF not found for ${GENOME}"
    MISSING=1
fi
if [[ ${MISSING} -eq 1 ]]; then
    echo "Some required files are missing. The generated YAML will have empty values for these."
    echo "Re-run the download script or provide files manually."
fi

# ---------------------------------------------------------------------------
# Helper: emit a YAML key-value pair (quoted path or null)
# ---------------------------------------------------------------------------
yaml_path() {
    local key="$1"
    local value="$2"
    if [[ -n "${value}" ]]; then
        echo "${key}: '${value}'"
    else
        echo "# ${key}: null  # not found — provide manually if needed"
    fi
}

yaml_value() {
    local key="$1"
    local value="$2"
    echo "${key}: ${value}"
}

# ---------------------------------------------------------------------------
# Write YAML
# ---------------------------------------------------------------------------
echo "Generating ${OUTFILE} for ${GENOME}..."

cat > "${OUTFILE}" << YAML_HEADER
# =============================================================================
# Pipeline parameters for ${GENOME}
# Generated by generate_params.sh on $(date -Iseconds)
# =============================================================================
# Usage:
#   nextflow run main.nf -params-file ${OUTFILE##*/} -profile singularity
# =============================================================================

YAML_HEADER

{
    # Run identification
    yaml_value "run_name" "'${RUN_NAME}'"
    yaml_value "outdir"   "'${PARENT_DIR}/results-${RUN_NAME}'"
    echo ""

    # Reference genome
    echo "# Reference genome (REQUIRED)"
    yaml_path "fasta" "${FASTA}"
    echo ""

    # Input BAM files
    echo "# Input BAM files (provide at least one)"
    yaml_path "illumina_wes_bam"      "${ILLUMINA_WES_BAM}"
    yaml_path "wes_sequencing_targets" "${WES_SEQ_TARGETS}"
    yaml_path "illumina_wgs_bam"      "${ILLUMINA_WGS_BAM}"
    yaml_path "pacbio_bam"            "${PACBIO_BAM}"
    yaml_path "ont_bam"               "${ONT_BAM}"
    echo ""

    # Benchmarking files
    echo "# Benchmarking files"
    yaml_value "skip_benchmarking" "false"
    yaml_path "benchmark_vcf"         "${BENCHMARK_VCF}"
    yaml_path "high_confidence_targets" "${HC_TARGETS}"
    yaml_path "gene_panel_targets"    "${GENE_PANEL}"
    yaml_path "wes_utr_targets"       "${WES_UTR}"
    echo ""

    # Tandem repeats
    echo "# Tandem repeat annotations (recommended for Sniffles)"
    yaml_path "tandem_repeats" "${TANDEM_REPEATS}"
    echo ""

    # Truvari parameters
    echo "# Truvari benchmarking parameters"
    yaml_value "truvari_refdist"  "500"
    yaml_value "truvari_pctsize"  "0.7"
    yaml_value "truvari_pctovl"   "0"
    yaml_value "truvari_pctseq"   "0"
    echo ""

    echo "# WES-specific Truvari parameters"
    yaml_value "truvari_wes_refdist"  "500"
    yaml_value "truvari_wes_pctsize"  "0"
    yaml_value "truvari_wes_pctovl"   "0"
    yaml_value "truvari_wes_pctseq"   "0"
    echo ""

    # Resource limits
    echo "# Resource limits (adjust to your cluster)"
    yaml_value "max_cpus"   "48"
    yaml_value "max_memory"  "'128.GB'"
    yaml_value "max_time"    "'48.h'"
    echo ""

    # Simulation
    echo "# Simulation (optional)"
    yaml_value "simulate_targets" "true"
    yaml_value "num_simulations"  "500"
    echo ""

    # Analysis
    echo "# Analysis"
    yaml_value "gather_statistics" "true"
} >> "${OUTFILE}"

echo "Done. Review ${OUTFILE} and adjust resource limits for your system."
