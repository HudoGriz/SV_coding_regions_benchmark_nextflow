#!/bin/bash
set -euo pipefail

# =============================================================================
# Download and prepare GRCh37 (hs37d5) data for SV benchmarking
# =============================================================================
#
# Downloads GIAB HG002 BAM files, truth sets, reference genome, and annotations
# for the GRCh37/hs37d5 reference build.
#
# Prerequisites:
#   - Singularity/Apptainer installed
#   - wget installed
#   - ~500 GB free disk space
#   - Stable internet connection (downloads may take several hours)
#
# Usage:
#   Called by prepare.sh, or directly:
#   bash download_and_prep_GRCh37.sh <build_directory> [<singularity_images_directory>] [--skip-download] [--skip-preprocessing]
#
# Arguments:
#   build_directory              Directory for GRCh37-specific data (creates data/ subdirectory)
#   singularity_images_directory Optional: Shared directory for Singularity images
#                                (defaults to <build_directory>/singularity_images)
#
# Output structure:
#   <build_directory>/
#     data/
#       Illumina_wes/bam/    Illumina WES BAM + index
#       Illumina_wgs/bam/    Illumina WGS BAM + index
#       Pacbio/bam/          PacBio HiFi BAM + index
#       ONT/bam/             ONT BAM + index
#       references/          Reference genome, truth sets, BED files
#   <singularity_images_directory>/
#     *.sif                    Singularity container images
# =============================================================================

if [ -z "${1:-}" ]; then
    echo "Usage: bash download_and_prep_GRCh37.sh <build_directory> [<singularity_images_directory>] [--skip-download] [--skip-preprocessing]"
    exit 1
fi

project_dir="$1"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
shift

singularity_dir="${project_dir}/singularity_images"
if [ $# -gt 0 ] && [[ "$1" != --* ]]; then
    singularity_dir="$1"
    shift
fi

SKIP_DOWNLOAD=false
SKIP_PREPROCESSING=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-download)
            SKIP_DOWNLOAD=true
            ;;
        --skip-preprocesing|--skip-preprocessing)
            SKIP_PREPROCESSING=true
            ;;
        *)
            echo "Unknown argument: $1"
            echo "Supported flags: --skip-download, --skip-preprocessing"
            exit 1
            ;;
    esac
    shift
done

echo "=== GRCh37 Data Download and Preparation ==="
echo "Build directory: ${project_dir}"
echo "Singularity images: ${singularity_dir}"
echo "Skip download: ${SKIP_DOWNLOAD}"
echo "Skip preprocessing: ${SKIP_PREPROCESSING}"
echo ""

data_dir="${project_dir}/data"
mkdir -p "${data_dir}"
references_dir="${data_dir}/references"

ensure_singularity_images() {
    echo "--- Ensuring Singularity images ---"
    mkdir -p "${singularity_dir}"
    cd "${singularity_dir}"
    [ -f samtools_latest.sif ] || singularity pull --force samtools_latest.sif docker://quay.io/biocontainers/samtools:1.19--h50ea8bc_0
    [ -f bedtools_latest.sif ] || singularity pull --force bedtools_latest.sif docker://quay.io/biocontainers/bedtools:2.31.1--h13024bc_3
    [ -f r-env_4-4-1.sif ] || singularity pull --force r-env_4-4-1.sif library://blazv/benchmark-sv/r-env:4-4-1
}

download_phase() {
    ensure_singularity_images

    # ---- BAM files ----

    # Illumina WES
    echo "--- Downloading Illumina WES BAM ---"
    mkdir -p "${data_dir}/Illumina_wes/bam"
    cd "${data_dir}/Illumina_wes/bam"
    wget -c https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam
    wget -c https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bai

    # Illumina WGS
    echo "--- Downloading Illumina WGS BAM ---"
    mkdir -p "${data_dir}/Illumina_wgs/bam"
    cd "${data_dir}/Illumina_wgs/bam"
    wget -c https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.60x.1.bam.bai
    wget -c https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.60x.1.bam

    # PacBio HiFi
    echo "--- Downloading PacBio HiFi BAM ---"
    mkdir -p "${data_dir}/Pacbio/bam"
    cd "${data_dir}/Pacbio/bam"
    wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh37.bam
    wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh37.bam.bai

    # ONT
    echo "--- Downloading ONT BAM ---"
    mkdir -p "${data_dir}/ONT/bam"
    cd "${data_dir}/ONT/bam"
    wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/HG002_GRCh37_ONT-UL_UCSC_20200508.phased.bam
    wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/HG002_GRCh37_ONT-UL_UCSC_20200508.phased.bam.bai


    # ---- References ----
    mkdir -p "${references_dir}"
    cd "${references_dir}"

    # GIAB SV truth set v0.6
    echo "--- Downloading GIAB SV truth set ---"
    wget -c https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed
    wget -c https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz

    # Reference genome (hs37d5)
    echo "--- Downloading reference genome (hs37d5) ---"
    wget -c -O human_hs37d5.fasta.gz ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz

    # Tandem repeat annotations
    echo "--- Downloading tandem repeat annotations ---"
    wget -c https://raw.githubusercontent.com/PacificBiosciences/pbsv/master/annotations/human_hs37d5.trf.bed

    # GENCODE GTF for coding regions
    echo "--- Downloading GENCODE v19 GTF ---"
    wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
}

require_file() {
    local file_path="$1"
    if [ ! -f "${file_path}" ]; then
        echo "ERROR: Required file not found: ${file_path}"
        echo "Run with default mode or --skip-preprocessing disabled first."
        exit 1
    fi
}

postprocess_phase() {
    ensure_singularity_images

    require_file "${data_dir}/Illumina_wes/bam/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam"
    require_file "${data_dir}/Illumina_wgs/bam/HG002.hs37d5.60x.1.bam"
    require_file "${data_dir}/Pacbio/bam/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh37.bam"
    require_file "${data_dir}/ONT/bam/HG002_GRCh37_ONT-UL_UCSC_20200508.phased.bam"
    require_file "${references_dir}/HG002_SVs_Tier1_v0.6.bed"
    require_file "${references_dir}/HG002_SVs_Tier1_v0.6.vcf.gz"
    require_file "${references_dir}/gencode.v19.annotation.gtf.gz"

    if [ ! -f "${references_dir}/human_hs37d5.fasta" ]; then
        require_file "${references_dir}/human_hs37d5.fasta.gz"
        gunzip -f "${references_dir}/human_hs37d5.fasta.gz"
    fi

    if [ ! -f "${references_dir}/HG002_SVs_Tier1_v0.6.vcf.gz.tbi" ]; then
        singularity exec \
            "${singularity_dir}/samtools_latest.sif" \
            tabix -p vcf "${references_dir}/HG002_SVs_Tier1_v0.6.vcf.gz"
    fi

    if [ ! -f "${references_dir}/human_hs37d5.fasta.fai" ]; then
        singularity exec \
            "${singularity_dir}/samtools_latest.sif" \
            samtools faidx "${references_dir}/human_hs37d5.fasta"
    fi

# Create exome+UTR BED file (--strip-chr for GRCh37/hs37d5 chromosome naming)
echo "--- Creating exome+UTR BED file ---"
singularity exec \
    -B "${project_dir}" \
    "${singularity_dir}/r-env_4-4-1.sif" \
    Rscript "${SCRIPT_DIR}/create_gencode_target_bed.R" \
        "${references_dir}/gencode.v19.annotation.gtf.gz" \
        "${references_dir}/exome_utr_gtf.bed" \
        --strip-chr

# Intersect exome+UTR with SV truth set
echo "--- Intersecting exome+UTR with truth set ---"
singularity exec \
    -B "${project_dir}" \
    "${singularity_dir}/bedtools_latest.sif" \
    bedtools intersect \
        -a "${references_dir}/HG002_SVs_Tier1_v0.6.bed" \
        -b "${references_dir}/exome_utr_gtf.bed" \
        > "${references_dir}/exome_utr_gtf.HG002_SVs_Tier1.bed"

# Copy BED files from repository data/ if they exist
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_DATA_DIR="${REPO_ROOT}/data"

if [ -f "${REPO_DATA_DIR}/Paediatric_disorders.HG002_SVs_Tier1.GRCh37.bed" ]; then
    echo "--- Copying Paediatric_disorders.HG002_SVs_Tier1.GRCh37.bed from repository ---"
    cp "${REPO_DATA_DIR}/Paediatric_disorders.HG002_SVs_Tier1.GRCh37.bed" "${references_dir}/Paediatric_disorders.HG002_SVs_Tier1.GRCh37.bed"
fi

if [ -f "${REPO_DATA_DIR}/agilent_sureselect_human_all_exon_v5_b37_targets.bed" ]; then
    echo "--- Copying Agilent SureSelect targets from repository ---"
    cp "${REPO_DATA_DIR}/agilent_sureselect_human_all_exon_v5_b37_targets.bed" "${references_dir}/agilent_sureselect_human_all_exon_v5_b37_targets.bed"
fi
}

if [ "${SKIP_DOWNLOAD}" = "true" ] && [ "${SKIP_PREPROCESSING}" = "true" ]; then
    echo "Both download and preprocessing are skipped. Nothing to do."
elif [ "${SKIP_DOWNLOAD}" = "true" ]; then
    postprocess_phase
elif [ "${SKIP_PREPROCESSING}" = "true" ]; then
    download_phase
else
    download_phase
    postprocess_phase
fi


echo ""
echo "=== GRCh37 download and preparation complete ==="
