#!/bin/bash
set -euo pipefail

# =============================================================================
# Download and prepare GRCh38 data for SV benchmarking
# =============================================================================
#
# Downloads GIAB HG002 BAM files, truth sets, reference genome, and annotations
# for the GRCh38 reference build.
#
# Prerequisites:
#   - Singularity/Apptainer installed
#   - wget installed
#   - ~500 GB free disk space
#   - Stable internet connection (downloads may take several hours)
#
# Usage:
#   Called by prepare.sh, or directly:
#   bash download_and_prep_GRCh38.sh <build_directory> [<singularity_images_directory>] [--skip-download] [--skip-preprocessing]
#
# Arguments:
#   build_directory              Directory for GRCh38-specific data (creates data/ subdirectory)
#   singularity_images_directory Optional: Shared directory for Singularity images
#                                (defaults to <build_directory>/singularity_images)
#
# Output structure:
#   <build_directory>/
#     data/
#       Illumina_wgs/bam_GRCh38/   Illumina WGS BAMs + indices
#       Pacbio/bam_GRCh38/          PacBio HiFi BAM + index
#       ONT/bam_GRCh38/             ONT BAM + index
#       filtered_bams/              Filtered BAMs + indices (GRCh38)
#       references/                  Reference genome, truth sets, BED files
#   <singularity_images_directory>/
#     *.sif                    Singularity container images
#
# Note: Illumina WES is not available for GRCh38 in this dataset.
# =============================================================================

if [ -z "${1:-}" ]; then
  echo "Usage: bash download_and_prep_GRCh38.sh <build_directory> [<singularity_images_directory>] [--skip-download] [--skip-preprocessing]"
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

echo "=== GRCh38 Data Download and Preparation ==="
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

  # Illumina WGS (60x only; skip 300x)
  echo "--- Downloading Illumina WGS BAMs (GRCh38) ---"
  mkdir -p "${data_dir}/Illumina_wgs/bam_GRCh38"
  cd "${data_dir}/Illumina_wgs/bam_GRCh38"
  wget -c https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam.bai
  wget -c https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam

  # PacBio HiFi
  echo "--- Downloading PacBio HiFi BAM (GRCh38) ---"
  mkdir -p "${data_dir}/Pacbio/bam_GRCh38"
  cd "${data_dir}/Pacbio/bam_GRCh38"
  wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam.bai
  wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam

  # ONT
  echo "--- Downloading ONT BAM (GRCh38) ---"
  mkdir -p "${data_dir}/ONT/bam_GRCh38"
  cd "${data_dir}/ONT/bam_GRCh38"
  wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/HG002_GRCh38_ONT-UL_UCSC_20200508.phased.bam.bai
  wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/HG002_GRCh38_ONT-UL_UCSC_20200508.phased.bam

  # ---- References ----
  mkdir -p "${references_dir}"
  cd "${references_dir}"

  # GRCh38 truth set (defrabb v0.012)
  echo "--- Downloading GRCh38 truth set ---"
  wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.012-20231107/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed
  wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.012-20231107/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz
  wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.012-20231107/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz.tbi

  # Reference genome (GRCh38 no_alt_analysis_set)
  echo "--- Downloading GRCh38 reference genome ---"
  wget -c -O human_GRCh38_no_alt_analysis_set.fasta.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

  # Tandem repeat annotations (GRCh38)
  echo "--- Downloading GRCh38 tandem repeat annotations ---"
  wget -c https://raw.githubusercontent.com/PacificBiosciences/pbsv/refs/heads/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed

  # GENCODE GTF for coding regions
  echo "--- Downloading GENCODE v49 GTF ---"
  wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v49.annotation.gtf.gz
}

require_file() {
  local file_path="$1"
  if [ ! -f "${file_path}" ]; then
    echo "ERROR: Required file not found: ${file_path}"
    echo "Run with default mode or --download-only first."
    exit 1
  fi
}

postprocess_phase() {
  ensure_singularity_images

  require_file "${data_dir}/Illumina_wgs/bam_GRCh38/HG002.GRCh38.60x.1.bam"
  require_file "${data_dir}/Pacbio/bam_GRCh38/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam"
  require_file "${data_dir}/ONT/bam_GRCh38/HG002_GRCh38_ONT-UL_UCSC_20200508.phased.bam"
  require_file "${references_dir}/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed"
  require_file "${references_dir}/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz"
  require_file "${references_dir}/gencode.v49.annotation.gtf.gz"

  if [ ! -f "${references_dir}/human_GRCh38_no_alt_analysis_set.fasta" ]; then
    require_file "${references_dir}/human_GRCh38_no_alt_analysis_set.fasta.gz"
    gunzip -f "${references_dir}/human_GRCh38_no_alt_analysis_set.fasta.gz"
  fi

  if [ ! -f "${references_dir}/human_GRCh38_no_alt_analysis_set.fasta.fai" ]; then
    singularity exec \
      "${singularity_dir}/samtools_latest.sif" \
      samtools faidx "${references_dir}/human_GRCh38_no_alt_analysis_set.fasta"
  fi

  if [ ! -f "${references_dir}/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz.tbi" ]; then
    singularity exec \
      "${singularity_dir}/samtools_latest.sif" \
      tabix -p vcf "${references_dir}/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz"
  fi

# ---- BAM filtering (GRCh38 only, independent of GRCh37 workflow) ----
echo "--- Filtering GRCh38 BAMs ---"

SAMTOOLS_IMAGE="${singularity_dir}/samtools_latest.sif"
REFERENCE="${references_dir}/human_GRCh38_no_alt_analysis_set.fasta"
filtered_dir="${data_dir}/filtered_bams"
mkdir -p "${filtered_dir}"

CHR_LIST=$(awk '{print $1}' "${REFERENCE}.fai" | tr '\n' ' ')

filter_bam_with_header() {
  local input_bam="$1"
  local output_bam="$2"
  local threads="${3:-4}"

  echo "Filtering paired-end BAM: ${input_bam}"

  (
    singularity exec "${SAMTOOLS_IMAGE}" samtools view -H "${input_bam}" | grep -v "^@SQ"
    for chr in ${CHR_LIST}; do
      singularity exec "${SAMTOOLS_IMAGE}" samtools view -H "${input_bam}" | grep "^@SQ" | grep -w "SN:${chr}"
    done
  ) > "${output_bam}.header.sam"

  singularity exec "${SAMTOOLS_IMAGE}" samtools view -@ "${threads}" -F 3852 -f 2 "${input_bam}" ${CHR_LIST} | \
    cat "${output_bam}.header.sam" - | \
    singularity exec "${SAMTOOLS_IMAGE}" samtools view -b -@ "${threads}" -o "${output_bam}"

  rm -f "${output_bam}.header.sam"
  singularity exec "${SAMTOOLS_IMAGE}" samtools index -@ "${threads}" "${output_bam}"
}

filter_bam_pacbio() {
  local input_bam="$1"
  local output_bam="$2"
  local threads="${3:-4}"

  echo "Filtering long-read BAM: ${input_bam}"

  (
    singularity exec "${SAMTOOLS_IMAGE}" samtools view -H "${input_bam}" | grep -v "^@SQ"
    for chr in ${CHR_LIST}; do
      singularity exec "${SAMTOOLS_IMAGE}" samtools view -H "${input_bam}" | grep "^@SQ" | grep -w "SN:${chr}"
    done
  ) > "${output_bam}.header.sam"

  singularity exec "${SAMTOOLS_IMAGE}" samtools view -@ "${threads}" -F 2308 -q 1 "${input_bam}" ${CHR_LIST} | \
    cat "${output_bam}.header.sam" - | \
    singularity exec "${SAMTOOLS_IMAGE}" samtools view -b -@ "${threads}" -o "${output_bam}"

  rm -f "${output_bam}.header.sam"
  singularity exec "${SAMTOOLS_IMAGE}" samtools index -@ "${threads}" "${output_bam}"
}

illumina_raw_bam="${data_dir}/Illumina_wgs/bam_GRCh38/HG002.GRCh38.60x.1.bam"
pacbio_raw_bam="${data_dir}/Pacbio/bam_GRCh38/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam"
ont_raw_bam="${data_dir}/ONT/bam_GRCh38/HG002_GRCh38_ONT-UL_UCSC_20200508.phased.bam"

illumina_filtered_bam="${filtered_dir}/HG002.Illumina.60.filtered.strict.bam"
pacbio_filtered_bam="${filtered_dir}/HG002.PacBio.filtered.header.strict.pacbiospec.bam"
ont_filtered_bam="${filtered_dir}/HG002.ONT.filtered.header.strict.longread.bam"

filter_bam_with_header \
  "${illumina_raw_bam}" \
  "${illumina_filtered_bam}" \
  30

filter_bam_pacbio \
  "${pacbio_raw_bam}" \
  "${pacbio_filtered_bam}" \
  30

filter_bam_pacbio \
  "${ont_raw_bam}" \
  "${ont_filtered_bam}" \
  30

# Create exome+UTR BED file (no --strip-chr: GRCh38 uses chr1,chr2,... natively)
echo "--- Creating exome+UTR BED file ---"
singularity exec \
  -B "${project_dir}" \
  "${singularity_dir}/r-env_4-4-1.sif" \
  Rscript "${SCRIPT_DIR}/create_gencode_target_bed.R" \
  "${references_dir}/gencode.v49.annotation.gtf.gz" \
  "${references_dir}/exome_utr_gtf_GRCh38.bed"

# Intersect exome+UTR with SV truth set benchmark regions
echo "--- Intersecting exome+UTR with truth set ---"
singularity exec \
  -B "${project_dir}" \
  "${singularity_dir}/bedtools_latest.sif" \
  bedtools intersect \
  -a "${references_dir}/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed" \
  -b "${references_dir}/exome_utr_gtf_GRCh38.bed" \
  >"${references_dir}/exome_utr_gtf.GRCh38_HG002-T2TQ100-V1.0_stvar.bed"

# Copy BED files from repository data/ if they exist
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_DATA_DIR="${REPO_ROOT}/data"

if [ -f "${REPO_DATA_DIR}/Paediatric_disorders.HG002_SVs_Tier1.GRCh38.bed" ]; then
  echo "--- Copying Paediatric_disorders.HG002_SVs_Tier1.GRCh38.bed from repository ---"
  cp "${REPO_DATA_DIR}/Paediatric_disorders.HG002_SVs_Tier1.GRCh38.bed" "${references_dir}/Paediatric_disorders.HG002_SVs_Tier1.GRCh38.bed"
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
echo "=== GRCh38 download and preparation complete ==="
