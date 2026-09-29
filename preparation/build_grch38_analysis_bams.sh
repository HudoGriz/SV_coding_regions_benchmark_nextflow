#!/bin/bash
set -euo pipefail

# =============================================================================
# Restrict the GRCh38 BAMs to the contigs of the analysis reference
# =============================================================================
#
# The GIAB GRCh38 BAMs were aligned against references with more contigs than
# the no-alt analysis set the pipeline calls against: the Illumina BAM carries
# 2,385 extra decoy and HLA contigs, the PacBio BAM 23 GIABv3 decoys. Manta and
# pbsv refuse a BAM whose header names contigs that are missing from the
# reference, so those contigs have to leave the header. Nothing else does.
#
# The output keeps every record on the analysis contigs whatever its flags:
# discordant pairs, supplementary alignments, reads with an unmapped mate,
# duplicates, secondary and low-MAPQ alignments all stay, because the callers
# use them as evidence and apply their own filters. The only other records
# removed are paired reads whose mate lies on a dropped contig, since their
# mate reference would be undefined once that contig leaves the header.
#
# This replaces the earlier filter (Illumina `-F 3852 -f 2`, long reads
# `-F 2308 -q 1`), which also discarded discordant pairs, mate-unmapped reads
# and supplementary alignments. Its outputs in data/filtered_bams/ are left in
# place. This script writes to data/analysis_bams/ and refuses to overwrite.
#
# A BAM whose header already matches the reference exactly (the ONT BAM) is not
# copied: it is used as distributed, and its manifest says so.
#
# The analysis contigs must form a prefix of the BAM header, in reference
# order. That keeps every record's reference and mate IDs valid when the
# trailing contigs are dropped, so the header can be replaced without decoding
# the reads. The script checks this and stops if it does not hold.
#
# A second stage removes the remaining references to the dropped contigs: SA-tag
# entries (the other parts of a split read) that point at a dropped contig. The
# SAM specification requires SA reference names to be header contigs, and pbsv
# call aborts on them. strip_absent_sa_entries.py removes those entries and
# nothing else; the records themselves are all kept. The restricted BAM of the
# first stage stays in place as the second stage's input.
#
# Each stage is skipped when its manifest already exists, and refuses to
# overwrite a partial output.
#
# Usage:
#   bash build_grch38_analysis_bams.sh <build_directory> [<singularity_images_directory>]
#        [--only illumina|pacbio|ont] [--threads N]
#
# Environment:
#   ANALYSIS_SIF   the python-r-analysis image, if not in the images directory
#
# Outputs, per technology:
#   <build_directory>/data/analysis_bams/HG002.<tech>.GRCh38.analysis_contigs.bam(.bai)
#   <build_directory>/data/analysis_bams/<tech>.manifest.tsv
#   <build_directory>/data/analysis_bams/HG002.<tech>.GRCh38.analysis_contigs.sa_filtered.bam(.bai)
#   <build_directory>/data/analysis_bams/<tech>.sa_filtered.manifest.tsv, .counts.tsv, .dropped_contigs.tsv
# The .sa_filtered BAMs are the pipeline inputs (generate_params.sh).
# =============================================================================

usage() {
  echo "Usage: bash build_grch38_analysis_bams.sh <build_directory> [<singularity_images_directory>] [--only illumina|pacbio|ont] [--threads N]"
  exit 1
}

[ -n "${1:-}" ] || usage
project_dir="$1"
shift

singularity_dir="${project_dir}/singularity_images"
if [ $# -gt 0 ] && [[ "$1" != --* ]]; then
  singularity_dir="$1"
  shift
fi

only=""
threads=8
while [[ $# -gt 0 ]]; do
  case "$1" in
    --only)    only="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    *) echo "Unknown argument: $1"; usage ;;
  esac
done

data_dir="${project_dir}/data"
references_dir="${data_dir}/references"
out_dir="${data_dir}/analysis_bams"
reference_fai="${references_dir}/human_GRCh38_no_alt_analysis_set.fasta.fai"
samtools_image="${singularity_dir}/samtools_latest.sif"
analysis_image="${ANALYSIS_SIF:-${singularity_dir}/python-r-analysis_py3.11-r4.4.1.sif}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
sa_script="${script_dir}/strip_absent_sa_entries.py"

declare -A input_bams=(
  [illumina]="${data_dir}/Illumina_wgs/bam_GRCh38/HG002.GRCh38.60x.1.bam"
  [pacbio]="${data_dir}/Pacbio/bam_GRCh38/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam"
  [ont]="${data_dir}/ONT/bam_GRCh38/HG002_GRCh38_ONT-UL_UCSC_20200508.phased.bam"
)
declare -A labels=([illumina]=Illumina [pacbio]=PacBio [ont]=ONT)

for required in "${reference_fai}" "${samtools_image}" "${analysis_image}" "${sa_script}"; do
  [ -f "${required}" ] || { echo "ERROR: required file not found: ${required}" >&2; exit 1; }
done

samtools() {
  singularity exec "${samtools_image}" samtools "$@"
}

# Print "name<TAB>length" for every @SQ line of a BAM, in header order.
header_contigs() {
  samtools view -H "$1" | awk -F'\t' '$1 == "@SQ" {
    name = ""; len = ""
    for (i = 2; i <= NF; i++) {
      if ($i ~ /^SN:/) name = substr($i, 4)
      if ($i ~ /^LN:/) len = substr($i, 4)
    }
    print name "\t" len
  }'
}

# Sum the records idxstats reports (mapped + unmapped-placed) over the first n contigs.
records_on_first_contigs() {
  samtools idxstats "$1" | awk -v n="$2" 'NR <= n { total += $3 + $4 } END { print total + 0 }'
}

# Value of a key in a two-column manifest.
manifest_value() {
  awk -F'\t' -v key="$2" '$1 == key { print $2 }' "$1"
}

restrict_contigs() {
  local tech="$1"
  local in_bam="${input_bams[$tech]}"
  local label="${labels[$tech]}"
  local out_bam="${out_dir}/HG002.${label}.GRCh38.analysis_contigs.bam"
  local manifest="${out_dir}/${label}.manifest.tsv"
  local work="${out_dir}/.${label}.work"

  if [ -f "${manifest}" ] && [ ! -e "${work}" ]; then
    echo "${label}: contig restriction already done (${manifest})"
    return 0
  fi
  [ -f "${in_bam}" ] || { echo "ERROR: input BAM not found: ${in_bam}" >&2; return 1; }
  for existing in "${out_bam}" "${out_bam}.bai" "${manifest}" "${work}"; do
    if [ -e "${existing}" ]; then
      echo "ERROR: refusing to overwrite ${existing}" >&2
      return 1
    fi
  done
  mkdir -p "${work}"

  local reference_contigs="${work}/reference.contigs.tsv"
  local bam_contigs="${work}/bam.contigs.tsv"
  cut -f1,2 "${reference_fai}" > "${reference_contigs}"
  header_contigs "${in_bam}" > "${bam_contigs}"
  local n_ref n_bam
  n_ref=$(wc -l < "${reference_contigs}")
  n_bam=$(wc -l < "${bam_contigs}")

  {
    printf 'key\tvalue\n'
    printf 'technology\t%s\n' "${label}"
    printf 'input_bam\t%s\n' "${in_bam}"
    printf 'input_bytes\t%s\n' "$(stat -c %s "${in_bam}")"
    printf 'reference_fai\t%s\n' "${reference_fai}"
    printf 'reference_contigs\t%s\n' "${n_ref}"
    printf 'input_header_contigs\t%s\n' "${n_bam}"
    printf 'samtools\t%s\n' "$(samtools --version | head -1)"
    printf 'samtools_image\t%s\n' "${samtools_image}"
    printf 'started_at\t%s\n' "$(date --iso-8601=seconds)"
  } > "${work}/manifest.partial.tsv"

  if cmp -s "${reference_contigs}" "${bam_contigs}"; then
    echo "${label}: header matches the reference exactly; the BAM is used as distributed"
    {
      cat "${work}/manifest.partial.tsv"
      printf 'action\tused_as_distributed\n'
      printf 'output_bam\t%s\n' "${in_bam}"
      printf 'records_removed\t0\n'
      printf 'finished_at\t%s\n' "$(date --iso-8601=seconds)"
    } > "${manifest}"
    rm -rf "${work}"
    return 0
  fi

  if ! head -n "${n_ref}" "${bam_contigs}" | cmp -s - "${reference_contigs}"; then
    echo "ERROR: ${label}: the reference contigs are not a prefix of the BAM header (same names, lengths and order)." >&2
    echo "       Records cannot be kept without renumbering; this script does not do that." >&2
    return 1
  fi

  local contig_list
  contig_list=$(cut -f1 "${reference_contigs}" | tr '\n' ' ')
  local tmp_bam="${work}/restricted.bam"
  local filter_expression="mrefid < ${n_ref}"

  echo "${label}: keeping records on ${n_ref} of ${n_bam} contigs, dropping mates on removed contigs"
  # shellcheck disable=SC2086
  samtools view -@ "${threads}" -b -e "${filter_expression}" -o "${tmp_bam}" "${in_bam}" ${contig_list}

  # Header: @HD, the kept @SQ lines, then every other line (read groups and the
  # full program history, including the samtools view step just run).
  samtools view --no-PG -H "${tmp_bam}" | awk -F'\t' -v n="${n_ref}" '
    $1 == "@HD" { print; next }
    $1 == "@SQ" { if (++sq <= n) print; next }
    { rest[++r] = $0 }
    END { for (i = 1; i <= r; i++) print rest[i] }
  ' > "${work}/header.sam"

  samtools reheader "${work}/header.sam" "${tmp_bam}" > "${work}/reheadered.bam"
  samtools quickcheck "${work}/reheadered.bam"
  samtools index -@ "${threads}" "${work}/reheadered.bam"

  local records_in records_out
  records_in=$(records_on_first_contigs "${in_bam}" "${n_ref}")
  records_out=$(records_on_first_contigs "${work}/reheadered.bam" "${n_ref}")

  header_contigs "${work}/reheadered.bam" > "${work}/output.contigs.tsv"
  if ! cmp -s "${work}/output.contigs.tsv" "${reference_contigs}"; then
    echo "ERROR: ${label}: rebuilt header does not match the reference contigs" >&2
    return 1
  fi

  mv "${work}/reheadered.bam" "${out_bam}"
  mv "${work}/reheadered.bam.bai" "${out_bam}.bai"

  {
    cat "${work}/manifest.partial.tsv"
    printf 'action\trestricted_to_reference_contigs\n'
    printf 'filter\tsamtools view -b -e '\''%s'\'' <input> <%s reference contigs>; samtools reheader\n' "${filter_expression}" "${n_ref}"
    printf 'flag_filters\tnone\n'
    printf 'output_bam\t%s\n' "${out_bam}"
    printf 'output_bytes\t%s\n' "$(stat -c %s "${out_bam}")"
    printf 'records_on_reference_contigs_in\t%s\n' "${records_in}"
    printf 'records_out\t%s\n' "${records_out}"
    printf 'records_removed\t%s\n' "$((records_in - records_out))"
    printf 'finished_at\t%s\n' "$(date --iso-8601=seconds)"
  } > "${manifest}"
  rm -rf "${work}"
  echo "${label}: wrote ${out_bam} (${records_out} of ${records_in} records kept)"
}

filter_sa_entries() {
  local tech="$1"
  local label="${labels[$tech]}"
  local restricted_manifest="${out_dir}/${label}.manifest.tsv"
  local out_bam="${out_dir}/HG002.${label}.GRCh38.analysis_contigs.sa_filtered.bam"
  local manifest="${out_dir}/${label}.sa_filtered.manifest.tsv"
  local counts="${out_dir}/${label}.sa_filtered.counts.tsv"
  local dropped="${out_dir}/${label}.sa_filtered.dropped_contigs.tsv"
  local work="${out_dir}/.${label}.sa_filtered.work"

  if [ "$(manifest_value "${restricted_manifest}" action)" = "used_as_distributed" ]; then
    # The header is the aligner's reference, so every SA entry names a header contig.
    echo "${label}: BAM used as distributed; no contig was dropped, so no SA entry can name one"
    return 0
  fi
  if [ -f "${manifest}" ] && [ ! -e "${work}" ]; then
    echo "${label}: SA filtering already done (${manifest})"
    return 0
  fi
  for existing in "${out_bam}" "${out_bam}.bai" "${manifest}" "${counts}" "${dropped}" "${work}"; do
    if [ -e "${existing}" ]; then
      echo "ERROR: refusing to overwrite ${existing}" >&2
      return 1
    fi
  done

  local in_bam
  in_bam=$(manifest_value "${restricted_manifest}" output_bam)
  [ -f "${in_bam}" ] || { echo "ERROR: restricted BAM not found: ${in_bam}" >&2; return 1; }
  mkdir -p "${work}"

  {
    printf 'key\tvalue\n'
    printf 'technology\t%s\n' "${label}"
    printf 'input_bam\t%s\n' "${in_bam}"
    printf 'input_bytes\t%s\n' "$(stat -c %s "${in_bam}")"
    printf 'script\t%s\n' "${sa_script}"
    printf 'script_sha256\t%s\n' "$(sha256sum "${sa_script}" | cut -d' ' -f1)"
    printf 'analysis_image\t%s\n' "${analysis_image}"
    printf 'analysis_image_sha256\t%s\n' "$(sha256sum "${analysis_image}" | cut -d' ' -f1)"
    printf 'pysam\t%s\n' "$(singularity exec "${analysis_image}" python3 -c 'import pysam; print(pysam.__version__)')"
    printf 'samtools\t%s\n' "$(samtools --version | head -1)"
    printf 'started_at\t%s\n' "$(date --iso-8601=seconds)"
  } > "${work}/manifest.partial.tsv"

  echo "${label}: removing SA entries that name contigs absent from the header"
  singularity exec -B "${out_dir}" -B "${script_dir}" "${analysis_image}" \
    python3 "${sa_script}" --input "${in_bam}" --chunk-dir "${work}/chunks" --workers "${threads}" \
      --counts "${work}/counts.tsv" --dropped "${work}/dropped_contigs.tsv"

  # Chunk names are zero-padded header indices, so the glob is in header order.
  samtools cat -o "${work}/filtered.bam" "${work}"/chunks/*.bam
  samtools quickcheck "${work}/filtered.bam"
  samtools index -@ "${threads}" "${work}/filtered.bam"

  samtools idxstats "${in_bam}" > "${work}/idxstats.in.tsv"
  samtools idxstats "${work}/filtered.bam" > "${work}/idxstats.out.tsv"
  if ! cmp -s "${work}/idxstats.in.tsv" "${work}/idxstats.out.tsv"; then
    echo "ERROR: ${label}: per-contig record counts differ between input and output" >&2
    return 1
  fi
  header_contigs "${work}/filtered.bam" > "${work}/output.contigs.tsv"
  cut -f1,2 "${reference_fai}" > "${work}/reference.contigs.tsv"
  if ! cmp -s "${work}/output.contigs.tsv" "${work}/reference.contigs.tsv"; then
    echo "ERROR: ${label}: output header does not match the reference contigs" >&2
    return 1
  fi

  mv "${work}/filtered.bam" "${out_bam}"
  mv "${work}/filtered.bam.bai" "${out_bam}.bai"
  mv "${work}/counts.tsv" "${counts}"
  mv "${work}/dropped_contigs.tsv" "${dropped}"

  local total_line
  total_line=$(awk -F'\t' '$1 == "total"' "${counts}")
  {
    cat "${work}/manifest.partial.tsv"
    printf 'action\tsa_entries_naming_absent_contigs_removed\n'
    printf 'output_bam\t%s\n' "${out_bam}"
    printf 'output_bytes\t%s\n' "$(stat -c %s "${out_bam}")"
    printf 'records\t%s\n' "$(cut -f2 <<< "${total_line}")"
    printf 'records_with_sa\t%s\n' "$(cut -f3 <<< "${total_line}")"
    printf 'sa_entries_removed\t%s\n' "$(cut -f4 <<< "${total_line}")"
    printf 'sa_tags_shortened\t%s\n' "$(cut -f5 <<< "${total_line}")"
    printf 'sa_tags_deleted\t%s\n' "$(cut -f6 <<< "${total_line}")"
    printf 'records_removed\t0\n'
    printf 'counts\t%s\n' "${counts}"
    printf 'dropped_contigs\t%s\n' "${dropped}"
    printf 'finished_at\t%s\n' "$(date --iso-8601=seconds)"
  } > "${manifest}"
  rm -rf "${work}"
  echo "${label}: wrote ${out_bam}"
}

build_one() {
  restrict_contigs "$1"
  filter_sa_entries "$1"
}

mkdir -p "${out_dir}"
if [ -n "${only}" ]; then
  [ -n "${input_bams[$only]:-}" ] || { echo "ERROR: --only must be illumina, pacbio or ont" >&2; exit 1; }
  build_one "${only}"
else
  for tech in illumina pacbio ont; do
    build_one "${tech}"
  done
fi
