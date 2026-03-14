#!/usr/bin/env Rscript

# Create exome + UTR BED file from GENCODE GTF annotation
#
# Usage:
#   Rscript create_gencode_target_bed.R <gtf_path> <output_bed_path> [--strip-chr]
#
# Arguments:
#   gtf_path        Path to GENCODE GTF file (can be .gtf.gz)
#   output_bed_path Path for output BED file
#   --strip-chr     Optional: strip 'chr' prefix from chromosome names (for GRCh37/hs37d5)
#
# Examples:
#   # GRCh37 (hs37d5 uses 1,2,3... without chr prefix):
#   Rscript create_gencode_target_bed.R gencode.v19.annotation.gtf.gz exome_utr_gtf.bed --strip-chr
#
#   # GRCh38 (uses chr1,chr2,chr3...):
#   Rscript create_gencode_target_bed.R gencode.v49.annotation.gtf.gz exome_utr_gtf_GRCh38.bed

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Usage: Rscript create_gencode_target_bed.R <gtf_path> <output_bed_path> [--strip-chr]")
}

gtf_path <- args[1]
output_bed_path <- args[2]
strip_chr <- "--strip-chr" %in% args

if (!file.exists(gtf_path)) {
    stop(paste("GTF file not found:", gtf_path))
}

cat("Reading GTF:", gtf_path, "\n")
gtf <- import(gtf_path, format = "gtf")

# Define valid chromosomes (with chr prefix, will strip later if needed)
valid_chr <- c(paste0("chr", c(1:22, "X", "Y")))

# Extract exonic regions
exons <- gtf[gtf$type == "exon"]
exons <- exons[as.character(seqnames(exons)) %in% valid_chr, ]
exons <- reduce(exons)

# Extract UTR regions
utr <- gtf[gtf$type == "UTR"]
utr <- utr[as.character(seqnames(utr)) %in% valid_chr, ]

# Combine exons and UTRs
exome_and_utrs <- c(exons, utr)
exome_and_utrs <- reduce(exome_and_utrs)

# Strip chr prefix if requested (for GRCh37/hs37d5 compatibility)
if (strip_chr) {
    cat("Stripping 'chr' prefix from chromosome names\n")
    seqlevels(exome_and_utrs) <- gsub("^chr", "", seqlevels(exome_and_utrs))
}

cat("Writing", length(exome_and_utrs), "regions to:", output_bed_path, "\n")
export(exome_and_utrs, output_bed_path, format = "bed")

cat("Done.\n")
