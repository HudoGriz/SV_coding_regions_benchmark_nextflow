#!/usr/bin/env Rscript


# Load necessary libraries
suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
})

chr <- c(1:24, "X", "Y")


gtf <- import("gencode.v19.annotation.gtf.gz", format = "gtf")

# Extract exonic regions
exons <- gtf[gtf$type == "exon"]
exons <- GRanges(
    seqnames = gsub("chr", "", seqnames(exons)),
    ranges = ranges(exons)
    )

exons <- exons[seqnames(exons) %in% chr, ]

exons <- reduce(exons)

# Extract UTR regions
utr <- gtf[gtf$type == "UTR"]
utr <- GRanges(
    seqnames = gsub("chr", "", seqnames(utr)),
    ranges = ranges(utr)
    )

# Combine exons and UTRs into one GenomeRanges object
exome_and_utrs <- c(exons, utr)

exome_and_utrs <- reduce(exome_and_utrs)
exome_and_utrs <- exome_and_utrs[seqnames(exome_and_utrs) != "M", ]

# Save to a file
export(exome_and_utrs, "exome_utr_gtf.bed", format = "bed")
