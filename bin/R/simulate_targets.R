library(parallel)
library(rtracklayer)
library(GenomicRanges)

# Import function from general_functions.R
import_bed <- function(file_path, mcols = FALSE) {
    bed <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    
    gr <- GRanges(
        seqnames = bed[[1]],
        ranges = IRanges(start = bed[[2]], end = bed[[3]]),
        strand = "*"
    )

    if (mcols) {
        mcols(gr) <- bed[, -c(1:4)]
    }

    return(gr)
}

# Set the seed for reproducibility
set.seed(42)


# Fetch command-line arguments
args <- commandArgs(trailingOnly = TRUE)

reps <- as.numeric(args[1])
out_file_path <- args[2]
wes_utr_targets <- args[3]
high_confidence_targets <- args[4]


simulate_targets <- function(repetition_id) {
    # Import exome and high confidence intervals
    exome_gr <- import_bed(wes_utr_targets)
    high_confidence_gr <- import_bed(high_confidence_targets)

    n_splits <- 10


    n_exons <- length(exome_gr)
    partial <- n_exons / n_splits

    exons_widths <- data.frame(
        widths = sort(width(exome_gr)),
        split = (1:n_exons %/% partial) + 1
    )

    exons_freq <- data.frame(table(seqnames(exome_gr)))
    exons_freq$partial <- exons_freq$Freq / sum(exons_freq$Freq)

    # Remove exons from high confidence intervals
    gr_filtered <- setdiff(high_confidence_gr, exome_gr, ignore.strand = TRUE)
    gr_filtered <- high_confidence_gr
    all_chunks <- GRanges()

    # Simulate targets
    for (i in n_splits:1) {
        exon_length <- median(exons_widths$widths[exons_widths$split == i])

        for (chr in as.character(unique(seqnames(gr_filtered)))) {

            random_chunks <- slidingWindows(
                keepSeqlevels(gr_filtered, chr, pruning.mode = "coarse"),
                width = exon_length,
                step = exon_length
            )
            random_chunks <- unlist(random_chunks)

            ok <- width(random_chunks) >= (exon_length / 1.25)
            random_chunks <- random_chunks[ok]

            selected <- sample(
                1:length(random_chunks),
                round(partial * exons_freq[exons_freq$Var1 == chr, "partial"]),
                replace = FALSE
                )
            random_chunks <- random_chunks[selected]

            all_chunks <- c(all_chunks, random_chunks)
        }
        # print(all_chunks)

        gr_filtered <- setdiff(gr_filtered, all_chunks, ignore.strand = TRUE)

        print(paste0("Rep ", repetition_id, " - Split ", i, " done"))
    }

    all_chunks <- sort(all_chunks)

    # Save the simulated targets for the current repetition
    output_file <- file.path(paste0(out_file_path, "/simulation", repetition_id, ".bed"))
    
    # Prepare output data frame
    output_df <- data.frame(
        seqnames = as.character(seqnames(all_chunks)),
        start = start(all_chunks),
        end = end(all_chunks),
        stringsAsFactors = FALSE
    )
    
    write.table(
        output_df,
        output_file,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
    )
}

# Use mclapply to run the repetitions in parallel
mclapply(1:reps, simulate_targets, mc.cores = detectCores() - 1)
