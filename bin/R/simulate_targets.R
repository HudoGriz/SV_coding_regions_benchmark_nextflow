#!/usr/bin/env Rscript
#
# Generate interval-matched noncoding target sets ("EX+UTR-like" BED files).
#
# History: the first published run of this script requested 500 repetitions but
# delivered only 111 (GRCh37) and 114 (GRCh38). The cause was `mc.cores =
# detectCores() - 1`: inside the Singularity container detectCores() reports the
# whole compute node (256 threads), not the 24 CPUs Nextflow requested, so 255 R
# workers were forked into a cgroup with a 64 GB memory request. slurmstepd
# logged 1414 and 3079 oom-kill events respectively, mclapply reported that
# ~200 of 255 workers "did not deliver results", and the task still exited 0 --
# so Nextflow published a partial set as if it were complete.
#
# The three changes that prevent a recurrence:
#   1. the worker count is passed in from Nextflow (task.cpus); there is no
#      detectCores() fallback;
#   2. the invariant setup (BED import, exon deciles, per-chromosome counts and
#      the HCI-minus-exome sampling pool) is computed once before forking, so the
#      workers share it copy-on-write instead of rebuilding it 500 times;
#   3. the script fails loudly. If any repetition errors, or if the number of
#      BED files written does not equal the number requested, it stops with a
#      non-zero exit status so the pipeline fails instead of silently publishing
#      an incomplete null.

library(parallel)
library(rtracklayer)
library(GenomicRanges)

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
    stop("usage: simulate_targets.R <reps> <outdir> <wes_utr.bed> <high_confidence.bed> [n_workers]")
}

reps                    <- as.integer(args[1])
out_file_path           <- args[2]
wes_utr_targets         <- args[3]
high_confidence_targets <- args[4]

# Never fall back to detectCores(): see the header note.
n_workers <- if (length(args) >= 5 && nzchar(args[5])) as.integer(args[5]) else 1L
if (is.na(n_workers) || n_workers < 1L) n_workers <- 1L

cat(sprintf("simulate_targets: %d repetitions, %d worker(s), output -> %s\n",
            reps, n_workers, out_file_path))

# L'Ecuyer-CMRG is the only generator for which mclapply's per-worker streams
# are reproducible. With R's default generator, set.seed() does not make forked
# results repeatable, so the previous run was not reproducible either.
RNGkind("L'Ecuyer-CMRG")
set.seed(42)

dir.create(out_file_path, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Invariant setup, computed once and shared with the workers copy-on-write
# ---------------------------------------------------------------------------
n_splits <- 10

exome_gr           <- import_bed(wes_utr_targets)
high_confidence_gr <- import_bed(high_confidence_targets)

n_exons <- length(exome_gr)
partial <- n_exons / n_splits

exons_widths <- data.frame(
    widths = sort(width(exome_gr)),
    split = (seq_len(n_exons) %/% partial) + 1
)

exons_freq <- data.frame(table(seqnames(exome_gr)))
exons_freq$partial <- exons_freq$Freq / sum(exons_freq$Freq)

# Remove exons from high confidence intervals: this is the sampling pool every
# repetition starts from.
gr_filtered_base <- setdiff(high_confidence_gr, exome_gr, ignore.strand = TRUE)

decile_lengths <- vapply(
    seq_len(n_splits),
    function(i) median(exons_widths$widths[exons_widths$split == i]),
    numeric(1)
)

cat(sprintf("simulate_targets: %d EX+UTR intervals, sampling pool of %d noncoding HCI intervals\n",
            n_exons, length(gr_filtered_base)))

# ---------------------------------------------------------------------------
# One repetition
# ---------------------------------------------------------------------------
simulate_targets <- function(repetition_id) {
    gr_filtered <- gr_filtered_base
    all_chunks <- GRanges()

    for (i in n_splits:1) {
        exon_length <- decile_lengths[i]

        for (chr in as.character(unique(seqnames(gr_filtered)))) {

            random_chunks <- slidingWindows(
                keepSeqlevels(gr_filtered, chr, pruning.mode = "coarse"),
                width = exon_length,
                step = exon_length
            )
            random_chunks <- unlist(random_chunks)

            ok <- width(random_chunks) >= (exon_length / 1.25)
            random_chunks <- random_chunks[ok]

            # Guards for the degenerate cases. The original code called
            # sample(1:length(x), n) directly, which samples from c(1, 0) when
            # the pool is empty and errors outright when n exceeds the pool.
            if (length(random_chunks) == 0L) next

            n_wanted <- exons_freq$partial[exons_freq$Var1 == chr]
            if (length(n_wanted) == 0L) next
            n_wanted <- round(partial * n_wanted)
            if (is.na(n_wanted) || n_wanted < 1L) next

            if (n_wanted > length(random_chunks)) {
                warning(sprintf(
                    "Rep %d, decile %d, %s: %d chunks wanted but only %d available; taking all",
                    repetition_id, i, chr, n_wanted, length(random_chunks)
                ))
                n_wanted <- length(random_chunks)
            }

            selected <- sample.int(length(random_chunks), n_wanted, replace = FALSE)
            all_chunks <- c(all_chunks, random_chunks[selected])
        }

        # Intervals already drawn leave the pool, so none is reused within a set.
        gr_filtered <- setdiff(gr_filtered, all_chunks, ignore.strand = TRUE)
    }

    all_chunks <- sort(all_chunks)

    output_df <- data.frame(
        seqnames = as.character(seqnames(all_chunks)),
        start = start(all_chunks),
        end = end(all_chunks),
        stringsAsFactors = FALSE
    )

    # Write to a temporary name and rename, so a worker killed mid-write can
    # never leave a truncated BED file that looks complete.
    final_file <- file.path(out_file_path, paste0("simulation", repetition_id, ".bed"))
    tmp_file <- paste0(final_file, ".partial")

    write.table(output_df, tmp_file, sep = "\t", row.names = FALSE,
                col.names = FALSE, quote = FALSE)
    if (!file.rename(tmp_file, final_file)) {
        stop(sprintf("Rep %d: could not rename %s to %s", repetition_id, tmp_file, final_file))
    }

    length(all_chunks)
}

# ---------------------------------------------------------------------------
# Run
#
# mc.preschedule = FALSE forks one worker per repetition, at most n_workers at a
# time. With prescheduling a killed worker takes its whole contiguous block of
# repetitions with it, which is what produced the scattered gaps in the first
# run.
# ---------------------------------------------------------------------------
results <- mclapply(
    seq_len(reps),
    function(i) tryCatch(simulate_targets(i), error = function(e) e),
    mc.cores = n_workers,
    mc.preschedule = FALSE
)

# ---------------------------------------------------------------------------
# Fail loudly
# ---------------------------------------------------------------------------
is_bad <- vapply(results, function(x) inherits(x, "error") || is.null(x), logical(1))

if (any(is_bad)) {
    bad <- which(is_bad)
    msgs <- vapply(results[bad], function(x) {
        if (inherits(x, "error")) conditionMessage(x) else "worker returned NULL (killed?)"
    }, character(1))
    for (k in seq_along(bad)) {
        cat(sprintf("simulate_targets: repetition %d FAILED: %s\n", bad[k], msgs[k]))
    }
}

written <- list.files(out_file_path, pattern = "^simulation[0-9]+\\.bed$")
stray <- list.files(out_file_path, pattern = "\\.partial$")
if (length(stray) > 0) {
    file.remove(file.path(out_file_path, stray))
}

cat(sprintf("simulate_targets: %d of %d repetitions wrote a BED file\n",
            length(written), reps))

if (length(written) != reps || any(is_bad)) {
    stop(sprintf(
        paste0("simulate_targets: expected %d simulated interval sets, got %d ",
               "(%d worker failures). Refusing to hand a partial empirical null to the ",
               "benchmarking step. If workers were killed by the cgroup out-of-memory ",
               "handler, lower the worker count or raise the memory request for ",
               "SIMULATE_TARGETS; check for oom-kill lines in .command.log."),
        reps, length(written), sum(is_bad)
    ))
}

cat("simulate_targets: complete\n")
