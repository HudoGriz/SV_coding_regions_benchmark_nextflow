# Load necessary libraries
suppressPackageStartupMessages({
    library(GenomicRanges)
    library(caret)
    library(vcfR)
    library(jsonlite)
})


calculate_overlap <- function(gr1, gr2) {
    # Merge overlapping intervals within each GRanges object
    gr1 <- reduce(gr1)
    gr2 <- reduce(gr2)

    # Find overlaps
    overlaps <- findOverlaps(gr1, gr2)

    # Calculate the total overlap length
    overlap_ranges <- pintersect(gr1[queryHits(overlaps)], gr2[subjectHits(overlaps)])
    total_overlap_length <- sum(width(overlap_ranges))

    print(paste0("gr1 overlap: ", total_overlap_length / sum(width(gr1))))
    print(paste0("gr2 overlap: ", total_overlap_length / sum(width(gr2))))

    # Print the total overlap length
    return(total_overlap_length)
}


import_bed <- function(file_path, mcols = FALSE) {
    bed <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    # chr <- ifelse(any(grepl("chr", bed[[1]])), "", "chr")

    gr <- GRanges(
        # seqnames = paste0(chr, bed[[1]]),
        seqnames = gsub("chr", "", bed[[1]]),
        ranges = IRanges(start = bed[[2]], end = bed[[3]]),
        strand = "*"
    )

    if (mcols) {
        mcols(gr) <- bed[, -c(1:4)]
    }

    return(gr)
}


show_overlap_distribution <- function(percentOverlap) {
    print(paste0("1 - 0.75: ", sum(percentOverlap <= 1 & percentOverlap > 0.75)))
    print(paste0("0.75 - 0.5: ", sum(percentOverlap < 0.75 & percentOverlap > 0.5)))
    print(paste0("0.5 - 0.25: ", sum(percentOverlap < 0.5 & percentOverlap > 0.25)))
    print(paste0("0.25 - 0: ", sum(percentOverlap < 0.25)))
}


# Import SV reference vcf
import_vcf <- function(vcf_file, add_sample_info_columns = c(), filter = "PASS") {

    vcf <- read.vcfR(vcf_file, verbose = TRUE)
    sample_info_df <- INFO2df(vcf)
    df <- vcf@fix[, c("CHROM", "POS", "FILTER")]
    df <- cbind(df, sample_info_df[, c("SVTYPE", "END")])
    df <- cbind(df, sample_info_df[, add_sample_info_columns])
    names(df) <- c("CHROM", "POS", "FILTER", "SVTYPE", "END", add_sample_info_columns)

    if (filter != FALSE) {
        df <- df[df$FILTER == filter, ]
    }

    # Handle missing END
    missing_end <- is.na(df$END)
    df$END[missing_end] <- df$POS[missing_end]

    # Create g ranges
    gr <- GRanges(
        seqnames = df$CHROM,
        ranges = IRanges(start = as.numeric(df$POS), end = as.numeric(df$END)),
        strand = "*"
    )

    if (length(add_sample_info_columns) == 0) {
        gr$SVTYPE <- df[, "SVTYPE"]
    } else {
        mcols(gr) <- df[, c("SVTYPE", add_sample_info_columns)]
    }

    return(gr)
}


get_target_stats <- function(list_of_paths, name_after_path = c(1, 2, 4)) {
    stat_means <- data.frame(
        precision = numeric(),
        recall = numeric(),
        f1 = numeric(),
        data_name = character()
    )
    stat_data <- data.frame(
        features = character(),
        DEL_expected = numeric(),
        INS_expected = numeric(),
        NORMAL_expected = numeric(),
        precision = numeric(),
        recall = numeric(),
        f1 = numeric(),
        data_name = character()
    )

    for (path in list_of_paths) {
        
        ### Get target benchmark statistics
        # Get data name
        path_split <- unlist(strsplit(path, "/"))
        data_name <- paste0(path_split[name_after_path], collapse = "_")
        data_name <- gsub("\\.tsv", "", data_name)

        # Import all .tsv files
        tsv_data <- read.table(path, header = TRUE, sep = "\t")
        tsv_data$data_name <- data_name
        stat_data <- rbind(stat_data, tsv_data)

        # Calculate average precision, recall, and f1
        col_means <- colMeans(tsv_data[, c("precision", "recall", "f1")])
        col_means <- c(col_means, data_name = data_name)
        stat_means <- rbind(stat_means, t(data.frame(col_means)))
    }

    # Set the right data types
    stat_means[, 1] <- as.numeric(stat_means[, 1])
    stat_means[, 2] <- as.numeric(stat_means[, 2])
    stat_means[, 3] <- as.numeric(stat_means[, 3])

    # Manually convert to long format
    data_long <- data.frame(
        data_name = rep(stat_means$data_name, times = 3),
        Metric = rep(c("precision", "recall", "f1"), each = nrow(stat_means)),
        Value = c(stat_means$precision, stat_means$recall, stat_means$f1)
    )

    return(
        list(
            stat_means = stat_means,
            stat_data = stat_data,
            data_long = data_long
        )
    )
}


name_files_after_path <- function(files, name_after_path = c(1, 2, 4), file_extension = "json") {
    names <- c()

    for (path in files) {
        # Get data name
        path_split <- unlist(strsplit(path, "/"))
        data_name <- paste0(path_split[name_after_path], collapse = "_")
        data_name <- gsub(paste0("\\.", file_extension), "", data_name)

        names <- c(names, data_name)
    }

    names(files) <- names

    return(files)
}


get_truvari_stats <- function(x) {
    wgs_stats <- data.frame()

    for (n_path in 1:length(json_files)) {

        json_data <- read_json_stats(json_files[n_path])
        json_data$data_name <- names(json_files)[n_path]

        wgs_stats <- rbind(wgs_stats, json_data)
    }

    # Manually convert to long format
    data_long <- data.frame(
        data_name = rep(wgs_stats$data_name, times = 3),
        Metric = rep(c("precision", "recall", "f1"), each = nrow(wgs_stats)),
        Value = c(wgs_stats$precision, wgs_stats$recall, wgs_stats$f1)
    )

    return(
        list(
            wgs_stats = wgs_stats,
            data_long = data_long
        )
    )
}


read_json_stats <- function(x) {
    json_data <- fromJSON(x)
    # Replace all NULL values with NA
    json_data <- lapply(json_data, function(x) if (is.null(x)) NA else x)
    df_main <- as.data.frame(json_data[!names(json_data) %in% "gt_matrix"])
    
    return(df_main)
}
