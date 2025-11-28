#!/usr/bin/env Rscript

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments was provided
if (length(args) < 2) {
    stop("Usage: Rscript general_statistics.R <run_dir> <output_file>", call. = FALSE)
}

# Extract arguments
run_dir <- args[1]
output_file <- args[2]

cat("Run directory:", run_dir, "\n")
cat("Output file:", output_file, "\n")

library(ggplot2)
library(jsonlite)
library(GenomicRanges)

# Import needed functions
read_json_stats <- function(x) {
    json_data <- fromJSON(x)
    # Replace all NULL values with NA
    json_data <- lapply(json_data, function(x) if (is.null(x)) NA else x)
    df_main <- as.data.frame(json_data[!names(json_data) %in% "gt_matrix"])
    
    return(df_main)
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

get_truvari_stats <- function(json_files) {
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


plot_matrices <- function(df, plot_path, return_plot = FALSE) {
    # Create the barplot
    p1 <- ggplot(df, aes(x = data_name, y = Value, fill = Metric)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme_minimal() +
        labs(
            title = "Performance Metrics by Data Name",
            x = "Data Name",
            y = "Metric Value",
            fill = "Metric"
        ) +
        ylim(0, 1) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # Adjust rotation and alignment
            plot.title = element_text(hjust = 0.5) # Center the title
        ) +
        geom_text(
            aes(label = round(Value, 2)), 
            color = "black", 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, # Adjust text position above bars
            size = 3
        )

    # Save plot
    ggsave(
        plot_path,
        p1, width = 10, height = 6, units = "in", dpi = 300, bg = 'white'
        )

    if (return_plot) {
        return(p1)
    }
}


# Initialize summary output
summary_lines <- c()
summary_lines <- c(summary_lines, "===============================================")
summary_lines <- c(summary_lines, "  GENERAL STATISTICS SUMMARY")
summary_lines <- c(summary_lines, "===============================================")
summary_lines <- c(summary_lines, "")
summary_lines <- c(summary_lines, paste("Generated on:", Sys.time()))
summary_lines <- c(summary_lines, paste("Run directory:", run_dir))
summary_lines <- c(summary_lines, "")

# Try to find and analyze real intervals data
path_real <- file.path(run_dir, "real_intervals")
if (dir.exists(path_real)) {
    summary_lines <- c(summary_lines, "=== Real Intervals Analysis ===")
    
    json_files <- list.files(path = path_real, pattern = "*summary.json$", 
                             full.names = TRUE, recursive = TRUE)
    
    if (length(json_files) > 0) {
        summary_lines <- c(summary_lines, paste("Found", length(json_files), "benchmark results"))
        
        # Parse file names to extract technology and caller info
        json_files_named <- name_files_after_path(json_files, 
                                                   name_after_path = c(3, 5, 6), 
                                                   file_extension = "json")
        truvari_stats <- get_truvari_stats(json_files_named)
        
        # Summary statistics by technology/caller
        for (i in 1:nrow(truvari_stats$wgs_stats)) {
            data_name <- truvari_stats$wgs_stats$data_name[i]
            precision <- round(truvari_stats$wgs_stats$precision[i], 4)
            recall <- round(truvari_stats$wgs_stats$recall[i], 4)
            f1 <- round(truvari_stats$wgs_stats$f1[i], 4)
            
            summary_lines <- c(summary_lines, "")
            summary_lines <- c(summary_lines, paste("  ", data_name))
            summary_lines <- c(summary_lines, paste("    Precision:", precision))
            summary_lines <- c(summary_lines, paste("    Recall:   ", recall))
            summary_lines <- c(summary_lines, paste("    F1 Score: ", f1))
        }
    } else {
        summary_lines <- c(summary_lines, "  No benchmark results found")
    }
    summary_lines <- c(summary_lines, "")
} else {
    summary_lines <- c(summary_lines, "=== Real Intervals Analysis ===")
    summary_lines <- c(summary_lines, "  Directory not found")
    summary_lines <- c(summary_lines, "")
}

# Try to find and analyze simulated intervals data
path_sim <- file.path(run_dir, "simulated_intervals", "benchmarks")
if (dir.exists(path_sim)) {
    summary_lines <- c(summary_lines, "=== Simulated Intervals Analysis ===")
    
    json_files <- list.files(path = path_sim, pattern = "*summary.json$", 
                             full.names = TRUE, recursive = TRUE)
    
    if (length(json_files) > 0) {
        summary_lines <- c(summary_lines, paste("Found", length(json_files), "simulated benchmark results"))
        
        # Parse and get basic statistics
        json_files_named <- name_files_after_path(json_files, 
                                                   name_after_path = c(4, 5, 7, 8), 
                                                   file_extension = "json")
        sim_stats <- get_truvari_stats(json_files_named)
        
        # Calculate mean metrics across all simulations
        mean_precision <- mean(sim_stats$wgs_stats$precision, na.rm = TRUE)
        mean_recall <- mean(sim_stats$wgs_stats$recall, na.rm = TRUE)
        mean_f1 <- mean(sim_stats$wgs_stats$f1, na.rm = TRUE)
        
        summary_lines <- c(summary_lines, "")
        summary_lines <- c(summary_lines, "  Mean Metrics Across All Simulations:")
        summary_lines <- c(summary_lines, paste("    Precision:", round(mean_precision, 4)))
        summary_lines <- c(summary_lines, paste("    Recall:   ", round(mean_recall, 4)))
        summary_lines <- c(summary_lines, paste("    F1 Score: ", round(mean_f1, 4)))
    } else {
        summary_lines <- c(summary_lines, "  No simulated benchmark results found")
    }
    summary_lines <- c(summary_lines, "")
} else {
    summary_lines <- c(summary_lines, "=== Simulated Intervals Analysis ===")
    summary_lines <- c(summary_lines, "  Directory not found")
    summary_lines <- c(summary_lines, "")
}

summary_lines <- c(summary_lines, "===============================================")
summary_lines <- c(summary_lines, "  END OF SUMMARY")
summary_lines <- c(summary_lines, "===============================================")

# Write summary to file
writeLines(summary_lines, output_file)

# Also print to console
cat(paste(summary_lines, collapse = "\n"))
cat("\n\nSummary written to:", output_file, "\n")
