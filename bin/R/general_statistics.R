#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript general_statistics.R <run_dir> <output_file>", call. = FALSE)

run_dir <- args[1]
output_file <- args[2]

cat(sprintf("Run directory: %s\nOutput file: %s\n", run_dir, output_file))

library(ggplot2)
library(jsonlite)
library(GenomicRanges)

source(file.path(dirname(sys.frame(1)$ofile), "functions.R"))

# Helper function to add summary lines
add_summary <- function(lines, ...) c(lines, sprintf(...))

# Initialize summary
summary_lines <- c(
    "===============================================",
    "  GENERAL STATISTICS SUMMARY",
    "===============================================",
    "",
    sprintf("Generated on: %s", Sys.time()),
    sprintf("Run directory: %s", run_dir),
    ""
)

# Analyze real intervals
summary_lines <- c(summary_lines, "=== Real Intervals Analysis ===")
path_real <- file.path(run_dir, "real_intervals")

if (dir.exists(path_real)) {
    truvari_stats <- load_truvari_data(path_real)
    
    if (!is.null(truvari_stats) && !is.null(truvari_stats$data_long)) {
        summary_lines <- add_summary(summary_lines, "Found %d benchmark results", 
                                    length(unique(truvari_stats$data_long$data_name)))
        
        # Get unique data names and their metrics
        unique_names <- unique(truvari_stats$data_long$data_name)
        for (dn in unique_names) {
            subset_data <- truvari_stats$data_long[truvari_stats$data_long$data_name == dn, ]
            prec <- subset_data$Value[subset_data$Metric == "precision"]
            rec <- subset_data$Value[subset_data$Metric == "recall"]
            f1 <- subset_data$Value[subset_data$Metric == "f1"]
            
            summary_lines <- c(summary_lines, "",
                             sprintf("  %s", dn),
                             sprintf("    Precision: %.4f", ifelse(length(prec) > 0, prec, NA)),
                             sprintf("    Recall:    %.4f", ifelse(length(rec) > 0, rec, NA)),
                             sprintf("    F1 Score:  %.4f", ifelse(length(f1) > 0, f1, NA)))
        }
    } else {
        summary_lines <- c(summary_lines, "  No valid benchmark results found")
    }
} else {
    summary_lines <- c(summary_lines, "  Directory not found")
}
summary_lines <- c(summary_lines, "")

# Analyze simulated intervals
summary_lines <- c(summary_lines, "=== Simulated Intervals Analysis ===")
path_sim <- file.path(run_dir, "simulated_intervals", "benchmarks")

if (dir.exists(path_sim)) {
    sim_stats <- load_truvari_data(path_sim)
    
    if (!is.null(sim_stats) && !is.null(sim_stats$data_long)) {
        summary_lines <- add_summary(summary_lines, "Found %d simulated benchmark results", 
                                    length(unique(sim_stats$data_long$data_name)))
        
        prec_vals <- sim_stats$data_long$Value[sim_stats$data_long$Metric == "precision"]
        rec_vals <- sim_stats$data_long$Value[sim_stats$data_long$Metric == "recall"]
        f1_vals <- sim_stats$data_long$Value[sim_stats$data_long$Metric == "f1"]
        
        summary_lines <- c(summary_lines, "",
                          "  Mean Metrics Across All Simulations:",
                          sprintf("    Precision: %.4f", mean(prec_vals, na.rm = TRUE)),
                          sprintf("    Recall:    %.4f", mean(rec_vals, na.rm = TRUE)),
                          sprintf("    F1 Score:  %.4f", mean(f1_vals, na.rm = TRUE)))
    } else {
        summary_lines <- c(summary_lines, "  No valid simulated benchmark results found")
    }
} else {
    summary_lines <- c(summary_lines, "  Directory not found")
}

summary_lines <- c(summary_lines, "",
                  "===============================================",
                  "  END OF SUMMARY",
                  "===============================================")

writeLines(summary_lines, output_file)
cat(paste(summary_lines, collapse = "\n"), "\n\nSummary written to:", output_file, "\n")
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
