#!/usr/bin/env Rscript

# Get script directory
script_dir <- tryCatch({
    dirname(sys.frame(1)$ofile)
}, error = function(e) {
    # Fallback: get directory from command line args
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    fileArg <- grep("^--file=", cmdArgs, value = TRUE)
    if (length(fileArg) > 0) {
        dirname(sub("^--file=", "", fileArg))
    } else {
        # Last fallback: assume functions.R is in same directory
        getwd()
    }
})

source(file.path(script_dir, "functions.R"))


run_dir <- getwd()
plots_dir <- "plots"
tables_dir <- "tables"

cat(sprintf("Run directory: %s\nPlots directory: %s\nTables directory: %s\n", 
            run_dir, plots_dir, tables_dir))

sapply(c(plots_dir, tables_dir), dir.create, showWarnings = FALSE, recursive = TRUE)

library(ggplot2)
library(jsonlite)
library(GenomicRanges)


# Load and process real intervals data
truvari_stats <- load_truvari_data(file.path(run_dir, "real_intervals"))

if (!is.null(truvari_stats)) {
    write.table(
        truvari_stats$raw_stats, file.path(tables_dir, "truvari_metrics_real_intervals.tsv"), 
        sep = "\t", row.names = FALSE
    )
} else {
    cat("Warning: No real intervals data found. Skipping real intervals analysis\n")
}

# Load and process simulated intervals data
truvari_stats_sim <- load_truvari_data(file.path(run_dir, "simulated_intervals/benchmarks"))

if (!is.null(truvari_stats_sim)) {
    # Write raw stats table
    write.table(
        truvari_stats_sim$raw_stats, 
        file.path(tables_dir, "truvari_metrics_simulated_intervals_raw.tsv"), 
        sep = "\t", row.names = FALSE
    )
} else {
    cat("Warning: No simulated intervals data found. Skipping simulated intervals analysis\n")
}

# Compare simulations with real data
if (!is.null(truvari_stats_sim$data_long) && !is.null(truvari_stats)) {
    # Helper function to calculate metrics for each tech/caller combination against a specific range
    calc_metrics <- function(tech, caller, range, metric) {
        # Get all simulation values for this tech/caller/metric combination (across all simulation runs)
        sim_vals <- truvari_stats_sim$data_long[
            truvari_stats_sim$data_long$tech == tech & 
            truvari_stats_sim$data_long$caller == caller & 
            truvari_stats_sim$data_long$Metric == metric, "Value"]
        
        # Get the real value for this tech/caller/range/metric combination
        real_val <- truvari_stats$data_long[
            truvari_stats$data_long$tech == tech & 
            truvari_stats$data_long$caller == caller & 
            truvari_stats$data_long$range == range & 
            truvari_stats$data_long$Metric == metric, "Value"]
        
        # Return NA if we don't have enough data
        if (length(sim_vals) < 2 || length(real_val) == 0) {
            return(c(percentile = NA, kde = NA, median = NA, diff = NA, sd = NA))
        }
        
        c(
            percentile = ecdf(sim_vals)(real_val) * 100,
            kde = kde(sim_vals, real_val),
            median = median(sim_vals),
            diff = median(sim_vals) - real_val,
            sd = sd(sim_vals)
        )
    }
    
    # Get unique combinations of tech/caller from simulations
    sim_combos <- unique(truvari_stats_sim$data_long[, c("tech", "caller")])
    
    # Get unique ranges from real data
    real_ranges <- unique(truvari_stats$data_long$range)
    
    # Create all combinations of sim tech/caller with real ranges
    all_combos <- expand.grid(
        tech = sim_combos$tech,
        caller = sim_combos$caller,
        range = real_ranges,
        stringsAsFactors = FALSE
    )
    
    match_stats <- t(sapply(seq_len(nrow(all_combos)), function(i) {
        tech <- all_combos$tech[i]
        caller <- all_combos$caller[i]
        range <- all_combos$range[i]
        
        metrics <- c("precision", "recall", "f1")
        stats <- lapply(metrics, function(m) calc_metrics(tech, caller, range, m))
        names(stats) <- metrics
        
        # Use unname() to remove the nested names from calc_metrics return values
        c(
            precision_percentile = unname(stats$precision[1]), 
            recall_percentile = unname(stats$recall[1]), 
            f1_percentile = unname(stats$f1[1]),
            precision_kde = unname(stats$precision[2]), 
            recall_kde = unname(stats$recall[2]), 
            f1_kde = unname(stats$f1[2]),
            precision_median = unname(stats$precision[3]), 
            recall_median = unname(stats$recall[3]), 
            f1_median = unname(stats$f1[3]),
            precision_diff = unname(stats$precision[4]), 
            recall_diff = unname(stats$recall[4]), 
            f1_diff = unname(stats$f1[4]),
            precision_sd = unname(stats$precision[5]), 
            recall_sd = unname(stats$recall[5]), 
            f1_sd = unname(stats$f1[5])
        )
    }))
    
    # Set rownames to tech-caller-range
    rownames(match_stats) <- paste(all_combos$tech, all_combos$caller, all_combos$range, sep = "-")
    
    # Filter out rows where all values are NA (non-existent combinations)
    valid_rows <- rowSums(!is.na(match_stats)) > 0
    match_stats <- match_stats[valid_rows, , drop = FALSE]
    
    write.table(
        format(data.frame(match_stats), scientific = FALSE), 
        file.path(tables_dir, "truvari_metrics_simulated_intervals.tsv"), 
        sep = "\t", row.names = TRUE
    )
} else {
    cat("Skipping simulation vs real comparison (missing data)\n")
    match_stats <- NULL
}

# Generate plots
# Handle different data availability scenarios
has_real <- !is.null(truvari_stats)
has_sim <- !is.null(truvari_stats_sim$data_long)
has_match <- !is.null(match_stats) && is.matrix(match_stats) && nrow(match_stats) > 0

if (!has_real && !has_sim) {
    cat("Skipping all plots: No data available\n")
} else {
    # Prepare data if available
    real_data <- if (has_real) annotate_data(truvari_stats$data_long) else NULL
    sim_data <- if (has_sim) annotate_data(truvari_stats_sim$data_long) else NULL
    
    # Create tech_caller columns
    if (!is.null(real_data)) {
        real_data$tech_caller <- paste(real_data$tech_clean, real_data$caller_display)
    }
    if (!is.null(sim_data)) {
        sim_data$tech_caller <- paste(sim_data$tech_clean, sim_data$caller)
    }
    
    # Get available colors (only for tech/caller combinations that exist in the data)
    available_techs <- character(0)
    if (!is.null(real_data)) available_techs <- c(available_techs, real_data$tech_caller)
    if (!is.null(sim_data)) available_techs <- c(available_techs, sim_data$tech_caller)
    available_techs <- unique(available_techs)
    
    plot_colors <- my_colors[names(my_colors) %in% available_techs]
    
    # Add missing colors with default gray
    missing_techs <- setdiff(available_techs, names(plot_colors))
    if (length(missing_techs) > 0) {
        new_colors <- setNames(rep("#808080", length(missing_techs)), missing_techs)
        plot_colors <- c(plot_colors, new_colors)
    }
    
    # Plot 1: Bar plot with simulation differences (only if match_stats has data)
    if (has_match) {
        tryCatch({
            # Extract diff and sd columns directly for wes_utr only
            wes_utr_rows <- grep("-wes_utr$", rownames(match_stats))
            
            if (length(wes_utr_rows) > 0) {
                diff_data <- data.frame(
                    tech_caller_raw = rownames(match_stats)[wes_utr_rows],
                    Precision = abs(match_stats[wes_utr_rows, "precision_diff"]),
                    Recall = abs(match_stats[wes_utr_rows, "recall_diff"]),
                    F1 = abs(match_stats[wes_utr_rows, "f1_diff"]),
                    Precision_sd = abs(match_stats[wes_utr_rows, "precision_sd"]),
                    Recall_sd = abs(match_stats[wes_utr_rows, "recall_sd"]),
                    F1_sd = abs(match_stats[wes_utr_rows, "f1_sd"]),
                    stringsAsFactors = FALSE
                )
                
                # Parse tech-caller-range format and create clean labels
                name_parts <- strsplit(diff_data$tech_caller_raw, "-")
                diff_data$tech <- sapply(name_parts, `[`, 1)
                diff_data$caller <- sapply(name_parts, `[`, 2)
                diff_data$range <- sapply(name_parts, `[`, 3)
                
                # Create a temporary dataframe to use annotate_data
                temp_df <- data.frame(
                    tech = diff_data$tech,
                    caller = diff_data$caller,
                    range = diff_data$range,
                    stringsAsFactors = FALSE
                )
                temp_annotated <- annotate_data(temp_df)
                diff_data$tech_caller <- paste(temp_annotated$tech_clean, temp_annotated$caller_display)
                
                # Reshape to long format
                long_df <- data.frame(
                    tech_caller = rep(diff_data$tech_caller, 3),
                    Metric = rep(c("Precision", "Recall", "F1"), each = nrow(diff_data)),
                    Value = c(diff_data$Precision, diff_data$Recall, diff_data$F1),
                    Value_sd = c(diff_data$Precision_sd, diff_data$Recall_sd, diff_data$F1_sd),
                    stringsAsFactors = FALSE
                )
                
                # Filter out NA values
                long_df <- long_df[!is.na(long_df$Value), ]
                
                if (nrow(long_df) > 0) {
                    ggsave(file.path(plots_dir, "bar_plot_sim_diff.png"),
                        ggplot(long_df, aes(x = tech_caller, y = Value, fill = tech_caller)) +
                        geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
                        geom_errorbar(
                            aes(ymin = Value - Value_sd, ymax = Value + Value_sd),
                            width = 0.2,
                            position = position_dodge(width = 0.8)
                        ) +
                        facet_grid(. ~ Metric, scales = "free_x", space = "free_x") +
                        scale_fill_manual(values = my_colors) +
                        theme_linedraw() +
                        theme(
                            strip.text = element_text(size = 10, face = "bold")
                        ) +
                        labs(
                            title = "Difference of metrics medians between simulated and real data (WES+UTR)",
                            y = "Median difference",
                            fill = "SV Caller"
                        ) +
                        guides(fill = guide_legend(title = "SV Caller")) +
                        scale_x_discrete(breaks = NULL) +
                        xlab(NULL),
                        width = 10, height = 6, units = "in", dpi = 300, bg = 'white')
                    cat("Plot 1 saved: bar_plot_sim_diff.png\n")
                } else {
                    cat("Skipping Plot 1: No valid data after filtering NAs\n")
                }
            } else {
                cat("Skipping Plot 1: No wes_utr data found\n")
            }
        }, error = function(e) {
            cat("Error creating Plot 1:", conditionMessage(e), "\n")
        })
    } else {
        cat("Skipping Plot 1: match_stats is empty\n")
    }
    
    # Plot 2: Boxplot for simulated and real data in EX+UTR region
    if (has_real && has_sim) {
        tryCatch({
            # Prepare real data subsets - check if ranges exist
            real_data_exutr <- data.frame()
            if ("EX+UTR" %in% real_data$ranges) {
                real_data_exutr <- rbind(real_data_exutr,
                    real_data[real_data$ranges == "EX+UTR" & real_data$caller_display != "WES Manta", ])
            }
            if ("HCI" %in% real_data$ranges) {
                real_data_exutr <- rbind(real_data_exutr,
                    real_data[real_data$ranges == "HCI" & real_data$caller_display != "WES Manta", ])
            }
        
        if (nrow(real_data_exutr) > 0 && nrow(sim_data) > 0) {
            real_data_exutr$tech_caller <- paste(real_data_exutr$tech_clean, real_data_exutr$caller)
            real_data_exutr$ranges_true <- "Simulated EX+UTR-like regions"
            
            ggsave(file.path(plots_dir, "facets_plot.png"),
                ggplot() +
                    geom_boxplot(data = sim_data, aes(x = factor(0), y = Value, fill = tech_caller),
                        outlier.size = 0.5, position = position_dodge(1), width = 0.7, alpha = 0.6) +
                    geom_point(data = real_data_exutr, 
                        aes(x = factor(0), y = Value, color = tech_caller, shape = ranges),
                        size = 2, position = position_dodge(1)) +
                    facet_grid(. ~ Metric, scales = "free_y", space = "free_x") +
                    scale_fill_manual(values = plot_colors, limits = force) +
                    scale_color_manual(values = plot_colors, limits = force) +
                    theme_linedraw() +
                    theme(strip.text = element_text(size = 10, face = "bold")) +
                    labs(title = "SV Metric Distributions in EX+UTR Region",
                        y = "Metric Value", shape = "Genome Range",
                        fill = "SV Caller (Simulated)", color = "SV Caller (True)") +
                    scale_x_discrete(breaks = NULL) + xlab(NULL),
                width = 10, height = 6, dpi = 300, bg = 'white')
            cat("Plot 2 saved: facets_plot.png\n")
        } else {
            cat("Skipping Plot 2: Insufficient data (need both real and simulated data)\n")
        }
        }, error = function(e) {
            cat("Error creating Plot 2:", conditionMessage(e), "\n")
        })
    } else {
        cat("Skipping Plot 2: Need both real and simulated data\n")
    }
    
    # Plot 3: Barplot by ranges (only needs real data)
    if (has_real) {
        tryCatch({
            real_data_wgs <- real_data[real_data$caller_display %in% 
                                    c("WES Manta", "WGS Manta", "Sniffles", "Pbsv", "CuteSV"), ]
        
        if (nrow(real_data_wgs) > 0) {
            real_data_wgs$tech_caller <- paste(real_data_wgs$tech_clean, real_data_wgs$caller_display)
            
            ggsave(file.path(plots_dir, "bar_plot.png"),
                ggplot(real_data_wgs, aes(x = ranges, y = Value, fill = tech_caller)) +
                    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
                    facet_grid(. ~ Metric, scales = "free_x", space = "free_x") +
                    scale_fill_manual(values = plot_colors, limits = force) +
                    theme_linedraw() +
                    theme(strip.text = element_text(size = 10, face = "bold")) +
                    labs(title = "Structural Variant Detection Metrics",
                        y = "Metric Value", fill = "SV Caller"),
                width = 10, height = 6, dpi = 300, bg = 'white')
            cat("Plot 3 saved: bar_plot.png\n")
        } else {
            cat("Skipping Plot 3: No real data for WGS/WES callers\n")
        }
        }, error = function(e) {
            cat("Error creating Plot 3:", conditionMessage(e), "\n")
        })
    } else {
        cat("Skipping Plot 3: No real data available\n")
    }
}

cat("\nAnalysis complete!\n")
