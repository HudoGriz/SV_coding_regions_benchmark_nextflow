# Shared functions for R analysis scripts
# This file contains common functions used by paper_plots.R and general_statistics.R

# Plot colors
my_colors <- c(
    "PacBio CuteSV" = "#ad2078", "PacBio Pbsv" = "#ea94c5ff",
    "ONT CuteSV" = "#007ba5", "ONT Sniffles" = "#9dcbe8ff",
    "Illumina WGS Manta" = "#ffb441", "Illumina WES Manta" = "#616161"
)

# Load and process real intervals data
load_truvari_data <- function(path) {
    # Check if path exists
    if (!dir.exists(path)) {
        cat(sprintf("Warning: Directory does not exist: %s\n", path))
        return(NULL)
    }
    
    json_files <- list.files(path, pattern = "*summary.json$", full.names = TRUE, recursive = TRUE)
    
    if (length(json_files) == 0) {
        cat(sprintf("Warning: No summary.json files found in: %s\n", path))
        return(NULL)
    }
    
    cat(sprintf("Found %d summary.json file(s) in: %s\n", length(json_files), path))
    
    named_files <- name_files_after_path(json_files, file_extension = "json")
    get_truvari_stats(named_files)
}

# Read JSON statistics from Truvari output
read_json_stats <- function(x) {
    tryCatch({
        json_data <- fromJSON(x)
        json_data[sapply(json_data, is.null)] <- NA  # Simplified NULL handling
        as.data.frame(json_data[!names(json_data) %in% "gt_matrix"])
    }, error = function(e) {
        cat(sprintf("Warning: Failed to read JSON file %s: %s\n", x, conditionMessage(e)))
        return(NULL)
    })
}

# Name files based on path components
name_files_after_path <- function(files, file_extension = "json") {
    names_res <- sapply(files, function(path) {
        bn <- basename(path)
        bn2 <- sub(paste0("\\.summary\\.", file_extension, "$"), "", bn)
        bn2 <- sub(paste0("\\.", file_extension, "$"), "", bn2)
        bn2
    }, USE.NAMES = FALSE)

    names(files) <- names_res
    files
}

# Get Truvari statistics from JSON files
get_truvari_stats <- function(json_files) {
    # Read all JSON files and filter out failed reads
    json_data_list <- lapply(seq_along(json_files), function(i) {
        json_data <- read_json_stats(json_files[i])
        if (!is.null(json_data)) {
            json_data$data_name <- names(json_files)[i]
        }
        json_data
    })
    
    # Remove NULL entries (failed reads)
    json_data_list <- json_data_list[!sapply(json_data_list, is.null)]
    
    if (length(json_data_list) == 0) {
        cat("Warning: No valid JSON files could be read\n")
        return(NULL)
    }
    
    raw_stats <- do.call(rbind, json_data_list)
    
    # Parse data_name into tech, caller, and range components
    # Handle cases where data_name might not follow expected format
    name_parts <- strsplit(raw_stats$data_name, "-")
    raw_stats$tech <- sapply(name_parts, function(x) ifelse(length(x) >= 1, x[1], NA))
    raw_stats$caller <- sapply(name_parts, function(x) ifelse(length(x) >= 2, x[2], NA))
    raw_stats$range <- sapply(name_parts, function(x) ifelse(length(x) >= 3, x[3], NA))
    
    # Warn about any missing components
    if (any(is.na(raw_stats$tech))) {
        cat("Warning: Some files have missing 'tech' component in their names\n")
    }
    if (any(is.na(raw_stats$caller))) {
        cat("Warning: Some files have missing 'caller' component in their names\n")
    }
    if (any(is.na(raw_stats$range))) {
        cat("Warning: Some files have missing 'range' component in their names\n")
    }
    
    # Convert to long format with parsed columns
    data_long <- data.frame(
        data_name = rep(raw_stats$data_name, times = 3),
        tech = rep(raw_stats$tech, times = 3),
        caller = rep(raw_stats$caller, times = 3),
        range = rep(raw_stats$range, times = 3),
        Metric = rep(c("precision", "recall", "f1"), each = nrow(raw_stats)),
        Value = c(raw_stats$precision, raw_stats$recall, raw_stats$f1)
    )
    
    return(
        list(
            raw_stats = raw_stats,
            data_long = data_long
        )
    )
}

# Get mean statistics and print differences
get_mean_stats <- function(stats_all, throw_away = NULL) {
    if (is.null(stats_all) || nrow(stats_all) == 0) {
        cat("Warning: No data available for statistics\n")
        return(invisible(NULL))
    }
    
    if (!is.null(throw_away)) {
        stats_all <- stats_all[!grepl(throw_away, stats_all$data_name), ]
        if (nrow(stats_all) == 0) {
            cat(sprintf("Warning: No data remaining after filtering out '%s'\n", throw_away))
            return(invisible(NULL))
        }
    }

    # Create subsets, handling missing ranges gracefully
    subsets <- list(
        high_confidence = stats_all[grepl("high_confidence", stats_all$data_name), ],
        gene_panel = stats_all[grepl("gene_panel", stats_all$data_name), ],
        wes_utr = stats_all[grepl("wes_utr", stats_all$data_name), ]
    )
    
    # Check which subsets have data
    available_subsets <- names(subsets)[sapply(subsets, nrow) > 0]
    
    if (length(available_subsets) == 0) {
        cat("Warning: No data found for any range (high_confidence, gene_panel, wes_utr)\n")
        return(invisible(NULL))
    }
    
    # Get unique metrics from available data
    metrics_available <- unique(stats_all$Metric)
    if (length(metrics_available) == 0) {
        cat("Warning: No metrics found in data\n")
        return(invisible(NULL))
    }

    for (m in metrics_available) {
        vals <- lapply(subsets, function(s) {
            if (nrow(s) > 0) {
                s$Value[s$Metric == m]
            } else {
                numeric(0)
            }
        })
        
        cat(sprintf("\nMetric: %s\n", m))
        
        # Print available statistics
        if (length(vals$high_confidence) > 0) {
            cat(sprintf("  High confidence: %.2f\n", mean(vals$high_confidence, na.rm = TRUE)))
        } else {
            cat("  High confidence: N/A (no data)\n")
        }
        
        if (length(vals$gene_panel) > 0) {
            cat(sprintf("  Gene panel: %.2f\n", mean(vals$gene_panel, na.rm = TRUE)))
        } else {
            cat("  Gene panel: N/A (no data)\n")
        }
        
        if (length(vals$wes_utr) > 0) {
            cat(sprintf("  WES UTR: %.2f\n", mean(vals$wes_utr, na.rm = TRUE)))
        } else {
            cat("  WES UTR: N/A (no data)\n")
        }
        
        # Print differences only if both datasets are available
        if (length(vals$high_confidence) > 0 && length(vals$gene_panel) > 0) {
            cat(sprintf("  HCI vs. GP: %.2f\n", 
                mean(vals$high_confidence - vals$gene_panel, na.rm = TRUE)))
        }
        
        if (length(vals$high_confidence) > 0 && length(vals$wes_utr) > 0) {
            cat(sprintf("  HCI vs. WES+UTR: %.2f\n", 
                mean(vals$high_confidence - vals$wes_utr, na.rm = TRUE)))
        }
        
        if (length(vals$gene_panel) > 0 && length(vals$wes_utr) > 0) {
            cat(sprintf("  GP vs. WES+UTR: %.2f\n", 
                mean(vals$gene_panel - vals$wes_utr, na.rm = TRUE)))
        }
    }
}

# Kernel Density Estimation for outlier detection
kde <- function(data, test_point) {
    kde_fit <- density(
        data, from = min(data, test_point) - 0.5, 
        to = max(data, test_point) + 0.5
    )
    test_density <- approx(kde_fit$x, kde_fit$y, xout = test_point)$y
    1 - (test_density / max(kde_fit$y))  # Outlier probability
}

# Annotate data with technology, ranges, and caller information
# Note: tech, caller, and range columns should already exist from get_truvari_stats
# This function now just maps range values to display labels and creates factor levels
annotate_data <- function(plot_data) {
    if (is.null(plot_data) || nrow(plot_data) == 0) {
        cat("Warning: annotate_data received empty data\n")
        return(plot_data)
    }
    
    # Check if required columns exist
    required_cols <- c("range", "caller", "tech")
    missing_cols <- setdiff(required_cols, names(plot_data))
    if (length(missing_cols) > 0) {
        cat(sprintf("Warning: Missing columns in data: %s\n", paste(missing_cols, collapse = ", ")))
        # Add missing columns with NA
        for (col in missing_cols) {
            plot_data[[col]] <- NA
        }
    }
    
    # Map range values to display labels
    plot_data$ranges <- ifelse(plot_data$range == "high_confidence", "HCI",
        ifelse(plot_data$range == "gene_panel", "GP",
        ifelse(plot_data$range == "wes_utr", "EX+UTR", 
        ifelse(is.na(plot_data$range), NA, plot_data$range))))
    
    # Clean tech names (remove _WGS and _WES suffixes)
    plot_data$tech_clean <- gsub("_WGS|_WES", "", plot_data$tech)
    
    # Map caller values to display labels (handle WES vs WGS Manta)
    plot_data$caller_display <- ifelse(
        plot_data$caller == "Manta" & grepl("_WES", plot_data$tech), "WES Manta",
        ifelse(plot_data$caller == "Manta" & grepl("_WGS", plot_data$tech), "WGS Manta",
        ifelse(is.na(plot_data$caller), NA, plot_data$caller)))
    
    # Create factor with available levels
    available_levels <- intersect(c("HCI", "GP", "EX+UTR"), unique(plot_data$ranges))
    if (length(available_levels) > 0) {
        plot_data$ranges <- factor(plot_data$ranges, levels = c("HCI", "GP", "EX+UTR"))
    }
    
    plot_data
}

# Get target statistics from TSV files
get_target_stats <- function(list_of_paths, name_after_path = c(1, 2, 4)) {
    if (length(list_of_paths) == 0) {
        cat("Warning: No paths provided to get_target_stats\n")
        return(NULL)
    }
    
    stat_list <- lapply(list_of_paths, function(path) {
        tryCatch({
            if (!file.exists(path)) {
                cat(sprintf("Warning: File does not exist: %s\n", path))
                return(NULL)
            }
            
            path_split <- strsplit(path, "/")[[1]]
            data_name <- gsub("\\.tsv", "", paste(path_split[name_after_path], collapse = "_"))
            
            tsv_data <- read.table(path, header = TRUE, sep = "\t")
            
            # Check if required columns exist
            required_cols <- c("precision", "recall", "f1")
            missing_cols <- setdiff(required_cols, names(tsv_data))
            if (length(missing_cols) > 0) {
                cat(sprintf("Warning: Missing columns in %s: %s\n", 
                    path, paste(missing_cols, collapse = ", ")))
                return(NULL)
            }
            
            tsv_data$data_name <- data_name
            
            list(
                data = tsv_data,
                means = c(
                    colMeans(tsv_data[, c("precision", "recall", "f1")], na.rm = TRUE), 
                    data_name = data_name
                )
            )
        }, error = function(e) {
            cat(sprintf("Error reading %s: %s\n", path, conditionMessage(e)))
            return(NULL)
        })
    })
    
    # Remove NULL entries (failed reads)
    stat_list <- stat_list[!sapply(stat_list, is.null)]
    
    if (length(stat_list) == 0) {
        cat("Warning: No valid TSV files could be read\n")
        return(NULL)
    }
    
    stat_data <- do.call(rbind, lapply(stat_list, `[[`, "data"))
    stat_means <- as.data.frame(do.call(rbind, lapply(stat_list, `[[`, "means")))
    stat_means[, 1:3] <- lapply(stat_means[, 1:3], as.numeric)
    
    list(
        stat_means = stat_means,
        stat_data = stat_data,
        data_long = data.frame(
            data_name = rep(stat_means$data_name, times = 3),
            Metric = rep(c("precision", "recall", "f1"), each = nrow(stat_means)),
            Value = c(stat_means$precision, stat_means$recall, stat_means$f1)
        )
    )
}

# Plot matrices as barplot
plot_matrices <- function(df, plot_path, return_plot = FALSE) {
    if (is.null(df) || nrow(df) == 0) {
        cat("Warning: No data available for plotting\n")
        return(invisible(NULL))
    }
    
    tryCatch({
        p1 <- ggplot(df, aes(x = data_name, y = Value, fill = Metric)) +
            geom_bar(stat = "identity", position = "dodge") +
            geom_text(aes(label = round(Value, 2)), 
                position = position_dodge(0.9), vjust = -0.5, size = 3) +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                plot.title = element_text(hjust = 0.5)) +
            labs(title = "Performance Metrics by Data Name",
                x = "Data Name", y = "Metric Value", fill = "Metric") +
            ylim(0, 1)

        ggsave(plot_path, p1, width = 10, height = 6, units = "in", dpi = 300, bg = 'white')
        cat(sprintf("Plot saved: %s\n", plot_path))
        
        if (return_plot) return(p1)
    }, error = function(e) {
        cat(sprintf("Error creating plot: %s\n", conditionMessage(e)))
        return(invisible(NULL))
    })
}
