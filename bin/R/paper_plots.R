#!/usr/bin/env Rscript

# Get the first command-line argument
args <- commandArgs(trailingOnly = TRUE)

# Check if at least one argument was provided
if (length(args) == 0) {
    stop("No arguments provided. Please provide at least one argument.", call. = FALSE)
}

# Extract and use the first argument
run_name <- args[1]
cat("The run name is:", run_name, "\n")

library(ggplot2)
library(jsonlite)
library(GenomicRanges)

# Import needed functions
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

read_json_stats <- function(x) {
    json_data <- fromJSON(x)
    # Replace all NULL values with NA
    json_data <- lapply(json_data, function(x) if (is.null(x)) NA else x)
    df_main <- as.data.frame(json_data[!names(json_data) %in% "gt_matrix"])
    
    return(df_main)
}


my_colors <- c(
    "Pacbio CuteSV"         = "#ad2078",
    "Pacbio Pbsv"           = "#ea94c5ff",
    "ONT CuteSV"            = "#007ba5",
    "ONT Sniffles"          = "#9dcbe8ff",
    "Illumina WGS Manta"    = "#ffb441",
    "Illumina WES Manta"    = "#616161"
)


# Differences
get_mean_stats <- function(stats_all, throw_away = c()) {

    if (!is.null(throw_away)) {
        stats_all <- stats_all[!grepl(throw_away, stats_all$data_name), ]
    }

    high_confidence <- stats_all[grepl("high_confidence", stats_all$data_name), ]
    gene_panel <- stats_all[grepl("gene_panel", stats_all$data_name), ]
    wes_utr <- stats_all[grepl("wes_utr", stats_all$data_name), ]

    met <- unique(wes_utr$Metric)
    for (m in met) {
        m1 <- high_confidence$Value[high_confidence$Metric == m]
        m2 <- gene_panel$Value[gene_panel$Metric == m]
        m3 <- wes_utr$Value[wes_utr$Metric == m]


        print(paste0("Metric: ", m))
        print(paste0("High confidence: ", round(mean(m1), 2)))
        print(paste0("Gene panel: ", round(mean(m2), 2)))
        print(paste0("WES UTR: ", round(mean(m3), 2)))

        print(paste0("HCI vs. GP: ", round(mean(m1 - m2), 2)))
        print(paste0("HCI vs. WES+UTR: ", round(mean(m1 - m3), 2)))
        print(paste0("GP vs. WES+UTR: ", round(mean(m2 - m3), 2)))
    }
}



# Get json data from Truvari
# path_truvari <- paste0(paths, "/truvari_benchmark")
path_truvari <- paste0(run_name, "/real_intervals")
json_files <- list.files(path = path_truvari, pattern = "*summary.json$", full.names = TRUE, recursive = TRUE)
json_files <- name_files_after_path(json_files, name_after_path = c(3, 5, 6), file_extension = "json")
truvari_stats <- get_truvari_stats(json_files)


get_mean_stats(truvari_stats$data_long)
get_mean_stats(truvari_stats$data_long, throw_away = "Illumina_wes")

# save as tsv
write.table(truvari_stats$wgs_stats, paste0(run_name, "/stats/truvari_metrics_real_intervals.tsv"), sep = "\t", row.names = FALSE)


### Simulated targets analysis
path_truvari <- paste0(run_name, "/simulated_intervals/benchmarks")
json_files <- list.files(path = path_truvari, pattern = "*summary.json$", full.names = TRUE, recursive = TRUE)
json_files <- name_files_after_path(json_files, name_after_path = c(4, 5, 7, 8), file_extension = "json")
truvari_stats_sim <- get_truvari_stats(json_files)

# Remove all simulations with no calls
simulations_long <- truvari_stats_sim$data_long
invalid_sims <- unique(simulations_long$data_name[is.na(simulations_long$Value)])
simulations_long <- simulations_long[!simulations_long$data_name %in% invalid_sims, ]
simulations_long$data_name <- sub("^[^_]+_", "", simulations_long$data_name)

simulations_long <- simulations_long[!grepl("Illumina_wes", simulations_long$data_name), ]


kde <- function(data, test_point) {
    # Fit KDE
    kde <- density(
        data,
        from = min(data, test_point) - 0.5, 
        to = max(data, test_point) + 0.5
        )

    # Approximate probability density at test_point
    test_density <- approx(kde$x, kde$y, xout = test_point)$y

    # Convert density to "outlier probability" (relative to min/max density)
    max_density <- max(kde$y)
    outlier_prob <- 1 - (test_density / max_density)  # Closer to 1 = more outlier-like

    return(outlier_prob)
}


data_names <- unique(simulations_long$data_name)

match_stats <- sapply(data_names, function(data_name) {
    simulations_long_sub <- simulations_long[simulations_long$data_name == data_name, ]
    real_long_sub <- truvari_stats$data_long[truvari_stats$data_long$data_name == data_name, ]

    precision <- ecdf(simulations_long_sub[simulations_long_sub$Metric == "precision", "Value"])(real_long_sub[real_long_sub$Metric == "precision", "Value"]) * 100
    precision_kde <- kde(
        simulations_long_sub[simulations_long_sub$Metric == "precision", "Value"],
        real_long_sub[real_long_sub$Metric == "precision", "Value"]
    )

    recall <- ecdf(simulations_long_sub[simulations_long_sub$Metric == "recall", "Value"])(real_long_sub[real_long_sub$Metric == "recall", "Value"]) * 100
    recall_kde <- kde(
        simulations_long_sub[simulations_long_sub$Metric == "recall", "Value"],
        real_long_sub[real_long_sub$Metric == "recall", "Value"]
    )

    f1 <- ecdf(simulations_long_sub[simulations_long_sub$Metric == "f1", "Value"])(real_long_sub[real_long_sub$Metric == "f1", "Value"]) * 100
    f1_kde <- kde(
        simulations_long_sub[simulations_long_sub$Metric == "f1", "Value"],
        real_long_sub[real_long_sub$Metric == "f1", "Value"]
    )

    precision_med = median(simulations_long_sub[simulations_long_sub$Metric == "precision", "Value"])
    recall_med = median(simulations_long_sub[simulations_long_sub$Metric == "recall", "Value"])
    f1_med = median(simulations_long_sub[simulations_long_sub$Metric == "f1", "Value"])

    # Difference in medians
    precision_diff = precision_med - real_long_sub[real_long_sub$Metric == "precision", "Value"]
    recall_diff = recall_med - real_long_sub[real_long_sub$Metric == "recall", "Value"]
    f1_diff = f1_med - real_long_sub[real_long_sub$Metric == "f1", "Value"]

    c(
        precision_per = precision,
        recall_per = recall,
        f1_per = f1,
        precision_kde = precision_kde,
        recall_kde = recall_kde,
        f1_kde = f1_kde,
        precision_med = precision_med,
        recall_med = recall_med,
        f1_med = f1_med,
        precision_diff = precision_diff,
        recall_diff = recall_diff,
        f1_diff = f1_diff,
        precision_sd = sd(simulations_long_sub[simulations_long_sub$Metric == "precision", "Value"]),
        recall_sd = sd(simulations_long_sub[simulations_long_sub$Metric == "recall", "Value"]),
        f1_sd = sd(simulations_long_sub[simulations_long_sub$Metric == "f1", "Value"])
    )
})

match_stats <- data.frame(t(match_stats))

formatted_match_stats <- format(match_stats, scientific=F)

# save as tsv
write.table(formatted_match_stats, paste0(run_name, "/stats/truvari_metrics_simulated_intervals.tsv"), sep = "\t", row.names = TRUE)


# New plots:
annotate_data <- function(plot_data) {
    plot_data$tech <- ""
    plot_data$tech[grepl("Illumina", plot_data$data_name)] <- "Illumina"
    plot_data$tech[grepl("Pacbio", plot_data$data_name)] <- "Pacbio"
    plot_data$tech[grepl("ONT", plot_data$data_name)] <- "ONT"

    plot_data$ranges <- ""
    plot_data$ranges[grepl("high_confidence", plot_data$data_name)] <- "HCI"
    plot_data$ranges[grepl("gene_panel", plot_data$data_name)] <- "GP"
    plot_data$ranges[grepl("wes_utr", plot_data$data_name)] <- "EX+UTR"

    plot_data$caller <- ""
    plot_data$caller[grepl("wes_Manta", plot_data$data_name)] <- "WES Manta"
    plot_data$caller[grepl("wgs_Manta", plot_data$data_name)] <- "WGS Manta"
    plot_data$caller[grepl("CuteSV", plot_data$data_name)] <- "CuteSV"
    plot_data$caller[grepl("Sniffles", plot_data$data_name)] <- "Sniffles"
    plot_data$caller[grepl("Pbsv", plot_data$data_name)] <- "Pbsv"

    plot_data$ranges <- factor(plot_data$ranges, levels = c("HCI", "GP", "EX+UTR"))

    return(plot_data)
}

real_data <- annotate_data(truvari_stats$data_long)
sim_data <- annotate_data(simulations_long)


# Box plot
real_data_exutr <- real_data[real_data$ranges == "EX+UTR", ]
real_data_exutr <- real_data_exutr[real_data_exutr$caller != "WES Manta", ]
real_data_exutr$tech_caller <- paste(real_data_exutr$tech, real_data_exutr$caller, sep = " ")
real_data_exutr$ranges_true <- "Simulated EX+UTR-like regions"

real_data_hci <- real_data[real_data$ranges == "HCI" & real_data$caller != "WES Manta", ]
# real_data_hci$ranges <- "EX+UTR"
real_data_hci$tech_caller <- paste(real_data_hci$tech, real_data_hci$caller, sep = " ")
real_data_hci$ranges_true <- "Simulated EX+UTR-like regions"
real_data_exutr <- rbind(real_data_exutr, real_data_hci)


sim_data$tech_caller <- paste(sim_data$tech, sim_data$caller, sep = " ")

# Barplot simulation differences
formatted_match_stats2 <- formatted_match_stats
formatted_match_stats2$data_name <- rownames(formatted_match_stats2)
formatted_match_stats2 <- annotate_data(formatted_match_stats2)
formatted_match_stats2$tech_caller <- paste(formatted_match_stats2$tech, formatted_match_stats2$caller, sep = " ")
long_df <- reshape(
    formatted_match_stats2,
    varying = list(c("precision_diff", "recall_diff", "f1_diff")),
    v.names = "Value",
    timevar = "Metric",
    times = c("Precision", "Recall", "F1"),
    direction = "long"
)
long_df_sd <- reshape(
    formatted_match_stats2,
    varying = list(c("precision_sd", "recall_sd", "f1_sd")),
    v.names = "Value_sd",
    timevar = "Metric",
    times = c("Precision", "Recall", "F1"),
    direction = "long"
)
long_df$Value_sd <- as.numeric(long_df_sd$Value_sd)
rownames(long_df) <- NULL
long_df$Value <- as.numeric(long_df$Value)

# Create the bar plot
p11 <- ggplot(long_df, aes(x = tech_caller, y = Value, fill = tech_caller)) +
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
        title = "Difference of metrics medians between simulated and real data",
        y = "Median difference",
        fill = "SV Caller"
    ) +
    guides(fill = guide_legend(title = "SV Caller")) +
    scale_x_discrete(breaks = NULL) +
    xlab(NULL)

    # Save the plot
ggsave(
    paste0(run_name, "/stats/bar_plot_sim_diff.png"),
    p11, width = 10, height = 6, units = "in", dpi = 300, bg = 'white'
)


# Faced boxplot for simulated data and real data in EX+UTR region
p_dodge <- 1
p2 <- ggplot() +
    # Boxplot for simulated data
    geom_boxplot(
        data = sim_data,
        aes(x = factor(0), y = Value, fill = tech_caller),
        outlier.size = 0.5,
        position = position_dodge(width = p_dodge),
        width = 0.7,
        alpha = 0.6
    ) +
    # Points for real (true) data
    geom_point(
        data = real_data_exutr,
        aes(x = factor(0), y = Value, color = tech_caller, shape = ranges),
        size = 2,
        position = position_dodge(width = p_dodge)
    ) +
    facet_grid(. ~ Metric, scales = "free_y", space = "free_x") +
    scale_fill_manual(values = my_colors) +
    scale_color_manual(values = my_colors) +
    theme_linedraw() +
    theme(
        strip.text = element_text(size = 10, face = "bold")
    ) +
    labs(
        title = "SV Metric Distributions in EX+UTR Region",
        y = "Metric Value",
        shape = "Genome Range",
        fill = "SV Caller (Simulated)",
        color = "SV Caller (True)",
    ) +
    scale_x_discrete(breaks = NULL) +
    xlab(NULL)

ggsave(
    paste0(run_name, "/stats/facets_plot.png"),
    p2, width = 10, height = 6, units = "in", dpi = 300, bg = 'white'
    )


# Facet plot WGS???
real_data_wgs <- real_data[real_data$caller %in% c("WES Manta", "WGS Manta", "Sniffles", "Pbsv", "CuteSV"), ]
real_data_wgs$tech_caller <- paste(real_data_wgs$tech, real_data_wgs$caller, sep = " ")

p3 <- ggplot(real_data_wgs, aes(x = ranges, y = Value, fill = tech_caller)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    facet_grid(. ~ Metric, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = my_colors) +
    theme_linedraw() +
    theme(
        strip.text = element_text(size = 10, face = "bold")
    ) +
    labs(
        title = "Structural Variant Detection Metrics",
        y = "Metric Value",
        fill = "SV Caller"
    ) +
    guides(fill = guide_legend(title = "SV Caller"))

ggsave(
    paste0(run_name, "/stats/bar_plot.png"),
    p3, width = 10, height = 6, units = "in", dpi = 300, bg = 'white'
    )
