source("scripts/R/general_functions.R")

library(ggplot2)


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


paths <- c(
    "Illumina_wes/exomes",
    "Illumina_wes/exomes_filtered",
    "Illumina_wes/utr_exomes",
    "Illumina_wes/utr_exomes_filtered",
    "Illumina_wgs/exomes_filtered",
    "Illumina_wgs/utr_exomes",
    "Pacbio/exomes_filtered",
    "Pacbio/utr_exomes",
    "ONT/exomes_filtered",
    "ONT/utr_exomes",
    "Illumina_wes/exome_filtered_on_wgs_data"
)

path_target <- paste0(paths, "/target_benchmark")
tsv_files <- list.files(path = path_target, pattern = "\\.tsv$", full.names = TRUE, recursive = TRUE)

target_metrics <- get_target_stats(tsv_files)

data_long <- target_metrics$data_long

data_long_utr <- data_long[grepl("utr", data_long$data_name), ]
data_long_utr <- data_long_utr[!grepl("Illumina_wes_utr_exomes_Manta", data_long_utr$data_name), ]
plot_matrices(data_long_utr, "plots/target_metrics_utr.png")

data_long_exome <- data_long[!grepl("utr", data_long$data_name), ]
data_long_exome <- data_long_exome[!grepl("Illumina_wes_exomes_Manta", data_long_exome$data_name), ]
plot_exome <- plot_matrices(data_long_exome, "plots/target_metrics_exons.png", return_plot = TRUE)


# Get json data from Truvari
path_truvari <- paste0(paths, "/truvari_benchmark")
json_files <- list.files(path = path_truvari, pattern = "*summary.json$", full.names = TRUE, recursive = TRUE)

wgs_stats <- get_wgs_stats(json_files)

wgs_stats_long <- wgs_stats$data_long

# It's the same 
# wgs_stats_long_utr <- wgs_stats_long[grepl("utr", wgs_stats_long$data_name), ]
# plot_matrices(wgs_stats_long_utr, "plots/wgs_metrics_utr_wgs.png")

wgs_stats_long_exome <- wgs_stats_long[!grepl("utr", wgs_stats_long$data_name), ]
plot_matrices(wgs_stats_long_exome, "plots/wgs_metrics_exons.png")





# Simulated targets analysis
sim_path <- "simulated_targets_benchmark/target_benchmark"
sim_files <- list.files(path = sim_path, pattern = "\\.tsv$", full.names = TRUE, recursive = TRUE)

sim_metrics <- get_target_stats(sim_files, name_after_path = c(3))
sim_data_long <- sim_metrics$data_long


split_data <- do.call(rbind, strsplit(sim_data_long$data_name, "_simulation"))
sim_data_long$data_name <- split_data[, 1]
sim_data_long$N <- as.integer(split_data[, 2])

head(sim_data_long)

df <- sim_data_long

unique(sim_data_long$data_name)

data_long_exome_sim <- data_long_exome[!grepl("Illumina_wes_exomes_filtered_Manta", data_long_exome$data_name), ]
data_long_exome_sim <- data_long_exome_sim[!grepl("Illumina_wes_exome_filtered_on_wgs_data_Manta", data_long_exome_sim$data_name), ]
data_long_exome_sim$data_name <- gsub("Illumina_wgs_exomes_filtered_Manta", "Illumina_Manta", data_long_exome_sim$data_name)
data_long_exome_sim$data_name <- gsub("ONT_exomes_filtered_CuteSV", "ONT_CuteSV", data_long_exome_sim$data_name)
data_long_exome_sim$data_name <- gsub("ONT_exomes_filtered_Sniffles", "ONT_Sniffles", data_long_exome_sim$data_name)
data_long_exome_sim$data_name <- gsub("Pacbio_exomes_filtered_CuteSV", "Pacbio_CuteSV", data_long_exome_sim$data_name)
data_long_exome_sim$data_name <- gsub("Pacbio_exomes_filtered_Pbsv", "Pacbio_Pbsv", data_long_exome_sim$data_name)

plot_sim <- ggplot(data_long_exome_sim, aes(x = data_name, y = Value, fill = Metric)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_violin(data = sim_data_long, aes(x = data_name, y = Value, fill = Metric),
            position = position_dodge(width = 0.9), alpha = 0.5) +
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


plot_path <- "plots/sim_metrics.png"
# Save plot
ggsave(
    plot_path,
    plot_sim, width = 10, height = 6, units = "in", dpi = 300, bg = 'white'
    )
