library(ggplot2)
library(patchwork)

# Define a consistent publication theme
theme_pub <- theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold")
  )

# Violin plot for TOTAL_READS
p1 <- ggplot(a, aes(x = svtype, y = TOTAL_READS, fill = svtype)) +
  geom_violin(trim = TRUE, scale = "width", color = "black") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  scale_y_log10() +
  labs(title = "Total Reads by SV Type", x = "SV Type", y = "Total Reads (log10)") +
  theme_pub +
  theme(legend.position = "none")

# Violin plot for AVG_READS_VALID_SAMPLE (log10 scale)
p2 <- ggplot(a, aes(x = svtype, y = AVG_READS_VALID_SAMPLE, fill = svtype)) +
  geom_violin(trim = TRUE, scale = "width", color = "black") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  scale_y_log10() +
  labs(title = "Average Reads per Sample", x = "SV Type", y = "Avg. Reads per Valid Sample (log10)") +
  theme_pub +
  theme(legend.position = "none")

# Violin plot for AVG_MAP_QUALITY
p3 <- ggplot(a, aes(x = svtype, y = AVG_MAP_QUALITY, fill = svtype)) +
  geom_violin(trim = TRUE, scale = "width", color = "black") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  labs(title = "Average Mapping Quality", x = "SV Type", y = "MAPQ") +
  theme_pub +
  theme(legend.position = "none")

# Violin plot for LD_SNP_RATE
p4 <- ggplot(a, aes(x = svtype, y = LD_SNP_RATE, fill = svtype)) +
  geom_violin(trim = TRUE, scale = "width", color = "black") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  labs(title = "LD SNP Rate", x = "SV Type", y = "LD SNP Rate") +
  theme_pub +
  theme(legend.position = "none")

# Violin plot for MAF (this is not MAF. It is a fraction of carriers (het or hom_alt) out of all valid genotype calls)
p5 <- ggplot(a, aes(x = svtype, y = MAF, fill = svtype)) +
  geom_violin(trim = TRUE, scale = "width", color = "black") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  labs(title = "Fraction of ALT allele carriers", x = "SV Type", y = "MAF") +
  theme_pub +
  theme(legend.position = "none")

# Violin plot for MISSING_RATE
p6 <- ggplot(a, aes(x = svtype, y = MISSING_RATE, fill = svtype)) +
  geom_violin(trim = TRUE, scale = "width", color = "black") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  labs(title = "Missing Rate", x = "SV Type", y = "Missing Rate") +
  theme_pub + theme(legend.position = "none")


# Combine all four plots
combined_plot <- (p1 + p2 + p3) / (p4 + p5 + p6)


# Save the final figure in high resolution
ggsave("SV_metrics_violin_plots.png", combined_plot, width = 12, height = 8, dpi = 300)
#ggsave("SV_metrics_violin_plots_supp.png", supp_plot, width = 8, height = 4, dpi = 300)
















library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

a = read.table("table_SV_QC.tsv", sep='\t', header=1)

vars_to_plot <- c(
  "TOTAL_READS", "AVG_READS_VALID_SAMPLE", "AVG_MAP_QUALITY",
  "LD_SNP_RATE", "MAF", "MISSING_RATE"
)

a_long <- a %>%
  select(svtype, all_of(vars_to_plot)) %>%
  pivot_longer(
    cols = all_of(vars_to_plot),
    names_to = "metric",
    values_to = "value"
  )

log_vars <- c("TOTAL_READS", "AVG_READS_VALID_SAMPLE")

a_long <- a_long %>%
  mutate(
    value_transformed = ifelse(metric %in% log_vars & value > 0, log10(value), value),
    metric_label = case_when(
      metric == "TOTAL_READS" ~ "Total Reads (log10)",
      metric == "AVG_READS_VALID_SAMPLE" ~ "Avg. Reads Valid GT (log10)",
      metric == "AVG_MAP_QUALITY" ~ "Average Mapping Quality",
      metric == "LD_SNP_RATE" ~ "LD SNP Rate",
      metric == "MAF" ~ "ALT Frequency in Valid GT",
      metric == "MISSING_RATE" ~ "Missing Rate",
      TRUE ~ metric
    )
  )

# Define a publication-quality ggplot theme
theme_pub <- theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold")
  )

# Create faceted histograms by metric
hist_plot <- ggplot(a_long, aes(x = value_transformed, fill = svtype)) +
  geom_histogram(alpha = 0.7, bins = 50, position = "identity") +
  facet_wrap(~metric_label, scales = "free", ncol = 3) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Distribution of SV Metrics by Type",
    x = "Value",
    y = "Count",
    fill = "SV Type"
  ) +
  theme_pub

# Save high-resolution figure
ggsave("SV_metrics_histograms.pdf", hist_plot, width = 12, height = 8, units='in')

# Display the combined figure
hist_plot

