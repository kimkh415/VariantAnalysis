library(qs)
library(dplyr)
library(data.table)

dat = readRDS("combined_results_with_SV_info.rds")
dat = dat %>% select(-`t-stat`)

df = dat[,1:38]
colnames(df)[4] = 'beta'
df$AF = as.numeric(gsub("[(),]", "", df$AF))
df$AF_PD = as.numeric(gsub("[(),]", "", df$AF_PD))
df$AF_ILBD = as.numeric(gsub("[(),]", "", df$AF_ILBD))
df$AF_HC = as.numeric(gsub("[(),]", "", df$AF_HC))

# correlate effect size (beta) with SV properties AF, SVLEN and caller support
cor_results <- df %>%
  group_by(CellType, Region, Test) %>%
  summarise(
    cor_beta_svlen   = cor(abs(beta), svlen, use = "complete.obs"),
    cor_beta_AF      = cor(abs(beta), AF, use = "complete.obs")
  ) %>%
  ungroup()

# plot res
library(tidyr)
cor_long <- cor_results %>%
  pivot_longer(
    cols = starts_with("cor_beta_"),
    names_to = "variable",
    values_to = "correlation"
  )

cor_long$variable = as.factor(cor_long$variable)
levels(cor_long$variable) = c("Allele frequency", "SV length")
write.table(cor_long, "table_corr_results_combined_SV_types.tsv", sep='\t', quote=F, row.names=F, col.names=T)

library(ggplot2)
ggplot(cor_long, aes(x = variable, y = CellType, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0, limit = c(-1, 1),
    name = "Correlation"
  ) +
  facet_grid(Region~Test) +
  labs(
    x = "Variable",
    y = "Cell Type",
    title = "Correlation between eQTL effect size and SV properties",
    subtitle = "MTG"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 14),     # Bold facet labels
    axis.text.x = element_text(angle = 45, hjust = 1, color = 'black', size=12),
    axis.text.y = element_text(color = 'black', size=12),
    axis.title.x = element_text(color='black', face='bold', size=12),
    axis.title.y = element_text(color='black', face='bold', size=12),
    plot.title = element_text(face = 'bold', size = 16),
    legend.position = "top"
  ) +
  coord_flip()

ggsave("fig5_f.pdf", height=7, width=12)

# correlation by SV type
cor_results <- df %>%
  group_by(CellType, svtype, Region, Test) %>%
  summarise(
    cor_beta_svlen   = cor(abs(beta), svlen, use = "complete.obs"),
    cor_beta_AF      = cor(abs(beta), AF, use = "complete.obs")
  ) %>%
  ungroup()

cor_long <- cor_results %>%
  pivot_longer(
    cols = starts_with("cor_beta_"),
    names_to = "variable",
    values_to = "correlation"
  )

cor_long$variable = as.factor(cor_long$variable)
levels(cor_long$variable) = c("Allele frequency", "SV length")
write.table(cor_long, "table_corr_results.tsv", sep='\t', quote=F, row.names=F, col.names=T)


library(tidyverse)
library(broom)

# Calculate correlations and p-values
cor_results <- df %>%
  group_by(CellType, svtype, Region, Test) %>%
  summarise(
    # Use tidy() to get estimate and p.value columns
    test_svlen = list(tidy(cor.test(abs(beta), svlen, use = "complete.obs"))),
    test_AF    = list(tidy(cor.test(abs(beta), AF, use = "complete.obs"))),
    .groups = "drop"
  ) %>%
  unnest(test_svlen, names_sep = "_") %>%
  unnest(test_AF, names_sep = "_") %>%
  select(
    CellType, svtype, Region, Test,
    cor_svlen = test_svlen_estimate, p_svlen = test_svlen_p.value,
    cor_AF    = test_AF_estimate,    p_AF    = test_AF_p.value
  )

# Pivot to long format for ggplot
cor_long <- cor_results %>%
  pivot_longer(
    cols = c(starts_with("cor_"), starts_with("p_")),
    names_to = c(".value", "variable"),
    names_sep = "_"
  ) %>%
  mutate(
    variable = factor(variable, levels = c("AF", "svlen"), labels = c("Allele frequency", "SV length")),
    # Create the asterisk labels based on p-value
    stars = case_when(
      p <= 0.001 ~ "****",
      p <= 0.001 ~ "***",
      p <= 0.01  ~ "**",
      p <= 0.05  ~ "*",
      TRUE      ~ ""
    )
  )

p <- ggplot(cor_long, aes(x = variable, y = CellType, fill = cor)) +
  geom_tile(color = "white") +
  # Add the asterisks here
  geom_text(aes(label = stars), color = "black", size = 5, vjust = 0.8) +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0, limit = c(-1, 1),
    name = "Correlation"
  ) +
  facet_grid(Region + Test ~ svtype) +
  labs(
    x = "Variable",
    y = "Cell Type",
    title = "Correlation between ASE/eQTL effect size and SV properties",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, color = 'black'),
    axis.text.y = element_text(color = 'black'),
    axis.title.x = element_text(color='black', face='bold', size=12),
    axis.title.y = element_text(color='black', face='bold', size=12),
    plot.title = element_text(face = 'bold', size = 16),
    legend.position = "top"
  ) +
  coord_flip()

ggsave("ext_data9_d_with_signif.pdf", plot=p, height=8, width=12)



# number of signif SVs across cell types
summary_dat <- dat %>%
  filter(is_sig == TRUE) %>%
  group_by(CellType, Region, Test) %>%
  summarise(n_variants = n_distinct(Variant), .groups = "drop")

ggplot(summary_dat, aes(x = CellType, y = n_variants, fill = Region)) +
  geom_col(position = position_dodge(preserve = "single")) + # 'dodge' puts regions side-by-side
  geom_text(aes(label = ifelse(n_variants > 0, n_variants, "")), position = position_dodge(width = 0.9),
    hjust = -0.2,
    size = 3.5,              # Adjust font size
    color = "black"
  ) +
  scale_fill_manual(values=c("MTG"='indianred',"Midbrain"='slategray4')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  facet_wrap(~Test, scales = "free_y", ncol = 2) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 14),     # Bold facet labels
    axis.text.x = element_text(color = 'black', size=12),
    axis.text.y = element_text(color = 'black', size=12),
    axis.title.x = element_text(color='black', face='bold', size=12),
    axis.title.y = element_text(color='black', face='bold', size=12),
    plot.title = element_text(face = 'bold', size = 16),
    legend.position = "top"
  ) +
  labs(
    title = "Significant SVs across cell type",
    x = "Cell Type",
    y = "Number of Unique Significant SVs",
    fill = "Region"
  ) +
  coord_flip()

ggsave("fig5_g.pdf", height=7, width=12)


