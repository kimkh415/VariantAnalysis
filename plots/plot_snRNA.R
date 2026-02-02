library(Seurat)
library(qs)
library(SCpubr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

egenes = scan("covar_adjusted_pb_count_matrices/pd_124_eGenes.txt",'')

seur <- qread("seur_snRNA_MTG.qs")

plot_data <- seur@meta.data %>%
  group_by(CellType) %>%
  summarise(count = n()) %>%
  mutate(percent = (count / sum(count)) * 100)

ggplot(plot_data, aes(x = CellType, y = count, fill = CellType)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.1f%%", percent)), vjust=-0.5, fontface="bold", size=4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
	labels=label_number(scale_cut = cut_short_scale())) +
  labs(x="CellType", y="Count") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, face='bold'),
	axis.text.y = element_text(face='bold'),
	axis.title = element_text(face='bold'),
        panel.grid.major.y = element_line(color = "grey90"),
        panel.grid.minor.y = element_line(color = "grey95"),
        plot.title = element_text(face = "bold", size = 14))

ggsave("cell_type_freq_MTG.pdf", height = 5, width = 5)

#SCpubr::do_BarPlot(sample=seur, group.by='CellType', legend.position='none',plot.title='Number of cells per MTG cell type') +
#	scale_y_log10(breaks = scales::log_breaks(n = 5), labels = scales::label_number(scale_cut = cut_short_scale())) +
#	theme(panel.grid.major.y = element_line(), panel.grid.minor.y = element_line())





li = qread('covar_adjusted_pb_count_matrices/pb_counts_MTG_from_Z_select.qs')
bulk = lapply(li, function(x) x %>% tibble::rownames_to_column('gene'))
mat = bind_rows(bulk, .id='CellType')
plot_data_long <- mat %>%
  pivot_longer(
    cols = starts_with("BN"), # Selects all columns starting with BN
    names_to = "Individual",
    values_to = "Expression"
  ) %>%
  mutate(Unique_Sample_ID = paste(CellType, Individual, sep = "_"))

heatmap_matrix <- plot_data_long %>%
  select(gene, Unique_Sample_ID, Expression) %>%
  pivot_wider(names_from = Unique_Sample_ID, values_from = Expression) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# replace NA to 0
heatmap_matrix[is.na(heatmap_matrix)] <- 0

annotation_df <- plot_data_long %>%
  select(Unique_Sample_ID, CellType) %>%
  distinct() %>%
  column_to_rownames("Unique_Sample_ID")

annotation_df <- annotation_df[colnames(heatmap_matrix), , drop = FALSE]

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
ct_cols = hue_pal()(length(unique(plot_data_long$CellType)))
names(ct_cols) = unique(plot_data_long$CellType)

pdf("heatmap_eGene_expr_MTG.pdf", height=24, width=7)
Heatmap(t(scale(t(heatmap_matrix))), # Row-scaling (Z-score) is usually best for gene expression
	cluster_columns = F,
        name = "Z-score",
        col = col_fun,
        column_split = annotation_df$CellType,
	column_title_rot = 45,
        top_annotation = HeatmapAnnotation(
            CellType = annotation_df$CellType,
	    col = list(CellType=ct_cols)
        ),
        show_column_names = FALSE, # Hide individual IDs if there are too many
        show_row_names = TRUE,
        cluster_column_slices = FALSE # Keeps CellTypes in alphabetical order (or factor order)
)
dev.off()

setwd("../midbrain")
seur <- qread("seur_snRNA_midbrain.qs")
use.ct = c('Astrocytes','DA neurons','Endothelial cells','Fibroblasts','GABA neurons','Glu neurons','Microglia','Monocytes',
	'Oligodendrocytes','OPCs','Pericytes')
seur <- subset(seur, CellType %in% use.ct)
SCpubr::do_BarPlot(sample=seur, group.by='CellType', legend.position='none',plot.title='Number of cells per midbrain cell type') +
        scale_y_log10(breaks = scales::log_breaks(n = 5), labels = scales::label_number(scale_cut = cut_short_scale())) +
        theme(panel.grid.major.y = element_line(), panel.grid.minor.y = element_line())
ggsave("cell_type_freq_Midbrain.pdf", height=5, width=5)

li = qread('covar_adjusted_pb_count_matrices/pb_counts_midbrain_from_Z_select.qs')
bulk = lapply(li, function(x) x %>% tibble::rownames_to_column('gene'))
mat = bind_rows(bulk, .id='CellType')
plot_data_long <- mat %>%
  pivot_longer(
    cols = starts_with("BN"), # Selects all columns starting with BN
    names_to = "Individual",
    values_to = "Expression"
  ) %>%
  mutate(Unique_Sample_ID = paste(CellType, Individual, sep = "_"))

heatmap_matrix <- plot_data_long %>%
  select(gene, Unique_Sample_ID, Expression) %>%
  pivot_wider(names_from = Unique_Sample_ID, values_from = Expression) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# replace NA to 0
heatmap_matrix[is.na(heatmap_matrix)] <- 0

annotation_df <- plot_data_long %>%
  select(Unique_Sample_ID, CellType) %>%
  distinct() %>%
  column_to_rownames("Unique_Sample_ID")

annotation_df <- annotation_df[colnames(heatmap_matrix), , drop = FALSE]

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
ct_cols = hue_pal()(length(unique(plot_data_long$CellType)))
names(ct_cols) = unique(plot_data_long$CellType)

pdf("heatmap_eGene_expr_Midbrain.pdf", height=24, width=8)
Heatmap(t(scale(t(heatmap_matrix))), # Row-scaling (Z-score) is usually best for gene expression
	cluster_columns = F,
        name = "Z-score",
        col = col_fun,
        column_split = annotation_df$CellType,
	column_title_rot = 45,
        top_annotation = HeatmapAnnotation(
            CellType = annotation_df$CellType,
	    col = list(CellType=ct_cols)
        ),
        show_column_names = FALSE, # Hide individual IDs if there are too many
        show_row_names = TRUE,
        cluster_column_slices = FALSE # Keeps CellTypes in alphabetical order (or factor order)
)
dev.off()

