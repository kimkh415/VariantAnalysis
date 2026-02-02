library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(qs)
library(patchwork)


# function to prepare ase results
chr_to_list <- function(df, convert.cols) {
    require(glue)
    for (col in convert.cols) {
        print(glue("converting {col}"))
        df[[col]] = sapply(strsplit(as.character(df[[col]]),"[][']|,\\s*"), function(x) x[nzchar(x)])
    }
    return(df)
}

read_ase_xlsx <- function(ase.file, map.file=NA) {
	require(readxl)
	sheet_names <- excel_sheets(ase.file)
	all_sheets <- lapply(sheet_names, function(sheet) {
			read_excel(ase.file, sheet = sheet)
		})
	names(all_sheets) <- sheet_names
	out <- all_sheets %>% bind_rows(.id="CellType") %>% as.data.frame()
	out <- out[,c('CellType','Estimate','SNP','padj','Gene', 'Std..Error')]
	out$SNP = gsub("^([^_]+)_", "", out$SNP)
        out$CellType = gsub(" ", "_", out$CellType)
	colnames(out) = c('Cell type','Effect size','Variant','padj','Gene', 'SE')
	if (!is.na(map.file)) {
		vmap = read.table(map.file, sep='\t', row.names=1, skip=1)
		out$Variant = vmap[out$Variant,]
	}
	out$Pair = paste0(out$Variant,'_',out$Gene)
	out$merge_key = paste0(out$Pair, '_', out[["Cell type"]])
	return(out)
}

read_eqtl_qs <- function(eqtl.file) {
	require(qs)
	res_list <- qread(eqtl.file)
	out <- res_list %>% bind_rows(.id="CellType")
	out <- out[,c('CellType',"beta","snps","FDR","gene","statistic")]
	colnames(out) = c('Cell type','Effect size','Variant','padj','Gene', "t-stat")
	out$SE = out$`Effect size` / out$`t-stat`
	out$Pair = paste0(out$Variant,'_',out$Gene)
	out$merge_key = paste0(out$Pair, '_', out[["Cell type"]])
	return(out)
}

prep_plot<- function(df, orderby=c('Chr_num', 'Position')) {
	# order SVs
	pair_order <- df %>% arrange(across(all_of(orderby))) %>% pull(Pair) %>% unique()
	df$Pair <- factor(df$Pair, levels = rev(pair_order))
	original_order <- levels(df$Pair)
	
	# SV-gene pair name formatting
	df$Chr <- formatC(df$Chr, width=5, format='d')
	df$Position <- formatC(as.numeric(df$Position), format='fg', big.mark=",", width=12)
	df$SV_Type <- formatC(df$SV_Type, width=5, format="s")
	df$SV_Length <- formatC(df$SV_Length, width=5, format = "d")
	df$variant_formatted = sprintf("%5s\t%-12s\t%-5s\t%-6s\t%12s",
	            df$Chr, df$Position, df$SV_Type, df$SV_Length, df$Gene)
	lookup <- setNames(df$variant_formatted, df$Pair)
	formatted_levels <- lookup[original_order]
	df$Pair <- factor(df$variant_formatted, levels=formatted_levels)
	
	# compute line positions for gridlines showing different SVs
	line_positions <- c()
	cur_sv = ""
	for(i in 1:length(pair_order)) {
	  parts <- stringr::str_split(pair_order[i], "_", simplify=T)[1,1:4]
	  this_sv = paste(parts, collapse='_')
	  if (cur_sv == "") {
	    cur_sv = this_sv
	    next
	  }
	  if(cur_sv != "" && this_sv != cur_sv) {
	    line_positions <- c(line_positions, i)
	    cur_sv = this_sv
	  }
	}
	line_positions = length(pair_order) - line_positions + 2

	return(list(df, line_positions))
}

make_bubble_plot <- function(dat, pairs, plot_height, filename, title, plot_width=15) {
#plot_dat = dat %>% filter(Variant %in% unique(sub("_[^_]*$", "", pairs)))
plot_dat = dat %>% filter(Pair %in% pairs)
plot_dat <- plot_dat %>% group_by(Pair) %>% filter(any(is_sig == TRUE)) %>% ungroup()

prepped = prep_plot(plot_dat)
plot_dat = prepped[[1]]
line_positions=prepped[[2]]

pdf(filename, width=plot_width, height=plot_height)
print(ggplot(plot_dat, aes(x = `Cell type`, y = Pair)) +
  geom_point(aes(color = `Effect size`, size = neg_log10_padj), alpha = 0.9) +
  geom_point(data = subset(plot_dat, is_sig == TRUE),
             aes(size = neg_log10_padj, shape="Significant"),
             color = "black", stroke = 0.8) +
  guides(color = guide_colorbar(order = 1), size  = guide_legend(order = 2), shape = guide_legend(order = 3)) +
  geom_hline(yintercept = line_positions - 0.5, color = "black", linetype = "dotted", size = 0.5) +
  facet_grid(~ Region + Dataset, scales = "free", space = "free") +  # rows ~ columns
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + # Diverging color scale
  scale_shape_manual(name = paste0("adj. p-value < ", 0.05), values = c("Significant" = 1)) +
  theme_bw() +
  labs(
    size = "-log10(padj)",
    color = "Effect Size",
    y = "Variant - Gene Pair",
    x = "Cell Type",
    title = title
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = 'black', face = 'bold', size=12),
    axis.text.y = element_text(color = 'black', face = 'bold', hjust=0, family='mono', size=12),
    plot.title = element_text(face = 'bold', size = 16),
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(face = "bold", size = 12),
    plot.margin = margin(l = 40, r = 10, t = 10, b = 10, unit = "pt"),
    axis.text.y.left = element_text(margin = margin(r = 30, unit = "pt"))
  ))
dev.off()
}



# Midbrain
ase1 = read_ase_xlsx("ASE.midbrain.MTG.VCF.xlsx")
eqtl1 = read_eqtl_qs("MatrixEQTL_PD_eGenes_Midbrain_cis_eQTLs.qs")

# MTG
ase2 = read_ase_xlsx("Resutls.MTG.Nov14.xlsx")
eqtl2 = read_eqtl_qs("MatrixEQTL_PD_eGenes_MTG_cis_eQTLs.qs")

# combine
li = list(ase_mb=ase1, ase_mtg=ase2, eqtl_mb=eqtl1, eqtl_mtg=eqtl2)
dat = li %>% bind_rows(.id='nam')
x = data.frame(do.call('rbind', strsplit(dat$nam, '_', fixed=T)))
colnames(x) = c('Test','Region')
x= x %>% mutate(Test = toupper(Test), Region=toupper(Region))
x$Region[x$Region=='MB'] = 'Midbrain'
dat = cbind(dat, x)

# Add SV info
dat = dat %>% separate(Pair, into = c("Chr", "Position", "SV_Type", "SV_Length", "Gene"), sep = "_", remove = F)
dat$Chr_num = as.numeric(gsub('chr','',dat$Chr))
dat$Position = as.numeric(gsub(',', '', dat$Position))
dat$pos = dat$Position

# Unify cell type naming
dat = dat %>% filter(!`Cell type` %in% c('CD8_T_Cells','Glu-GABA_neurons','T_Cells','Unknown_Cluster_67'))
dat$`Cell type`[dat$`Cell type`=='Endothelial_cells'] = 'Endothelial_Cells'
dat$`Cell type`[dat$`Cell type`=='GABA_neurons'] = 'GABA_Neurons'
dat$`Cell type`[dat$`Cell type`=='GLU_Neurons'] = 'Glu_neurons'

# transform padj
dat <- dat %>% mutate(neg_log10_padj = -log10(padj), is_sig = padj < 0.05)

# convert to factors
dat$Test = as.factor(dat$Test)
dat$Region = as.factor(dat$Region)

# sig only
plot_dat <- dat %>% group_by(Pair) %>% filter(any(is_sig == TRUE)) %>% ungroup()

# filter out HLA region
dat_hla = plot_dat %>% filter(Chr_num==6 & pos > 28510120 & pos < 33480577)
hla.svs = unique(dat_hla$Variant)
plot_dat = plot_dat %>% filter(!Variant %in% hla.svs)

# Indicate those SV-gene pairs showing differential ASE in cell types
svs = scan("../conditional_eQTL/differential_ASE_across_cellTypes.tsv", '')
plot_dat$Pair[plot_dat$Pair %in% svs] = paste0(plot_dat$Pair[plot_dat$Pair %in% svs], '*')
plot_dat = plot_dat %>% separate(Pair, into = c("Chr", "Position", "SV_Type", "SV_Length", "Gene"), sep = "_", remove = F)
plot_dat$Chr_num = as.numeric(gsub('chr','',plot_dat$Chr))
plot_dat$Position = as.numeric(gsub(',', '', plot_dat$Position))
plot_dat$pos = plot_dat$Position

# prep df
prepped = prep_plot(plot_dat)
plot_dat = prepped[[1]]
line_positions=prepped[[2]]

#pdf("bubbleplot_signif_pairs_wide.pdf", width=15, height=40)
#pdf("bubbleplot_signif_pairs_no_HLA_wide.pdf", width=15, height=30)
pdf("bubbleplot_signif_pairs_no_HLA_wide_filtered_with_asterisk.pdf", width=15, height=25)
print(ggplot(plot_dat, aes(x = `Cell type`, y = Pair)) +
  geom_point(aes(color = `Effect size`, size = neg_log10_padj), alpha = 0.9) +
  geom_point(data = subset(plot_dat, is_sig == TRUE), 
             aes(size = neg_log10_padj, shape="Significant"), 
             color = "black", stroke = 0.8) +
  guides(color = guide_colorbar(order = 1), size  = guide_legend(order = 2), shape = guide_legend(order = 3)) +
  geom_hline(yintercept = line_positions - 0.5, color = "black", linetype = "dotted", size = 0.5) +
  facet_grid(~ Region + Test, scales = "free", space = "free") +  # rows ~ columns
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + # Diverging color scale
  scale_shape_manual(name = paste0("adj. p-value < ", 0.05), values = c("Significant" = 1)) +
  theme_bw() +
  labs(
    size = "-log10(padj)", 
    color = "Effect Size",
    y = "Variant - Gene Pair",
    x = "Cell Type",
    title = "ASE/eQTL Significant SV-gene Pairs in Midbrain and MTG"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = 'black', face = 'bold', size=12),
    axis.text.y = element_text(color = 'black', face = 'bold', hjust=0, family='mono', size=12),
    plot.title = element_text(face = 'bold', size = 16),
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(face = "bold", size = 12),
    plot.margin = margin(l = 40, r = 10, t = 10, b = 10, unit = "pt"),
    axis.text.y.left = element_text(margin = margin(r = 30, unit = "pt"))
  ))
dev.off()

prepped_hla = prep_plot(dat_hla)
dat_hla = prepped_hla[[1]]
line_positions=prepped_hla[[2]]

pdf("bubbleplot_signif_pairs_HLA_only_wide.pdf", width=15, height=25)
print(ggplot(dat_hla, aes(x = `Cell type`, y = Pair)) +
  geom_point(aes(color = `Effect size`, size = neg_log10_padj), alpha = 0.9) +
  geom_point(data = subset(dat_hla, is_sig == TRUE),
             aes(size = neg_log10_padj, shape="Significant"),
             color = "black", stroke = 0.8) +
  guides(color = guide_colorbar(order = 1), size  = guide_legend(order = 2), shape = guide_legend(order = 3)) +
  geom_hline(yintercept = line_positions - 0.5, color = "black", linetype = "dotted", size = 0.5) +
  facet_grid(~ Region + Test, scales = "free", space = "free") +  # rows ~ columns
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + # Diverging color scale
  scale_shape_manual(name = paste0("adj. p-value < ", 0.05), values = c("Significant" = 1)) +
  theme_bw() +
  labs(
    size = "-log10(padj)",
    color = "Effect Size",
    y = "Variant - Gene Pair",
    x = "Cell Type",
    title = "ASE/eQTL Significant SV-gene Pairs in Midbrain and MTG (HLA region)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = 'black', face = 'bold', size=12),
    axis.text.y = element_text(color = 'black', face = 'bold', hjust=0, family='mono', size=12),
    plot.title = element_text(face = 'bold', size = 16),
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(face = "bold", size = 12),
    plot.margin = margin(l = 40, r = 10, t = 10, b = 10, unit = "pt"),
    axis.text.y.left = element_text(margin = margin(r = 30, unit = "pt"))
  ))
dev.off()


# SV-gene pairs significant in both regions and those uniquely significant in one
pair_classification <- dat %>%
  group_by(Pair, Region) %>%
  summarise(is_sig_in_region = any(padj < 0.05), .groups = "drop") %>%
  pivot_wider(names_from = Region, values_from = is_sig_in_region, values_fill = FALSE) %>%
  mutate(
    Category = case_when(
      Midbrain & MTG  ~ "Common",
      Midbrain & !MTG ~ "Midbrain_Only",
      !Midbrain & MTG ~ "MTG_Only",
      TRUE            ~ "Neither"
    )
  )

write.table(pair_classification, "table_sv-gene_pair_signif_regions.tsv", sep='\t', quote=F, row.names=F, col.names=T)

common_pairs   <- pair_classification %>% filter(Category == "Common") %>% pull(Pair)
mb_only_pairs  <- pair_classification %>% filter(Category == "Midbrain_Only") %>% pull(Pair)
mtg_only_pairs <- pair_classification %>% filter(Category == "MTG_Only") %>% pull(Pair)

make_bubble_plot(dat, common_pairs, 25, filename="bubbleplot_signif_pairs_common_bet_regions.pdf", 
	title="SV-gene pairs significant in both regions")
make_bubble_plot(dat, mtg_only_pairs, 20, filename="bubbleplot_signif_pairs_only_in_MTG.pdf",
	title="SV-gene pairs significant in MTG")
make_bubble_plot(dat, mb_only_pairs, 20, filename="bubbleplot_signif_pairs_only_in_Midbrain.pdf",
	title="SV-gene pairs significant in midbrain")

#################################################################################################################
# confidence interval (cell type specificity)
dat <- readRDS("combined_results.rds")
pairs.plot = c('chr10_15531400_INS_60_ITGA8', 'chr8_22625880_INS_330_BIN3', 'chr14_55489789_DEL_-65_KTN1',
	'chr4_76150170_DEL_-109_SCARB2','chr17_45778945_INS_55_MAPT-AS1','chr12_123581929_INS_974_DDX55','chr12_123521586_DEL_-108_DDX55',
	'chr17_46237501_DEL_-724_ARL17B','chr17_46252361_INS_828_ARL17B','chr4_1324644_DEL_-94_UVSSA', 'chr4_1324644_DEL_-94_CRIPAK',
	'chr4_1372768_INS_323_UVSSA','chr4_1372768_INS_323_CRIPAK','chr4_1357015_INS_410_UVSSA', 'chr19_2349089_INS_215_LSM7',
	'chr12_123453161_INS_3316_RILPL2','chr10_119953464_DEL_-87_SEC23IP')
dir.create('ci_plots')
for (pair in pairs.plot) {
  print(pair)
  cursv = sub("(.*)_.*", "\\1", pair)
  curgene = sub(".*_", "", pair)
  sdat = dat %>% filter(Pair==pair)
  sdat$beta_ci_low = sdat$`Effect size` - sdat$SE
  sdat$beta_ci_high = sdat$`Effect size` + sdat$SE
  ggplot(sdat, aes(x=CellType, y=`Effect size`, color=is_sig)) +
  geom_point() +
  geom_errorbar(aes(ymin = beta_ci_low, ymax = beta_ci_high), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(title=paste0("SV:\t",cursv,"\nGene:\t",curgene), x = "Cell Type", y = "Effect Size") +
  facet_grid(Region ~ Test, scale='free') +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    title = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, color='black'),
    axis.text.y = element_text(size = 12, color='black')
  ) +
  coord_flip()
  ggsave(paste0('ci_plots/plot_ci_', pair,'.pdf'))
}


#################################################################################################################




plot_effect_size_scatter2 <- function(dat, file.name=NULL, plot.title="", ncol=NULL, wid=7, hei=10) {
        dat = dat[complete.cases(dat),]
        num_cell_types <- length(unique(dat$Cell.type))
        dat = dat %>% mutate(
                padj.combo = case_when(
                        padj < 0.05 & padj_select < 0.05 ~ "Both",
                        padj < 0.05 ~ "ASE",
                        padj_select < 0.05 ~ "ASE_select",
                        TRUE ~ NA_character_)) %>%
                filter(!is.na(padj.combo))
        dat$padj.combo = factor(dat$padj.combo, levels=c("Both", "ASE", "ASE_select"))
        scatter_plot <- ggplot(dat, aes(x = Effect.size, y = Effect.size_select)) +
                geom_point(aes(color = padj.combo), size = 3) +
                geom_hline(yintercept=0) +
                geom_vline(xintercept=0) +
                facet_wrap(~ Cell.type, ncol=ncol) +
                labs(title = "ASE effect size comparison",
			subtitle = "Hybrid selection vs not",
                        x = "Effect Size",
                        y = "Effect Size (selection)",
                        color = "Significant in") +
                theme_classic() +
                theme(
                        plot.title = element_text(size=16, face='bold'),
                        axis.title = element_text(size=14, face='bold'),
                        axis.text = element_text(size=12),
                        legend.title = element_text(size=14, face='bold'),
                        legend.text = element_text(size = 12),
                        strip.text = element_text(size = 14, face = "bold")  # This controls the facet labels
                )
        if (!is.null(file.name)) {
                ggsave(file.name, scatter_plot, width = wid, height = hei, dpi = 300)
                cat(paste("Plot saved to", file.name, "\n"))
        }

}

plot_effect_size_scatter2(dat, "scatter_ASE_hybrid_sel.pdf", ncol=2)


plot_padj_scatter <- function(dat, file.name=NULL, plot.title="", ncol=NULL, wid=7, hei=10) {
        dat = dat[complete.cases(dat),]
        num_cell_types <- length(unique(dat$Cell.type))
        dat = dat %>% mutate(
                padj.combo = case_when(
                        padj < 0.05 & padj_select < 0.05 ~ "Both",
                        padj < 0.05 ~ "ASE",
                        padj_select < 0.05 ~ "ASE_select",
                        TRUE ~ NA_character_)) %>%
                filter(!is.na(padj.combo))
        dat$padj.combo = factor(dat$padj.combo, levels=c("Both", "ASE", "ASE_select"))
        scatter_plot <- ggplot(dat, aes(x = -log10(padj), y = -log10(padj_select))) +
                geom_point(aes(color = padj.combo), size = 3) +
                geom_hline(yintercept=0) +
                geom_vline(xintercept=0) +
                facet_wrap(~ Cell.type, ncol=ncol) +
                labs(title = "ASE p-value comparison",
			subtitle = "Hybrid selection vs not",
                        x = "-log10(Padj)",
                        y = "-log10(Padj) (selection)",
                        color = "Significant in") +
                theme_classic() +
                theme(
                        plot.title = element_text(size=16, face='bold'),
                        axis.title = element_text(size=14, face='bold'),
                        axis.text = element_text(size=12),
                        legend.title = element_text(size=14, face='bold'),
                        legend.text = element_text(size = 12),
                        strip.text = element_text(size = 14, face = "bold")  # This controls the facet labels
                )
        if (!is.null(file.name)) {
                ggsave(file.name, scatter_plot, width = wid, height = hei, dpi = 300)
                cat(paste("Plot saved to", file.name, "\n"))
        }

}

plot_padj_scatter(dat, "scatter_ASE_padj_hybrid_sel.pdf", ncol=2)


#########################################################################################################
# Compare ASE (hybrid selection vs not)
ase1 = read_ase_xlsx("ASE.midbrain.MTG.VCF.xlsx")
ase2 = read_ase_xlsx("Resutls.MTG.Nov14.xlsx")
ase3 = read_ase_xlsx("ASE.midbrain.selection.MTG.VCF.xlsx")
ase4 = read_ase_xlsx("ASE.MTG.selection.xlsx")

li = list(ase1, ase2, ase3, ase4)
names(li) = c('Midbrain', 'MTG', 'Midbrain_select','MTG_select')

dat = bind_rows(li, .id='Dataset')
dat$Dataset = as.factor(dat$Dataset)
dat$`Cell type` = as.factor(dat$`Cell type`)
dat = dat %>% filter(!`Cell type` %in% c('CD8_T_Cells','Glu-GABA_neurons','T_Cells','Unknown_Cluster_67','Unknown_Immune_Cells'))
dat = droplevels(dat)
levels(dat$`Cell type`) = c('Astrocytes','Endothelial_cells','Endothelial_cells','Fibroblast_Like_Cells','Fibroblasts','GABA_neurons',
	'GABA_neurons','Glu_neurons','Glu_neurons','Microglia','Monocytes','Oligodendrocytes','OPCs','Pericytes')

dat$Region = as.factor(gsub("_select", "", dat$Dataset))
dat$HybridSel = as.factor(ifelse(grepl("select", dat$Dataset), T, F))

plot_data <- dat %>%
  select(Pair, Region, `Cell type`, HybridSel, `Effect size`, padj) %>%
  # Reshape to wide format to have columns for both methods
  pivot_wider(
    id_cols = c(Pair, Region, `Cell type`),
    names_from = HybridSel,
    values_from = c(`Effect size`, padj),
    names_sep = "_"
  ) %>%
  # Remove incomplete pairs (NAs)
  filter(!is.na(`Effect size_TRUE`) & !is.na(`Effect size_FALSE`)) %>%
  # Create the color category based on significance
  mutate(Significance = case_when(
    padj_TRUE < 0.05 & padj_FALSE < 0.05 ~ "Both",
    padj_TRUE < 0.05  ~ "Hybrid Selection Only",
    padj_FALSE < 0.05 ~ "Non-selected Only",
    TRUE ~ NA_character_ # Matches pairs significant in neither
  )) %>%
  # Filter: keep only if significant in at least one (remove NAs we just made)
  filter(!is.na(Significance))

regions <- unique(plot_data$Region)

# Effect size scatter plot
for (reg in regions) {  
  # Filter data for the current region
  region_data <- plot_data %>% filter(Region == reg)
  p <- ggplot(region_data, aes(x = `Effect size_TRUE`, y = `Effect size_FALSE`)) +
    geom_point(aes(color = Significance), alpha = 0.6) +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    #geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    # Facet by Cell type as requested
    facet_wrap(~`Cell type`) +
    theme_classic() +
    theme(
      plot.title = element_text(size=16, face='bold'),
      axis.title = element_text(size=14, face='bold'),
      axis.text = element_text(size=12),
      legend.title = element_text(size=14, face='bold'),
      legend.text = element_text(size = 12),
      strip.text = element_text(size = 14, face = "bold")  # This controls the facet labels
    ) +
    labs(
      title = paste("ASE Effect Size Comparison:", reg),
      x = "Effect Size (Hybrid Selection)",
      y = "Effect Size (Non-selected)",
      color = "Significance Source"
    )
  ggsave(paste0("scatter_ASE_vs_ASE-select_", reg, ".pdf"), height=10, width=10)
}

# number of unique SVs per cell type
for (is_sel in c(T, F)) {
nam = ""
if (is_sel) nam="_select"
sdat = dat %>% filter(HybridSel == is_sel)
print(table(sdat$HybridSel))

sv_counts <- sdat %>%
  filter(padj < 0.05) %>%
  mutate(SV = sub("_[^_]+$", "", Pair)) %>%
  group_by(Region, `Cell type`) %>%
  summarise(n_SVs = n_distinct(SV), .groups = "drop")

regions <- unique(sv_counts$Region)

for (reg in regions) {
  
  # Filter for current region
  plot_data <- sv_counts %>% 
    filter(Region == reg)
  
  # Create plot
  p <- ggplot(plot_data, aes(x = reorder(`Cell type`, n_SVs), y = n_SVs)) +
    geom_bar(stat = "identity", fill = "#1E90FF") + # Dodger blue
    geom_text(aes(label = n_SVs), hjust = -0.2, size = 3.5) + # Add count labels
    coord_flip() + # Horizontal bars
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 10)
    ) +
    labs(
      title = paste0("ASE", nam, ": Significant SVs per Cell Type"),
      subtitle = reg,
      x = "Cell Type",
      y = "Number of SVs"
    ) 
  ggsave(paste0("bar_n_unique_SVs_", reg, nam, ".pdf"), height=4, width=6)
}
}



sel = readRDS("sel.rds")
notsel = readRDS("notsel.rds")
dat = readRDS("initial_ase_results.rds")
sel.big = rbind(sel, dat %>% filter(Gene %in% unique(sel$Gene)))
notsel.big = rbind(notsel, dat %>% filter(Gene %in% unique(notsel$Gene)))

make_bubble_plot(sel.big, unique(sel.big$Pair), 27, 'bubble_selected_genes_only.pdf', "Genes included in the hybrid selection")
make_bubble_plot(notsel.big, unique(notsel.big$Pair), 25, 'bubble_genes_not_selected.pdf', "Genes not included in the hybrid selection")


# nSV in HLA
x = dat %>% separate(Pair, into=c('chr','pos','svtype','svlen','gene'))
x$pos = as.numeric(x$pos)
dat.hla = x %>% filter(chr=='chr6' & pos > 28510120 & pos < 33480577)

for (is_sel in c(T, F)) {
nam = ""
if (is_sel) nam="_select"
sdat = dat.hla %>% filter(HybridSel == is_sel)
print(table(sdat$HybridSel))

sv_counts <- sdat %>%
  filter(padj < 0.05) %>%
  group_by(Region) %>%
  summarise(n_SVs = n_distinct(Variant), .groups = "drop")

print(sv_counts)
}



