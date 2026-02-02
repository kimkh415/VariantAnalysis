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
	colnames(out) = c('CellType','Effect size','Variant','padj','Gene', 'SE')
	if (!is.na(map.file)) {
		vmap = read.table(map.file, sep='\t', row.names=1, skip=1)
		out$Variant = vmap[out$Variant,]
	}
	out$Pair = paste0(out$Variant,'_',out$Gene)
	out$merge_key = paste0(out$Pair, '_', out[["CellType"]])
	return(out)
}

read_eqtl_qs <- function(eqtl.file) {
	require(qs)
	res_list <- qread(eqtl.file)
	out <- res_list %>% bind_rows(.id="CellType")
	out <- out[,c('CellType',"beta","snps","FDR","gene","statistic")]
	colnames(out) = c('CellType','Effect size','Variant','padj','Gene', "t-stat")
	out$SE = out$`Effect size` / out$`t-stat`
	out$Pair = paste0(out$Variant,'_',out$Gene)
	out$merge_key = paste0(out$Pair, '_', out[["CellType"]])
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
print(ggplot(plot_dat, aes(x = `CellType`, y = Pair)) +
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
dat = dat %>% filter(!`CellType` %in% c('CD8_T_Cells','Glu-GABA_neurons','T_Cells','Unknown_Cluster_67'))
dat$`CellType`[dat$`CellType`=='Endothelial_cells'] = 'Endothelial_Cells'
dat$`CellType`[dat$`CellType`=='GABA_neurons'] = 'GABA_Neurons'
dat$`CellType`[dat$`CellType`=='GLU_Neurons'] = 'Glu_neurons'
dat$CellType = as.factor(dat$CellType)
levels(dat$CellType)=c("Astrocytes","Dopaminergic neurons","Endothelial cells","Fibroblast Like Cells","Fibroblasts","GABAergic neurons",
		"Glutamatergic neurons","Microglia","Monocytes","Oligodendrocytes","Oligodendrocyte precursors","Pericytes")

# transform padj
dat <- dat %>% mutate(neg_log10_padj = -log10(padj), is_sig = padj < 0.05)

# convert to factors
dat$Test = as.factor(dat$Test)
dat$Region = as.factor(dat$Region)

# save
saveRDS(dat, "combined_results.rds")


# load mask (binary; is the gene expressed?)
mtg.mask = readRDS("../MatrixEQTL/MTG/plot_expr/confident_expr_bool_mat.rds")
colnames(mtg.mask)[2] = "Endothelial cells"
mb.mask = readRDS("../MatrixEQTL/Midbrain/plot_expr/confident_expr_bool_mat_MB.rds")

# apply mask
dat_masked <- dat %>% select(-`t-stat`)
is_mb <- dat$Region == "Midbrain"
dat_masked$is_expressed_mb <- mapply(function(g, ct) {
  if(g %in% rownames(mb.mask) && ct %in% colnames(mb.mask)) {
    return(mb.mask[g, ct])
  } else {
    return(FALSE) # Assume not expressed if not in mask
  }
}, dat$Gene, as.character(dat$CellType))

is_mtg <- dat$Region == "MTG"
dat_masked$is_expressed_mtg <- mapply(function(g, ct) {
  if(g %in% rownames(mtg.mask) && ct %in% colnames(mtg.mask)) {
    return(mtg.mask[g, ct])
  } else {
    return(FALSE)
  }
}, dat$Gene, as.character(dat$CellType))

dat_masked$padj <- ifelse(
  (is_mb & !dat_masked$is_expressed_mb) | (is_mtg & !dat_masked$is_expressed_mtg),
  NA, 
  dat_masked$padj
)

dat_masked = dat_masked[complete.cases(dat_masked),]
dat_masked <- dat_masked %>% mutate(neg_log10_padj = -log10(padj), is_sig = padj < 0.05)

# Remove cell types based on number of cells filter
rm.ct.mb = c('Dopaminergic neurons','Fibroblasts','Glutamatergic neurons')
rm.idx.mb = which(dat_masked$CellType %in% rm.ct.mb & dat_masked$Region=='Midbrain')
dat_masked = dat_masked[-rm.idx.mb,]

# sig SV-gene pairs only
plot_dat <- dat_masked %>% group_by(Pair) %>% filter(any(is_sig == TRUE)) %>% ungroup()
plot_dat = plot_dat %>% separate(Pair, into = c("Chr", "Position", "SV_Type", "SV_Length", "Gene"), sep = "_", remove = F)
plot_dat$Chr_num = as.numeric(gsub('chr','',plot_dat$Chr))
plot_dat$Position = as.numeric(gsub(',', '', plot_dat$Position))
plot_dat$pos = plot_dat$Position

# filter out HLA region
dat_hla = plot_dat %>% filter(Chr_num==6 & pos > 28510120 & pos < 33480577)
hla.svs = unique(dat_hla$Variant)
plot_dat = plot_dat %>% filter(!Variant %in% hla.svs)

# processing for visual
prepped = prep_plot(plot_dat)
plot_dat = prepped[[1]]
line_positions=prepped[[2]]

# plot without HLA region
pdf("bubbleplot_signif_pairs_noHLA_with_expr_mask.pdf", width=14, height=15)
print(ggplot(plot_dat, aes(x = `CellType`, y = Pair)) +
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

# process HLA data for plotting
prepped_hla = prep_plot(dat_hla)
dat_hla = prepped_hla[[1]]
line_positions=prepped_hla[[2]]

# plot HLA region
pdf("bubbleplot_signif_pairs_HLA_only_with_expr_mask.pdf", width=12, height=15)
print(ggplot(dat_hla, aes(x = `CellType`, y = Pair)) +
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




# Count number of significant SVs
summary_dat <- dat_masked %>%
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
  facet_wrap(~Test, scales = "free_y", ncol = 1) +           # Separate panels for ASE and EQTL
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
    title = "Significant SVs per Cell Type",
    subtitle = "affecting genes confidently expressed",
    x = "Cell Type",
    y = "Number of Unique Significant SVs",
    fill = "Region"
  ) +
  coord_flip()

ggsave("barchart_nSVs_confident.pdf", height=10, width=10)



# expression-driven specificity vs potential regulatory specificity
baseline_specific_pairs <- dat %>%
  filter(is_sig == TRUE) %>%
  group_by(Region, Test, Variant, Gene) %>%
  mutate(n_sig_cts = n_distinct(CellType)) %>%
  filter(n_sig_cts == 1) %>%
  ungroup()

baseline_specific_pairs <- baseline_specific_pairs %>%
  rowwise() %>%
  mutate(
    # 1. Total count (all cell types in the mask)
    total_cts_expr = if(Region == "MTG") {
      sum(mtg.mask[Gene, ], na.rm = TRUE)
    } else {
      sum(mb.mask[Gene, ], na.rm = TRUE)
    },
    
    # 2. Reliable count (excluding low-cell-count types in Midbrain)
    reliable_cts_expr = if(Region == "MTG") {
      total_cts_expr # Or apply a similar filter for MTG if needed
    } else {
      # Only sum columns that are NOT in the remove list
      valid_cols <- setdiff(colnames(mb.mask), rm.ct.mb)
      sum(mb.mask[Gene, valid_cols], na.rm = TRUE)
    },
    
    # 3. 4-Way Classification
    specificity_type = case_when(
      total_cts_expr < 1 ~ "Low-power",
      reliable_cts_expr == 1 & total_cts_expr == 1 ~ "Expression-Driven Specificity",
      reliable_cts_expr <= 1 & total_cts_expr > 0 ~ "Low-Cell-Count Driven",
      reliable_cts_expr > 1 ~ "Regulatory Specificity",
      TRUE ~ "Other" # Fallback for edge cases
    )
  ) %>%
  ungroup()

specificity_breakdown_ct <- baseline_specific_pairs %>%
  group_by(Region, Test, CellType, specificity_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = specificity_type, values_from = count, values_fill = 0)

write.table(specificity_breakdown_ct, "table_specificity_breakdown_CT.tsv", sep='\t', quote=F, row.names=F, col.names=T)

specificity_breakdown <- baseline_specific_pairs %>%
  group_by(Region, Test, specificity_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = specificity_type, values_from = count, values_fill = 0)

write.table(specificity_breakdown, "table_specificity_breakdown.tsv", sep='\t', quote=F, row.names=F, col.names=T)

# Summarize for plotting
pdat <- baseline_specific_pairs %>%
  group_by(Region, Test, CellType, specificity_type) %>%
  summarise(count = n(), .groups = "drop")
pdat$specificity_type = factor(pdat$specificity_type)
levels(pdat$specificity_type) = c('Low Expression','Low Expression\nand Cell Number','Potential\nRegulatory specificity')

ggplot(pdat, aes(x = CellType, y = count, fill = Region)) +
  geom_col(position = position_dodge(width = 0.8, preserve = "single")) +
  
  # Add counts to ends of bars
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.8),
    hjust = -0.3,
    size = 3.5
  ) +
  
  # Facet by Test (Rows) and Specificity Type (Columns)
  facet_grid(Test ~ specificity_type, scales = "free_x") +
  
  coord_flip(clip = "off") +
  scale_fill_manual(values = c("MTG" = 'indianred', "Midbrain" = 'slategray4')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  
  # Apply your requested theme
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(color = 'black', size = 12),
    axis.text.y = element_text(color = 'black', size = 12),
    axis.title.x = element_text(color = 'black', face = 'bold', size = 12),
    axis.title.y = element_text(color = 'black', face = 'bold', size = 12),
    plot.title = element_text(face = 'bold', size = 16),
    legend.position = "top",
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    title = "Cell Type-Specific Effects Breakdown",
    x = "Cell type",
    y = "Number of significant SV-gene pairs",
    fill = "Region"
  )

ggsave("barchart_specificity_breakdown.pdf", width=11, height=7)

# separate region
for (reg in unique(pdat$Region)) {
pdat1 = pdat %>% filter(Region==reg)
ggplot(pdat1, aes(x = CellType, y = count, fill = Test)) +
  geom_col(position = position_dodge(width = 0.8, preserve = "single")) +
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.8),
    hjust = -0.3,
    size = 3.5
  ) +
  # one panel per Region × specificity_type
  facet_grid(Region ~ specificity_type, scales = "free_y") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c("ASE" = "indianred", "EQTL" = "slategray4")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title.x = element_text(color = "black", face = "bold", size = 12),
    axis.title.y = element_text(color = "black", face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "top",
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    title = paste0(reg, ": Potential Cell Type-Specific Effects"),
    x = "Cell type",
    y = "Number of significant SV-gene pairs",
    fill = "Test"
  )

ggsave(paste0("barchart_specificity_breakdown_", reg, "_v2.pdf"), width=10, height=5)
}

