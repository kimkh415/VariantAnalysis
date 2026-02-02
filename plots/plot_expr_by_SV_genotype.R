library(Seurat)
library(ggplot2)
library(qs)
library(dplyr)
library(tidyr)

# convert character to list
chr_to_list <- function(df, convert.cols) {
    require(glue)
    for (col in convert.cols) {
        print(glue("converting {col}"))
        df[[col]] = sapply(strsplit(as.character(df[[col]]),"[][']|,\\s*"), function(x) x[nzchar(x)])
    }
    return(df)
}


# get genotype
get_genotype <- function(data, chrom, pos, svtype, svlen) {
    found = data %>% filter(CHROM==chrom & POS==pos & SVTYPE==svtype & SIZE==svlen)
    if (nrow(found) != 1) print("FOUND DUPLICATED VARIANT!!!")
    gt = found$GENOTYPES[[1]]
    samps = found$SAMPLES[[1]]
    names(gt) = samps
    return(gt)
}


# get aggregate expression data
prep_expr <- function(seur, gene, cond.map, genotypes) {
    ae <- AverageExpression(seur, features=gene, group.by=c('CellType', 'sample_id'), slot='counts')[["RNA"]]
    out <- t(ae)
    colnames(out) = "expr"
    out <- as.data.frame(out)
    arr = stringr::str_split(rownames(out), '_', simplify=T)
    out$CellType = arr[,1]
    out$sample_id = arr[,2]
    rownames(out) = NULL
    out$condition = cond.map[out$sample_id]
    out$genotype = genotypes[out$sample_id]
    return(out)
}


boxplot_expr_all <- function(data, gene = "", title="", valid_genotypes, p_values) {
  valid_genotypes = intersect(valid_genotypes, unique(data$genotype))
  data$genotype <- factor(data$genotype, levels = valid_genotypes)

  p <- ggplot(data, aes(x = genotype, y = expr)) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    geom_jitter(aes(color=condition), width = 0.2, size = 1) +
    scale_color_brewer(palette="Dark2") +
    facet_wrap(~ CellType_txt, scales = "free_y", nrow=1) + 
    labs(title = glue("Gene: {gene}"),
         subtitle = glue("SV: {title}"),
         x = "Genotype",
         y = "Expression",
         color = "Condition") +
    theme_classic() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 16, face = 'bold'),
      plot.subtitle = element_text(size = 16, face = 'bold'),
      axis.title = element_text(size = 14, colour = 'black', face='bold'),
      axis.text = element_text(size = 12, colour = 'black', face='bold'),
      strip.text = element_text(size = 14, face = 'bold'),
      strip.background = element_rect(fill = "white"),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.key.size = unit(0.8, "cm"),
      panel.spacing = unit(0, "lines")
    )

  if (!is.null(p_values)) {
    get_stars <- function(p) {
      if (is.null(p) || is.na(p)) return("")
      if (p <= 0.0001) return("****")
      if (p <= 0.001)  return("***")
      if (p <= 0.01)   return("**")
      if (p <= 0.05)   return("*")
      return("")
    }

    # Calculate x_mid as the middle position of the genotypes
    num_genotypes <- length(valid_genotypes)
    x_mid <- (1 + num_genotypes) / 2

    label_df <- data %>%
      group_by(CellType_txt) %>%
      summarise(
        # Get the original CellType to look up the p-value
        orig_ct = first(CellType), 
        y_pos = max(expr, na.rm = TRUE) * 1.05, 
        .groups = 'drop'
      ) %>%
      mutate(
        p_val = p_values[orig_ct],
        label = sapply(p_val, get_stars)
      )

    p <- p + geom_text(data = label_df,
                       aes(x = x_mid, y = y_pos, label = label),
                       hjust = 0.5, vjust = 0.5,
                       inherit.aes = FALSE,
                       fontface = "bold", size = 6)
  }
  
  return(p)
}


# genotype data
#sv_file = "../../figures/merged/merged_variant_info_filtered.csv"
sv_file = "PD_eGenes_variant_info.csv"
sv = data.table::fread(sv_file, header=T)
li.cols = c("SAMPLES","CONDITION","SAMP_COND","GENE","GENOTYPES")
sv = chr_to_list(sv, li.cols)
sv = sv %>% unnest(GENE) %>% as.data.frame()

# load sample conditions
cond.map = readRDS("/stanley/levin_asap_storage/kwanho/Revio_gDNA/QC/sample_to_condition_mapping.rds")

# load snRNA
seur <- qread("seur_snRNA_MTG.qs")

# select variant
#res <- qread("results_Z/stats_result_merged.qs")
res.list <- qread("../results/MatrixEQTL_PD_eGenes_MTG_cis_eQTLs.qs")
res <- do.call(rbind, res.list)
colnames(res) = c('variant','gene',"statistic","pvalue", 'padj','effect_size', 'cell_type')
res$cell_type = factor(res$cell_type)
levels(res$cell_type) = levels(seur)

# Scatter plot: padj and eGene expression in each cell type
library(tidyverse)
expr_df <- as.data.frame(a) %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(cols = -gene, names_to = "cell_type", values_to = "expression")

plot_data <- res %>%
  left_join(expr_df, by = c("gene", "cell_type")) %>%
  mutate(
    log10_padj = -log10(padj),
    is_significant = padj < 0.05,
    # Normalize expression to range 0 to 30
    expression_norm = ((expression - min(expression, na.rm = TRUE)) / 
                       (max(expression, na.rm = TRUE) - min(expression, na.rm = TRUE))) * 30
  )
ggplot(plot_data, aes(x = expression_norm, y = log10_padj)) +
  geom_point(aes(color = is_significant), alpha = 0.7) +
  scale_color_manual(
    values = c("FALSE" = "grey", "TRUE" = "black"),
    labels = c("FALSE" = "padj >= 0.05", "TRUE" = "padj < 0.05"),
    name = "Significance"
  ) +
  facet_wrap(~ cell_type) +
  labs(
    x = "Normalized Expression",
    y = "-log10(Adjusted P-value)"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

#variant = "chr17_45904200_DEL_139_MAPT"
#variant = "chr10_15531400_INS_59_ITGA8"
#variant = "chr9_17738514_DEL_72_SH3GL2"
#variant = "chr17_46256090_INS_55_KANSL1"
#variant = "chr4_76150171_DEL_109_SCARB2"
#variant = "chr17_46237501_DEL_724_ARL17B"
#variant = "chr17_46238115_INS_230_ARL17B"
#variant = "chr10_15531400_INS_60_ITGA8"
#variant = "chr17_46237501_DEL_724_ARL17B"

res = readRDS("../../../integrate_ASE/combined_results.rds")
svs <- c("chr17_46231903_INV_45106_ARL17B", "chr8_22625880_INS_330_BIN3", "chr4_1372768_INS_323_CRIPAK", "chr4_1357015_INS_410_UVSSA",
	"chr14_55489789_DEL_-65_KTN1", "chr17_45778945_INS_55_MAPT-AS1")
#svs <- c("chr10_15531400_INS_60_ITGA8","chr4_903606_INS_59_SLC26A1","chr17_46237501_DEL_-724_ARL17B","chr17_46238115_INS_230_LRRC37A","chr12_123581929_INS_974_DDX55","chr8_22625880_INS_330_BIN3",
#         "chr12_123521586_DEL_-108_DDX55", "chr17_46076599_INS_2336_MAPT", "chr17_46231903_INV_45106_ARL17B", "chr4_932428_DEL_-84_TMEM175","chr4_15662168_INS_99_CC2D2A","chr17_46277763_DEL_-4456_ARL17B","chr17_46256090_INS_56_ARL17B","chr6_28222417_DEL_-337_ZSCAN26","chr6_28222417_DEL_-337_ZSCAN16-AS1","chr6_28176185_INS_2034_ZSCAN26","chr6_28176185_INS_2034_ZSCAN16-AS1")
#svs <- c('chr4_1372768_INS_323_CRIPAK','chr14_55489789_DEL_-65_KTN1','chr12_123521586_DEL_-108_DDX55','chr12_123581929_INS_974_DDX55')
cts = c("GABAergic neurons","Glutamatergic neurons","Astrocytes","Microglia","Oligodendrocytes","Oligodendrocyte precursors")
#cts = c('GABA neurons','Glu neurons','Astrocytes','Microglia','Oligodendrocytes','OPCs')
#VALID_GENOTYPES <- c("0|0", "0|1", "1|0", "1|1", "0/0","0/1","1/0","1/1")
VALID_GENOTYPES <- c("0|0", "0|1", "1|1")
VALID_GENOTYPES_INV <- c("./.","0/0","0/1","1/1")  # INV genotypes are not phased

snrna_samples = names(table(seur$sample_id))

for (variant in svs) {
print(variant)

arr = stringr::str_split(variant, '_', simplify=T)
chrom = arr[1]
pos = arr[2]
svtype = arr[3]
svlen = sub("^-", "", arr[4])
gene = arr[5]

genotypes = get_genotype(sv, chrom, pos, svtype, svlen)  # get variant genotypes
if (any(genotypes=="1|0")) {
  genotypes = gsub("1\\|0", "0|1", genotypes)
} else if (any(genotypes=="1/0")) {
  genotypes = gsub("1/0", "0/1", genotypes)
}
na_samps = setdiff(snrna_samples, names(genotypes))
na_gts = rep("./.", length(na_samps))
new_nams = c(names(genotypes), na_samps)
genotypes = c(genotypes, na_gts)
names(genotypes) = new_nams
dat <- prep_expr(seur, gene, cond.map, genotypes)

use_gt = NULL
if (svtype=='INV') {
  use_gt=VALID_GENOTYPES_INV
} else {
  use_gt=VALID_GENOTYPES
}
# plot
sdat <- dat %>% filter(CellType %in% cts)
sdat <- sdat[complete.cases(sdat),]
filt.res = res %>% filter(CellType %in% cts & Pair==variant & Region=='MTG' & Test=='EQTL')
p_val_vector <- setNames(filt.res$padj, filt.res$CellType)

sdat$CellType_txt = stringr::str_wrap(sdat$CellType, width=14)
p = boxplot_expr_all(sdat, gene=gene, title=glue("{chrom} {pos} {svtype} {svlen}"), valid_genotypes=use_gt, p_values=p_val_vector)
ggsave(plot=p, filename=glue("boxplot_v2_{gene}_{chrom}_{pos}_{svtype}.pdf"), height=4, width=15)

}

