library(MatrixEQTL)
library(qs)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(grid)


print("Set up!")
# gene expr file pattern
gex_file_pattern = "/stanley/levin_asap_storage/kwanho/Revio_gDNA/variant_calls/stats/MTG/covar_adjusted_pb_count_matrices/PD_124_gene_scaled_expr_"
# variant genotype
variant_gt_file = "../input_MatrixEQTL_genotype.tsv"
# variant position
variant_pos_file = "../input_MatrixEQTL_variant_positions.tsv"
# gene position
gene_pos_file = "../eGene_positions.tsv"
# covariates (adding genotype PCs not considered during snRNA-seq data preparation)
covariates_file_name = "../genomic_pcs_1-6.tsv"

outdir = "results"
dir.create(outdir, showWarnings=F)

output_prefix = paste0(outdir, "/MatrixEQTL_PD_eGenes_with_6genotypePCs")

useModel = modelLINEAR # modelANOVA or modelLINEAR or modelLINEAR_CROSS
pvOutputThreshold = 1

plot_only=F

cell_types = c("Astrocytes", "Endothelial_Cells", "Fibroblast_Like_Cells", "GABA_Neurons", "GLU_Neurons", "Microglia", 
		 "Oligodendrocytes", "OPCs", "Pericytes")

# load covar
cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"
cvrt$fileOmitCharacters = "NA"
cvrt$fileSkipRows = 1
cvrt$fileSkipColumns = 1
cvrt$fileSliceSize = 2000
cvrt$LoadFile( covariates_file_name )

create_bubble_plot <- function(data, sig_threshold = 0.05, color_low = "green4", color_mid = "#FFEA00", color_high = "red3",
                             size_range = c(0, 5), max_abs_effect_size = 3, plot_title = NULL, order_by = NULL, rm.pattern = NULL,
                             cols = list(gene = "gene", cell_type = "cell_type", effect_size = "beta", padj = "FDR")) {
  data <- data[complete.cases(data), ]
  log10_padj <- -log10(data[[cols$padj]])
  size_breaks = seq(min(log10_padj), max(log10_padj), length.out = 5)
  size_labels = round(size_breaks)  #sprintf("%.2f", size_breaks)

  if (!is.null(rm.pattern)) {
    data = data %>% filter(!grepl(rm.pattern, pair))
  }

  if (!is.null(order_by)) {
    var_order = data %>% arrange(across(all_of(order_by))) %>% pull(pair) %>% unique()
    data$variant = factor(data$pair, levels=rev(var_order))
  } else {
    data$variant = data$pair
  }

  data = data %>% separate(variant, into = c("Chr", "Position", "SV_Type", "SV_Length", "Gene"), sep = "_", remove = F)
  original_order <- levels(data$variant)
  if(is.null(original_order)) {
    original_order <- unique(data$variant)
  }
  data$Chr <- formatC(data$Chr, width=5, format='d')
  data$Position <- formatC(as.numeric(data$Position), format='fg', big.mark=",", width=11)
  data$SV_Type <- formatC(data$SV_Type, width=5, format="s")
  data$SV_Length <- formatC(data$SV_Length, width=5, format = "d")
  data$variant_formatted = sprintf("%5s\t%-11s\t%-5s\t%-6s\t%12s", data$Chr, data$Position, data$SV_Type, data$SV_Length, data$Gene)
  lookup <- setNames(data$variant_formatted, data$variant)
  formatted_levels <- lookup[original_order]
  data$variant <- factor(data$variant_formatted, levels=formatted_levels)

  line_positions <- c()
  cur_sv = ""
  for(i in 1:length(var_order)) {
    parts <- str_split(var_order[i], "_", simplify=T)[1,1:4]
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
  line_positions = length(var_order) - line_positions + 2

  # Create the plot
  p <- ggplot(data) +
    geom_point(aes(
        x = !!sym(cols$cell_type),
        y = variant,
        size = -log10(!!sym(cols$padj))
      ),
      color = "gray80"
    ) +
    geom_point(data = . %>% filter(!!sym(cols$padj) <= sig_threshold),  # mark signif
      aes(
        x = !!sym(cols$cell_type),
        y = variant,
        size = -log10(!!sym(cols$padj)),
        color = pmin(pmax(!!sym(cols$effect_size), -max_abs_effect_size), max_abs_effect_size)
      )
    ) +
    geom_hline(yintercept = line_positions - 0.5, color = "black", linetype = "dotted", size = 0.5) +
    scale_color_gradient2(
      low = color_low,
      mid = color_mid,
      high = color_high,
      midpoint = 0,
      limits = c(-max_abs_effect_size, max_abs_effect_size),
      name = "Effect Size"
    ) +
    scale_size_continuous(
      name = "-log10(adj. p-value)",
      range = size_range,
      breaks = size_breaks,
      labels = size_labels
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color='black', face='bold', size=12),
      axis.text.y = element_text(color='black', face='bold', hjust=0, family='mono', size=12),
      plot.caption = element_text(hjust=0, size=8),
      plot.title = element_text(face='bold', size=16),
      plot.margin = margin(l = 40, r = 10, t = 10, b = 10, unit = "pt"),
      axis.text.y.left = element_text(margin = margin(r = 30, unit = "pt"))
    ) +
    labs(
      x = "Cell Type",
      y = "Variant - Gene Information",
      #y = "{Chromosome}_{position}_{SV type}_{SV size}_{Gene}",
      caption = paste0("significance threshold = ", sig_threshold),
      title = plot_title
    )
  return(p)
}

if (!plot_only) {
svdf_big = read.table(variant_gt_file, sep='\t', row.names=1, header=T)

cis.full = list()
trans.full = list()
for (ct in cell_types) {
print(ct)
# gene expression
gene_expr_file = paste0(gex_file_pattern, ct, ".tsv")

gene = SlicedData$new()
gene$fileDelimiter = "\t"      # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values
gene$fileSkipRows = 1          # one row of column labels
gene$fileSkipColumns = 1       # one column of row labels
gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
gene$LoadFile( gene_expr_file )

svdf = svdf_big[, intersect(colnames(svdf_big), colnames(gene))]
svs = SlicedData$new()
svs$CreateFromMatrix(as.matrix(svdf))

sv_pos = read.table(variant_pos_file, sep='\t', header=T, row.names=1)
sv_pos = sv_pos %>% tibble::rownames_to_column("svid")
sv_pos = sv_pos[, c(1,2,5)]
sv_pos$chr = gsub("chr", "", sv_pos$chr)

gene_pos = read.table(gene_pos_file, sep='\t', header=T)
gene_pos$chr = gsub("chr", "", gene_pos$chr)

# subset cvrt
cvrt_subset = cvrt$Clone()
sample_indices <- match(colnames(gene), cvrt_subset$columnNames)
cvrt_subset$ColumnSubsample(sample_indices)

output_file = paste0(output_prefix, '_', ct, '.tsv')

res <- Matrix_eQTL_main(
    snps = svs,
    gene = gene,
    cvrt = cvrt_subset,
    output_file_name.cis = output_file,
    pvOutputThreshold.cis = pvOutputThreshold,
    useModel = useModel,
    snpspos = sv_pos,
    genepos = gene_pos,
    verbose = TRUE,
    pvalue.hist = TRUE,
    cisDist = 100000,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)

cis.eqtl = res$cis$eqtls
trans.eqtl = res$trans$eqtls
cis.eqtl$cell_type = ct
trans.eqtl$cell_type = ct
cis.full[[ct]] = cis.eqtl
trans.full[[ct]] = trans.eqtl
}

qsave(cis.full, sprintf("%s_MTG_cis_eQTLs.qs", output_prefix))
qsave(trans.full, sprintf("%s_MTG_trans_eQTLs.qs", output_prefix))
} else {
cis.full = qread(sprintf("%s_MTG_cis_eQTLs.qs", output_prefix))
trans.full = qread(sprintf("%s_MTG_trans_eQTLs.qs", output_prefix))
}
df = do.call(rbind, cis.full)
df$pair = paste0(df$snps, '_', df$gene)

sig_cutoff = 0.05
sig_pairs = df %>% filter(FDR<sig_cutoff) %>% pull(pair) %>% unique()
write.table(sig_pairs, sprintf("%s_MTG_cis_significant_pairs.tsv", output_prefix), sep='\t', quote=F, row.names=F, col.names=F)
sig_vars = df %>% filter(FDR<sig_cutoff) %>% pull(snps) %>% unique()
write.table(sig_vars, sprintf("%s_MTG_cis_significant_variants.tsv", output_prefix), sep='\t', quote=F, row.names=F, col.names=F)
sres <- df %>% filter(pair %in% sig_pairs)

sres$chrom = stringr::str_split(sres$snps, '_', simplify=T)[,1]
sres$pos = stringr::str_split(sres$snps, '_', simplify=T)[,2]
sres$chrom_num = as.numeric(gsub("chr", "", sres$chrom))

#for (chr in unique(sres$chrom)) {
#print(chr)
#ssres <- sres %>% filter(chrom==chr)
#hei = length(unique(ssres$snps))/4
#hei = ifelse(hei<60, hei+5, hei)
#pdf(sprintf("%s/bubbleplot_matrixEQTL_MTG_cis_%s.pdf", outdir, chr), height=hei, width=6.2)
#print(create_bubble_plot(ssres, sig_threshold=sig_cutoff, plot_title=paste0("MTG - ",chr," cis eQTL results"), 
#			 order_by=c('chrom_num','gene')))
#dev.off()
#}

pdf(sprintf("%s/bubbleplot_Kanpig_MatrixEQTL_MTG_cis_full.pdf", outdir), height=32, width=9)
print(create_bubble_plot(sres, sig_threshold=sig_cutoff,plot_title="MTG - cis eQTL results",order_by=c('chrom_num','pos','gene')))
dev.off()

pdf(sprintf("%s/bubbleplot_Kanpig_MatrixEQTL_MTG_cis_no_HLA.pdf", outdir), height=18, width=9)
print(create_bubble_plot(sres, sig_threshold=sig_cutoff, plot_title="MTG - cis eQTL results", order_by=c('chrom_num','pos','gene'), rm.pattern="HLA-"))
dev.off()

print("DONE!")


