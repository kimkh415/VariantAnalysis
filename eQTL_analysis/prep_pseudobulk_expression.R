li = lapply(list.files(pattern="^covariances_adjusted"), read.table)
ct = list.files(pattern="^covariances_adjusted")
ct=gsub("covariances_adjusted.chr1.gexp.cell_type.|\\.pooled.norm.align_to_geno$", "", ct)
names(li) = ct

all_samps = scan("../../all_samples.txt", '')

# match sample names with SV genotype samples
li = lapply(li, function(x) {colnames(x)=str_split(colnames(x), '_', simplify=T)[,2]
                             colnames(x)[which(colnames(x)=="BN01160")] = "BN1160"
                             colnames(x)[which(colnames(x)=="BN1815")] = "BN1805"
			     x = x[,intersect(all_samps, colnames(x))]
			     return(x)})

saveRDS(li, "pb_counts_MTG_from_Z.rds")





###



library(qs)

# prep pb matrix focusing on 124 causal PD genes
li = readRDS("pb_counts_MTG_from_Z.rds")
x = li[-6]
names(x)[5] = "GLU_Neurons"
#gene_file = "/stanley/levin_asap_storage/kwanho/Revio_gDNA/refs/no_alt_GRCh38_ref/pd_124_genes.txt"
gene_file = "/stanley/levin_asap_storage/kwanho/Revio_gDNA/refs/no_alt_GRCh38_ref/pd_271_genes.txt"
genes_focus = scan(gene_file, '')
a = lapply(x, function(z) return(z[intersect(rownames(z), genes_focus), ]))
#qsave(a, "pb_counts_MTG_from_Z_select_124_PD_genes.qs")
qsave(a, "pb_counts_MTG_from_Z_select_271_PD_genes.qs")

for (i in 1:length(a)) {
nam = names(a)[i]
mat = a[[i]]
mat = t(scale(t(mat)))
#write.table(mat, paste0("PD_124_gene_scaled_expr_", nam, ".tsv"), sep='\t', quote=F, row.names=T, col.names=NA)
write.table(mat, paste0("PD_271_gene_scaled_expr_", nam, ".tsv"), sep='\t', quote=F, row.names=T, col.names=NA)
}

# prep gene positions
extract_gene_positions <- function(gtf_file, genes_of_interest = NULL) {
  # Read GTF file, specifying column names and handling comments
  df <- read.table(gtf_file, sep="\t", quote="", comment.char="#",
                   col.names=c("chr", "source", "feature", "start", "end",
                               "score", "strand", "frame", "attributes"))

  # Filter for gene features only
  genes <- df[df$feature == "gene", ]

  # Extract gene IDs from attributes column
  gene_names <- gsub(".*gene_name \"([^\"]+)\".*", "\\1", genes$attributes)

  # Create initial data frame with required columns
  result <- data.frame(
    geneid = gene_names,
    chr = genes$chr,
    left = genes$start,
    right = genes$end
  )

  # Filter for genes of interest if provided
  if (!is.null(genes_of_interest)) {
    result <- result[result$geneid %in% genes_of_interest, ]
  }

  return(result)
}

gtf_file = "/stanley/levin_asap_storage/kwanho/Revio_gDNA/refs/no_alt_GRCh38_ref/gencode.v46.annotation.gtf"
genepos = extract_gene_positions(gtf_file, genes_focus)
#write.table(genepos, "table_124_gene_positions.tsv", sep='\t', quote=F, row.names=F, col.names=T)
write.table(genepos, "table_271_gene_positions.tsv", sep='\t', quote=F, row.names=F, col.names=T)

