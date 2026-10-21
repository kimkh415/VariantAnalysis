# SV-eQTL lead-SNP conditioning: refit each significant SV-eQTL with the eGene's single lead eSNP dosage as a covariate (MatrixEQTL, per gene, 100 kb cis window, modelLINEAR) to classify independent regulator vs coincident with the lead-SNP haplotype; writes SNPcond_<region>.tsv.
library(MatrixEQTL)
library(dplyr)

# CONFIG
ROOT     <- "/Volumes/stanley_levin_asap_storage/kwanho/Revio_gDNA"
BCFTOOLS <- "/Users/kwanho/miniconda3/envs/svrev-cli/bin/bcftools"
FF       <- file.path(ROOT, "analysis/variant_calls/finalfinal")
STATS    <- file.path(ROOT, "analysis/variant_calls/stats")
snp_vcf_file     <- file.path(ROOT, "MEGA/comb_new.no.chr.vcf.gz")
variant_gt_file  <- file.path(FF, "MatrixEQTL/input_MatrixEQTL_genotype.tsv")
variant_pos_file <- file.path(FF, "MatrixEQTL/input_MatrixEQTL_variant_positions.tsv")
gene_pos_file    <- file.path(FF, "MatrixEQTL/eGene_positions.tsv")
lookup_file      <- "/Users/kwanho/Work/SV/analysis/SV-eQTL/lead_esnp_lookup_perGene.tsv"
rds_dir          <- "/Users/kwanho/Work/SV/analysis/SV-eQTL/rerun_inputs"   # baseline .qs -> .rds
outdir           <- "/Users/kwanho/Work/SV/analysis/SV-eQTL/results_rerun"

REGIONS <- list(
  MTG = list(
    gex_pattern = file.path(STATS, "MTG/covar_adjusted_pb_count_matrices/PD_124_gene_scaled_expr_"),
    unconditioned = file.path(rds_dir, "uncond_MTG.rds"),
    cell_types = c("Astrocytes","Endothelial_Cells","Fibroblast_Like_Cells","GABA_Neurons","GLU_Neurons",
                   "Microglia","Oligodendrocytes","OPCs","Pericytes")),
  Midbrain = list(
    gex_pattern = file.path(STATS, "midbrain/covar_adjusted_pb_count_matrices/PD_124_gene_scaled_expr_"),
    unconditioned = file.path(rds_dir, "uncond_Midbrain.rds"),
    cell_types = c("Astrocytes","DA_neurons","Endothelial_cells","Fibroblasts","GABA_neurons","Glu_neurons",
                   "Microglia","Monocytes","Oligodendrocytes","OPCs","Pericytes"))
)
cisDist <- 100000; fdr_cutoff <- 0.05
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

# ---- SNP dosage via bcftools, fetched BY POSITION (chr:bp), cached ----
snp_cache <- new.env(hash = TRUE, parent = emptyenv()); snp_samples <- NULL
parse_gt <- function(gt) {
  gt <- sub(":.*$", "", gt)
  vapply(strsplit(gt, "[/|]"), function(a)
    if (length(a) == 0 || any(a == ".")) NA_real_ else sum(as.numeric(a)), numeric(1))
}
vcf_samples <- function(vcf) {
  if (!is.null(snp_samples)) return(snp_samples)
  snp_samples <<- system(sprintf("%s query -l %s", BCFTOOLS, shQuote(vcf)), intern = TRUE); snp_samples
}
fetch_snp_pos <- function(chrom, bp, vcf) {
  key <- sprintf("%s:%d", chrom, bp)
  if (exists(key, snp_cache)) return(get(key, snp_cache))
  smp <- vcf_samples(vcf)
  region <- sprintf("%s:%d-%d", sub("^chr", "", chrom), bp, bp)
  raw <- tryCatch(system(sprintf("%s query -r %s -f '%%POS[\\t%%GT]\\n' %s",
                                 BCFTOOLS, region, shQuote(vcf)), intern = TRUE), error = function(e) character(0))
  res <- NULL
  for (line in raw) {
    p <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(p) - 1 == length(smp) && as.integer(p[1]) == bp) {
      d <- parse_gt(p[-1]); names(d) <- smp; res <- d; break
    }
  }
  assign(key, res, envir = snp_cache); res
}
# single lead-eSNP dosage row (1 x donors) for a gene; mean-imputes sporadic missing; NULL if absent/monomorphic
lead_snp_matrix <- function(g, lut, donors) {
  sr <- lut[lut$gene == g, , drop = FALSE]; if (!nrow(sr)) return(list(mat = NULL, ids = ""))
  v <- fetch_snp_pos(sr$chr[1], sr$bp[1], snp_vcf_file)
  if (is.null(v)) return(list(mat = NULL, ids = sr$snp_id[1]))       # eSNP not fetchable
  row <- v[donors]
  if (length(unique(row[!is.na(row)])) < 2) return(list(mat = NULL, ids = sr$snp_id[1]))  # monomorphic in cohort
  na <- is.na(row); if (any(na)) row[na] <- mean(row[!na])
  mat <- matrix(row, nrow = 1, dimnames = list("condSNP1", donors))
  list(mat = mat, ids = sr$snp_id[1])
}

# per-gene conditioned MatrixEQTL: ALL SVs, single gene, cvrt = lead eSNP; cisDist pairs them.
cond_one_gene <- function(g, svs_slice, gene_slice, cvd, sv_pos, gp_g) {
  res <- tryCatch(Matrix_eQTL_main(
    snps = svs_slice, gene = gene_slice, cvrt = cvd, output_file_name.cis = tempfile(),
    pvOutputThreshold.cis = 1, pvOutputThreshold = 0, useModel = modelLINEAR,
    snpspos = sv_pos, genepos = gp_g, cisDist = cisDist,
    verbose = FALSE, pvalue.hist = FALSE, min.pv.by.genesnp = FALSE, noFDRsaveMemory = FALSE),
    error = function(e) NULL)
  if (is.null(res) || nrow(res$cis$eqtls) == 0) return(NULL)
  e <- res$cis$eqtls[, c("snps","beta","pvalue")]; e$gene <- g; e
}

# ---- shared inputs ----
svdf_big <- read.table(variant_gt_file, sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)
sv_pos <- read.table(variant_pos_file, sep = "\t", header = TRUE, row.names = 1) |>
  tibble::rownames_to_column("svid"); sv_pos <- sv_pos[, c(1,2,5)]; sv_pos$chr <- gsub("chr","",sv_pos$chr)
gene_pos <- read.table(gene_pos_file, sep = "\t", header = TRUE); gene_pos$chr <- gsub("chr","",gene_pos$chr)

# lookup: one lead eSNP per eGene, use==TRUE only; map to lut cols gene/chr/bp/snp_id
lraw <- read.table(lookup_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "#")
lraw <- lraw[tolower(as.character(lraw$use)) == "true", , drop = FALSE]
lut  <- data.frame(gene = lraw$gene, chr = lraw$chr, bp = as.integer(lraw$bp),
                   snp_id = lraw$lead_eSNP, stringsAsFactors = FALSE)

vsamp <- vcf_samples(snp_vcf_file); sv_vcf_samples <- intersect(colnames(svdf_big), vsamp)
log_msg(sprintf("SVs=%d  eGenes(lookup,use)=%d  shared samples=%d", nrow(svdf_big), length(unique(lut$gene)), length(sv_vcf_samples)))

for (region in names(REGIONS)) {
  cfg <- REGIONS[[region]]; log_msg("=== region:", region, "===")
  after <- list()
  for (ct in cfg$cell_types) {
    ef <- paste0(cfg$gex_pattern, ct, ".tsv")
    if (!file.exists(ef)) { log_msg("  missing expr:", ef); next }
    expr <- read.table(ef, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
    base <- intersect(sv_vcf_samples, colnames(expr)); if (length(base) < 10) next
    svs_slice <- SlicedData$new(); svs_slice$CreateFromMatrix(as.matrix(svdf_big[, base, drop = FALSE]))
    per_gene <- list()
    for (g in intersect(unique(lut$gene), rownames(expr))) {
      sm <- lead_snp_matrix(g, lut, base); if (is.null(sm$mat)) next
      gs <- SlicedData$new(); gs$CreateFromMatrix(as.matrix(expr[g, base, drop = FALSE]))
      cvd <- SlicedData$new(); cvd$CreateFromMatrix(sm$mat)
      e <- cond_one_gene(g, svs_slice, gs, cvd, sv_pos, gene_pos[gene_pos$geneid == g, , drop = FALSE])
      if (!is.null(e)) { e$conditioning_snps <- sm$ids; per_gene[[g]] <- e }
    }
    if (length(per_gene)) {
      cb <- do.call(rbind, per_gene); cb$cell_type <- ct
      cb$FDR_combined <- p.adjust(cb$pvalue, "BH")
      after[[ct]] <- cb
    }
  }
  if (!length(after)) { log_msg("  no conditioned results"); next }
  aft <- do.call(rbind, after)
  aft$pair <- paste0(aft$snps, "_", aft$gene)

  # BEFORE = submitted cis-eQTL (RDS converted from submitted .qs) -> exact Fig.4 significance
  uncond <- dplyr::bind_rows(readRDS(cfg$unconditioned), .id = "cell_type")
  uncond$pair <- paste0(uncond$snps, "_", uncond$gene)

  cmp <- dplyr::inner_join(
    uncond |> dplyr::transmute(pair, snps, gene, cell_type, FDR_uncond = FDR, beta_uncond = beta),
    aft   |> dplyr::transmute(pair, cell_type, FDR_combined, beta_cond = beta, conditioning_snps),
    by = c("pair","cell_type"))
  cmp$region <- region
  cmp$still_sig <- !is.na(cmp$FDR_combined) & cmp$FDR_combined <= fdr_cutoff
  write.table(cmp, file.path(outdir, sprintf("SNPcond_%s_comparison.tsv", region)), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(aft, file.path(outdir, sprintf("SNPcond_%s.tsv", region)), sep = "\t", quote = FALSE, row.names = FALSE)
  orig <- cmp[cmp$FDR_uncond <= fdr_cutoff, ]
  log_msg(sprintf("  [%s] %d cell types; %d submitted-sig eQTL tests conditioned; %d survive (%.1f%%) = candidate regulators; %d tag the haplotype",
                  region, length(after), nrow(orig), sum(orig$still_sig), 100*mean(orig$still_sig), sum(!orig$still_sig)))
}
log_msg("Done.")
