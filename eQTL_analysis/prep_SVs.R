library(vcfR)

extract_info_field <- function(info_field, field_name) {
  pattern <- paste0(field_name, "=([^;]+)")
  values <- rep(NA_character_, length(info_field))

  matches <- regexpr(pattern, info_field, perl = TRUE)
  has_match <- matches != -1

  if (any(has_match)) {
    # Extract the capture group directly using regmatches and regex capture groups
    for (i in which(has_match)) {
      match_text <- regmatches(info_field[i], regexec(pattern, info_field[i], perl = TRUE))[[1]][2]
      values[i] <- match_text
    }
  }

  return(values)
}


vcf_analysis <- function(vcf_file) {
  if (!requireNamespace("vcfR", quietly = TRUE))
    stop("Package 'vcfR' needed. Install with: install.packages('vcfR')")

  # Read VCF file once
  vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)

  # Extract genotype information
  gt <- vcfR::extract.gt(vcf, element = "GT", as.numeric = FALSE)

  # Extract variant information
  fix <- vcfR::getFIX(vcf)
  info_field <- vcf@fix[, "INFO"]

  gt_numeric <- matrix(3, nrow = nrow(gt), ncol = ncol(gt))
  rownames(gt_numeric) <- rownames(gt)
  colnames(gt_numeric) <- colnames(gt)

  gt_numeric[gt == "0/0" | gt == "0|0"] <- 0
  gt_numeric[gt == "0/1" | gt == "1/0" | gt == "0|1" | gt == "1|0"] <- 1
  gt_numeric[gt == "1/1" | gt == "1|1"] <- 2

  # Extract SVTYPE and SVLEN
  svtype <- extract_info_field(info_field, "SVTYPE")
  svlen <- abs(as.numeric(extract_info_field(info_field, "SVLEN")))

  # Create variant positions dataframe
  pos_start <- as.numeric(fix[, "POS"])
  pos_end <- pos_start
  pos_end <- pos_start + svlen

  variant_df <- data.frame(
    svid = fix[, "ID"],
    chr = fix[, "CHROM"],
    pos_start = pos_start,
    pos_end = pos_end,
    pos_mid = (pos_start + pos_end) / 2
  )

  return(list(
    genotype_matrix = gt_numeric,
    variant_positions = variant_df
  ))
}

# Use the function
vcf_file = "../PD_eGenes/final.pd.sv.vcf.gz"
res_list <- vcf_analysis(vcf_file)

write.table(res_list[[1]], "input_MatrixEQTL_genotype.tsv", sep='\t', quote=F, row.names=T, col.names=NA)
write.table(res_list[[2]], "input_MatrixEQTL_variant_positions.tsv", sep='\t', quote=F, row.names=F, col.names=T)

