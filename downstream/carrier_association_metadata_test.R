# SV carrier status vs donor/QC metadata (diagnosis, ancestry PCs, age, sex, PMI, RIN, batch) by Fisher/Wilcoxon tests on existing genotypes; writes carrier_associations.tsv.
library(data.table)

# CONFIG
ROOT <- Sys.getenv("SV_ROOT", "/stanley/levin_asap_storage")
FF <- file.path(ROOT, "kwanho/Revio_gDNA/analysis/variant_calls/finalfinal")
geno_file <- file.path(FF, "MatrixEQTL/input_MatrixEQTL_genotype.tsv")  # SVs x samples (0/1/2)
pcs_file  <- file.path(FF, "MatrixEQTL/genomic_pcs_1-6.tsv")
meta_file <- "cohort_metadata.tsv"     # sample_id + diagnosis,age,sex,PMI,RIN,batch (provide)
lead_file <- "lead_svs.txt"            # one SV_id per line
cat_vars  <- c("diagnosis","sex","batch")
num_vars  <- c("age","PMI","RIN")

geno <- as.data.frame(fread(geno_file), row.names=1)
lead <- readLines(lead_file); lead <- intersect(lead, rownames(geno))
meta <- fread(meta_file); setnames(meta, 1, "sample_id")
pcs <- fread(pcs_file); setnames(pcs, 1, "sample_id")
meta <- merge(meta, pcs, by="sample_id")
pc_vars <- grep("^PC", names(pcs), value=TRUE)

test_one <- function(carrier, var, values) {
  df <- data.frame(carrier=factor(carrier), v=values)
  df <- df[complete.cases(df),]
  if (nrow(df) < 5 || length(unique(df$carrier)) < 2) return(c(NA, NA_real_))
  if (var %in% cat_vars) {
    p <- tryCatch(fisher.test(table(df$carrier, df$v), simulate.p.value=TRUE)$p.value, error=function(e) NA_real_)
    c("fisher", p)
  } else {
    p <- tryCatch(wilcox.test(v ~ carrier, data=df)$p.value, error=function(e) NA_real_)
    c("wilcoxon", p)
  }
}

rows <- list()
for (sv in lead) {
  dos <- as.numeric(geno[sv, meta$sample_id])
  carrier <- as.integer(dos > 0)
  for (v in c(cat_vars, num_vars, pc_vars)) {
    if (!(v %in% names(meta))) next
    r <- test_one(carrier, v, meta[[v]])
    rows[[length(rows)+1]] <- data.table(SV_id=sv, variable=v, test=r[1],
                                         p=as.numeric(r[2]), n_carriers=sum(carrier, na.rm=TRUE))
  }
}
res <- rbindlist(rows)
res[, FDR := p.adjust(p, "BH")]
fwrite(res, "carrier_associations.tsv", sep="\t")
cat(sprintf("Wrote carrier_associations.tsv (%d lead SVs x %d variables)\n",
            length(lead), length(unique(res$variable))))
