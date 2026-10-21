# SV carrier-status association tests (Table S9): per significant SV, (1) three Fisher contrasts of carrier status vs diagnosis (3-way HC/ILBD/PD, PD vs rest, disease vs HC) and (2) a confounding check of carrier status vs donor/QC metadata (age, sex, PMI, RIN, genotype PCs, sequencing batch; Wilcoxon/Fisher, BH-FDR). Writes carrier_by_diagnosis.tsv and carrier_associations.tsv.
library(data.table)

ROOT <- Sys.getenv("SV_ROOT", "/stanley/levin_asap_storage")
FF   <- file.path(ROOT, "kwanho/Revio_gDNA/analysis/variant_calls/finalfinal")
W    <- file.path(ROOT, "kwanho/Revio_gDNA/analysis/revision")
geno_file  <- file.path(FF, "MatrixEQTL/input_MatrixEQTL_genotype.tsv")   # SVs x samples, 0/1/2
pheno_file <- file.path(W,  "conditional_eQTL/cohort_phenotype.tsv")      # Donor_ID <tab> diagnosis (no header)
outdir     <- file.path(W,  "cohort_SVs/SV_category")
prev_file  <- file.path(outdir, "carrier_by_diagnosis.tsv")              # source of the SV set + genes/role

g0    <- fread(geno_file); svcol <- names(g0)[1]
geno  <- as.data.frame(g0[, -1, with=FALSE]); rownames(geno) <- g0[[svcol]]
pheno <- fread(pheno_file, header=FALSE, col.names=c("sample_id","dx"))
prev  <- fread(prev_file)                                                # SV_id, genes, role, ...

dx3 <- c("HC","ILBD","PD")
pheno <- pheno[dx %in% dx3]                                              # drop the non-HC/ILBD/PD donor
n <- pheno[, .N, by=dx][match(dx3, dx), N]; names(n) <- dx3             # group sizes (34/31/34)
samp <- intersect(pheno$sample_id, colnames(geno))
dxv  <- setNames(pheno$dx, pheno$sample_id)[samp]

rows <- lapply(prev$SV_id, function(sv) {
  dos <- suppressWarnings(as.numeric(geno[sv, samp])); carr <- !is.na(dos) & dos > 0
  cc  <- sapply(dx3, function(g) sum(carr & dxv==g))                     # carriers per dx
  p3  <- fisher.test(rbind(carrier=cc, non=n-cc))$p.value
  pd_c <- cc["PD"]; il_c <- cc["ILBD"]; hc_c <- cc["HC"]
  ppr <- fisher.test(rbind(PD  =c(pd_c, n["PD"]-pd_c),
                           rest=c(hc_c+il_c, (n["HC"]+n["ILBD"])-(hc_c+il_c))))$p.value
  pdz <- fisher.test(rbind(disease=c(il_c+pd_c, (n["ILBD"]+n["PD"])-(il_c+pd_c)),
                           HC     =c(hc_c,       n["HC"]-hc_c)))$p.value
  data.table(SV_id=sv, HC=cc["HC"], ILBD=cc["ILBD"], PD=cc["PD"], n_carriers=sum(cc),
             fisher_p_3way=p3, fisher_p_PDvsRest=ppr, fisher_p_diseaseVsHC=pdz)
})
res <- rbindlist(rows)
res[, FDR_3way        := p.adjust(fisher_p_3way, "BH")]
res[, FDR_PDvsRest    := p.adjust(fisher_p_PDvsRest, "BH")]
res[, FDR_diseaseVsHC := p.adjust(fisher_p_diseaseVsHC, "BH")]
out <- merge(prev[, .(SV_id, genes, role)], res, by="SV_id")
setorder(out, n_carriers, fisher_p_diseaseVsHC)
fwrite(out, file.path(outdir, "carrier_by_diagnosis.tsv"), sep="\t")

cat(sprintf("group sizes HC/ILBD/PD = %d/%d/%d | SVs=%d\n", n["HC"], n["ILBD"], n["PD"], nrow(out)))
for (col in c("fisher_p_PDvsRest","fisher_p_diseaseVsHC")) {
  fdr <- sub("fisher_p", "FDR", col)
  cat(sprintf("%-22s: min p=%.4f, nominal<0.05=%d, FDR<0.05=%d (min FDR=%.3f)\n",
      col, min(out[[col]]), sum(out[[col]]<0.05), sum(out[[fdr]]<0.05), min(out[[fdr]])))
}
rare <- out[n_carriers<=10]
cat(sprintf("rare (<=10 carriers): n=%d | PD-vs-rest nominal<0.05=%d (min %.3f) | disease-vs-HC nominal<0.05=%d (min %.3f)\n",
    nrow(rare), sum(rare$fisher_p_PDvsRest<0.05), min(rare$fisher_p_PDvsRest),
    sum(rare$fisher_p_diseaseVsHC<0.05), min(rare$fisher_p_diseaseVsHC)))
cat("RILPL2 (chr12_123453161_INS_3316): PD-vs-rest p=", round(out[SV_id=="chr12_123453161_INS_3316"]$fisher_p_PDvsRest,3),
    " disease-vs-HC p=", round(out[SV_id=="chr12_123453161_INS_3316"]$fisher_p_diseaseVsHC,3), "\n")

## ---- carrier-status confounding check vs donor/QC metadata : carrier_associations.tsv ----
# Carrier status vs age, sex, PMI, RIN, genotype PCs, and sequencing batch. Donor-level covariates
# are unioned from the per-cell-type phenofiles; genotype PCs from the MatrixEQTL PC file.
# Continuous -> Wilcoxon; categorical -> Fisher (simulated for multi-level batch, as in the former
phe_dir <- file.path(W, "saige_percelltype")
pf  <- list.files(phe_dir, pattern="^phenofile_.*\\.tsv$", full.names=TRUE)
donor_cols <- c("Donor_ID","age","sex","rin","pmi","batch")            # 'batch' = sequencing batch
cov <- unique(rbindlist(lapply(pf, function(f)
  unique(fread(f, select=intersect(donor_cols, names(fread(f, nrows=0)))))), fill=TRUE))
setnames(cov, "Donor_ID", "sample_id")
# sequencing batch: taken from the phenofiles if present, else from META_FILE (the
# cohort_metadata.tsv used by the former carrier_association_metadata_test.R).
if (!"batch" %in% names(cov)) {
  mf <- Sys.getenv("META_FILE", file.path(W, "cohort_metadata.tsv"))
  if (file.exists(mf)) {
    bt <- fread(mf); setnames(bt, 1, "sample_id")
    if ("batch" %in% names(bt)) cov <- merge(cov, unique(bt[, .(sample_id, batch)]), by="sample_id", all.x=TRUE)
  }
}
pcs0 <- fread(file.path(FF, "MatrixEQTL/genomic_pcs_1-6.tsv"))          # id + donor cols; rows PC1..PC6
pcm  <- t(as.matrix(pcs0[, -1, with=FALSE])); colnames(pcm) <- pcs0$id
pcdt <- data.table(sample_id=rownames(pcm), as.data.table(pcm))
meta <- merge(cov, pcdt, by="sample_id")
num_vars <- c("age","pmi","rin", grep("^PC", names(meta), value=TRUE))
cat_vars <- intersect(c("sex","batch"), names(meta))                   # batch tested when available
metav <- meta[match(samp, meta$sample_id)]                             # align to genotype sample order

crows <- list()
for (sv in prev$SV_id) {
  dos <- suppressWarnings(as.numeric(geno[sv, samp])); carr <- as.integer(!is.na(dos) & dos > 0)
  for (v in c(num_vars, cat_vars)) {
    val <- metav[[v]]; ok <- !is.na(val) & !is.na(carr)
    if (sum(ok) < 5 || length(unique(carr[ok])) < 2) { p <- NA_real_; tst <- NA_character_
    } else if (v %in% cat_vars) {
      sim <- length(unique(val[ok])) > 2                               # multi-level (batch) -> simulate, as before
      p <- fisher.test(table(carr[ok], val[ok]), simulate.p.value=sim, B=1e5)$p.value
      tst <- if (sim) "fisher_sim" else "fisher"
    } else { p <- wilcox.test(val[ok] ~ carr[ok])$p.value; tst <- "wilcoxon" }
    crows[[length(crows)+1]] <- data.table(SV_id=sv, variable=v, test=tst, p=p, n_carriers=sum(carr))
  }
}
dxrow <- res[, .(SV_id, variable="diagnosis", test="fisher", p=fisher_p_3way, n_carriers)]
ca <- rbind(rbindlist(crows), dxrow)
ca[, FDR := p.adjust(p, "BH"), by=variable]                            # per-variable BH across SVs
fwrite(ca, file.path(outdir, "carrier_associations.tsv"), sep="\t")
cat(sprintf("\ncovariate check: %d SVs x %d variables | any FDR<0.05: %d\n",
    length(prev$SV_id), uniqueN(ca$variable), sum(ca$FDR < 0.05, na.rm=TRUE)))
print(ca[is.finite(FDR), .(minP=round(min(p),4), minFDR=round(min(FDR),3)), by=variable][order(minFDR)])
