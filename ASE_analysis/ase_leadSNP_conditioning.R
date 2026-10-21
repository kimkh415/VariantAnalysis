# SV-ASE lead-SNP conditioning: classify each ASE-significant SV as an independent regulator vs coincident with the lead-SNP haplotype by adding the lead eSNP as a 3-level per-donor phase covariate (homozygous/cis/trans) to the beta-binomial ASE test; writes per-(pair, cell type) conditional results.
library(qs)
library(aod)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
T   <- args[1]      # scratch: combvcf/, mega_leadsnps_all.tsv, all_gene_snp.tsv, base_ase_all.tsv
ASED<- args[2]      # phasedResults dir
OUT <- args[3]

## ---- lead SNP: gene -> chr,bp ; MEGA phased GT ----
gs <- read.table(file.path(T,"all_gene_snp.tsv"), sep="\t",
                 col.names=c("chr","bp","gene","rsid"), stringsAsFactors=FALSE)
gs$key <- paste0(gs$chr,":",gs$bp)
mega <- read.table(file.path(T,"mega_leadsnps_all.tsv"), sep="\t", header=TRUE,
                   check.names=FALSE, comment.char="", stringsAsFactors=FALSE)
colnames(mega)[1:2] <- c("chr","bp"); mega$key <- paste0(mega$chr,":",mega$bp)
mega <- mega[!duplicated(mega$key),]
megaSamp <- setdiff(colnames(mega), c("chr","bp","key"))
megaGT <- as.matrix(mega[, megaSamp]); rownames(megaGT) <- mega$key

## ---- comb.vcf: SVid x sample GT ----
cf <- list.files(file.path(T,"combvcf"), pattern="\\.comb\\.vcf$", full.names=TRUE)
csamp <- sub("\\.comb\\.vcf$","",basename(cf))
combList <- setNames(lapply(cf, function(f){
  x <- read.table(f, sep="\t", stringsAsFactors=FALSE, quote="", comment.char="#", colClasses="character")
  gt <- sub(":.*","",x$V10); names(gt) <- x$V3; gt}), csamp)
allSV <- unique(unlist(lapply(combList, names)))
combGT <- matrix("./.", nrow=length(allSV), ncol=length(csamp), dimnames=list(allSV, csamp))
for (s in csamp) { g <- combList[[s]]; combGT[names(g), s] <- unname(g) }

bn <- function(x) sub(".*?(BN[0-9]+).*","\\1", x)
phase_for_pair <- function(svid, gene){
  key <- gs$key[gs$gene==gene][1]
  if (is.na(key) || !(key %in% rownames(megaGT)) || !(svid %in% rownames(combGT))) return(NULL)
  donors <- intersect(colnames(combGT), megaSamp)
  sv  <- combGT[svid, donors]; snp <- megaGT[key, donors]
  ph <- rep(NA_character_, length(donors)); names(ph) <- donors
  svhet <- sv %in% c("1/0","0/1"); homo <- snp %in% c("0|0","1|1"); snphet <- snp %in% c("1|0","0|1")
  sv_alt1 <- sv=="1/0"; snp_alt1 <- snp=="1|0"
  ph[svhet & homo] <- "homo"
  ph[svhet & snphet & (sv_alt1==snp_alt1)] <- "cis"
  ph[svhet & snphet & (sv_alt1!=snp_alt1)] <- "trans"
  data.frame(BN=donors, phase=ph, stringsAsFactors=FALSE)
}
fit_int <- function(df, form){
  df$Num <- df$alt + df$ref
  fit <- tryCatch(suppressWarnings(betabin(form, random=~1, data=df, link="logit")), error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  co <- tryCatch(summary(fit)@Coef, error=function(e) NULL); if (is.null(co)) return(NULL)
  i <- match("(Intercept)", rownames(co)); if (is.na(i)) return(NULL)
  list(est=co[i,1], p=co[i,4])
}

## ---- ASE-significant (pair, region, cellType) from published base calls ----
base <- read.table(file.path(T,"base_ase_all.tsv"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
sig <- base[!is.na(base$padj) & base$padj<0.05, ]
cat("ASE-significant cells to condition:", nrow(sig), " unique pairs:",
    nrow(unique(sig[,c("gene","svid")])), "\n")

res <- list()
for (region in c("MTG","MB")){
  s <- sig[sig$region==region,]; if (nrow(s)==0) next
  dat <- qread(file.path(ASED, if (region=="MTG") "ASE.MTG.qs" else "ASE.MB.qs"))
  dat$BN <- bn(dat$Sample)
  phcache <- list()
  for (j in seq_len(nrow(s))){
    gene <- s$gene[j]; svid <- s$svid[j]; ct <- s$cellType[j]
    pk <- paste(gene, svid)
    if (is.null(phcache[[pk]])) phcache[[pk]] <- phase_for_pair(svid, gene)
    ph <- phcache[[pk]]
    row <- data.frame(region=region, gene=gene, svid=svid, cellType=ct,
                      base_est=s$Estimate[j], base_padj=s$padj[j], base_NumSamp=s$NumSamp[j],
                      n_ASE_donors=NA, nsamp=NA, n_homo=NA, n_cis=NA, n_trans=NA, nlev=NA,
                      # signed (S10d): 3-value phase covariate, all phased donors
                      signed_ok=FALSE, signed_est=NA, signed_p=NA,
                      # homozygous-only (S10e): imbalance among lead-SNP-homozygous donors
                      homo_ok=FALSE, homo_est=NA, homo_p=NA, note="", stringsAsFactors=FALSE)
    # total donors with ASE reads for this pair x cell (before phase selection)
    d0 <- dat[dat$Gene==gene & dat$SNP==svid & dat$CellType==ct,]
    if (nrow(d0)>0){
      w0 <- d0 %>% group_by(BN) %>% summarise(tot=sum(Count), .groups="drop")
      row$n_ASE_donors <- sum(w0$tot>0)
    } else row$n_ASE_donors <- 0
    if (is.null(ph)){ row$note<-"no_leadSNP_or_SV_missing"; res[[length(res)+1]]<-row; next }
    d <- merge(d0, ph, by="BN")
    if (nrow(d)==0){ row$note<-"no_ASE_rows"; res[[length(res)+1]]<-row; next }
    w <- d %>% group_by(BN, phase) %>%
      summarise(alt=sum(Count[Allele=="alt"]), ref=sum(Count[Allele=="ref"]), .groups="drop") %>% as.data.frame()
    w <- w[!is.na(w$phase) & (w$alt+w$ref)>0,]
    row$nsamp<-nrow(w); row$n_homo<-sum(w$phase=="homo"); row$n_cis<-sum(w$phase=="cis"); row$n_trans<-sum(w$phase=="trans")
    lv <- names(which(table(w$phase)>0)); row$nlev<-length(lv)
    if (nrow(w) < 3){ row$note<-"too_few_donors"; res[[length(res)+1]]<-row; next }
    ## --- S10e: homozygous-only (test imbalance among lead-SNP-homozygous donors) ---
    if (row$n_homo >= 1){
      wh <- w[w$phase=="homo",]
      if (nrow(wh) >= 1){ fh <- fit_int(wh, cbind(alt, Num-alt) ~ 1); if(!is.null(fh)){ row$homo_ok<-TRUE; row$homo_est<-fh$est; row$homo_p<-fh$p } }
    }
    if (length(lv) < 2){ row$note<-paste0("single_phase_",lv[1]); res[[length(res)+1]]<-row; next }
    ## --- S10d: signed lead-SNP phase covariate (all phased donors): cis=+1, trans=-1, homo=0 ---
    w$phasenum <- c(homo=0, cis=1, trans=-1)[as.character(w$phase)]
    fc <- fit_int(w, cbind(alt, Num-alt) ~ phasenum)
    if (!is.null(fc)){ row$signed_ok<-TRUE; row$signed_est<-fc$est; row$signed_p<-fc$p }
    else row$note<-"signed_fit_failed"
    res[[length(res)+1]]<-row
  }
  cat(region, "done\n")
}
all <- do.call(rbind, res)
# conditional FDR per region x cellType, for each approach
all$signed_padj <- NA; all$homo_padj <- NA
for (k in unique(paste(all$region, all$cellType))){
  m <- paste(all$region, all$cellType)==k
  si <- m & all$signed_ok; if (sum(si)>0) all$signed_padj[si] <- p.adjust(all$signed_p[si], "fdr")
  hi <- m & all$homo_ok;   if (sum(hi)>0) all$homo_padj[hi] <- p.adjust(all$homo_p[hi], "fdr")
}
write.table(all, OUT, sep="\t", quote=FALSE, row.names=FALSE)
cat("wrote", nrow(all), "cell-rows ->", OUT, "\n")
