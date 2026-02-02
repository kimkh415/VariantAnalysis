library(GenomicRanges)
library(regioneR)
library(GenomeInfoDb)
library(ggplot2)

# Permutation test
# Simulate SV position overlapping with regulatory elements by chance


# Read SV pos
# columns = chr, start, end
df <- readRDS('sv_pos_only.rds')
mask <- readRDS('mask.rds')  # gaps, centromeres, and segemental duplications

sv.gr <- makeGRangesFromDataFrame(df,
	seqnames.field = "chr",
	start.field = "start",
	end.field = "end")

mask.gr <- makeGRangesFromDataFrame(mask, 
	seqnames.field="chr", 
	start.field="start", 
	end.field="end")

# Set boundaries
hg38_lengths <- getChromInfoFromUCSC("hg38")
hg38_lengths <- hg38_lengths[hg38_lengths$chrom %in% seqlevels(sv.gr), ]
seqlengths(sv.gr) <- setNames(hg38_lengths$size, hg38_lengths$chrom)
seqlengths(mask.gr) <- setNames(hg38_lengths$size, hg38_lengths$chrom)
mask.gr = trim(mask.gr)

genome(sv.gr) <- "hg38"
genome(mask.gr) <- "hg38"

# Set universe to sample new positions from
universe <- filterChromosomes(getGenome("hg38"), keep.chr = seqlevels(sv.gr))
allowed_universe <- setdiff(universe, mask.gr)

# Function to resample SV position from allowed universe
myResample <- function(A, universe, ...) {
  w <- width(A)  # SV size
  n <- length(A) # num SVs
  # pick random location from universe
  # sample with replacement
  # sampling probability weighted by the size of each block in universe
  selected_blocks <- universe[sample(seq_along(universe), size = n, replace=T, prob=width(universe))]
  # pick a random starting pos
  max_starts <- width(selected_blocks) - w
  max_starts[max_starts < 0] <- 0 # in case SV is wider than the block
  pos <- floor(runif(n, min = 0, max = max_starts))
  # randomized SV pos
  new_starts <- start(selected_blocks) + pos
  new_svs <- GRanges(seqnames = seqnames(selected_blocks),
                     ranges = IRanges(start = new_starts, width = w))
  seqinfo(new_svs) <- seqinfo(A)
  return(new_svs)
}

# Regulatory locations
#nams = c('all_RE','brain_active_RE','cpg','DNase','H1hescH3k27ac','H1hescH3k4me1','tfbs', 'phastCons100way')
nams = c('repeatMasker')

for (nam in nams) {
print(nam)
reg <- readRDS(paste0('region_',nam,'.rds'))
reg.gr <- makeGRangesFromDataFrame(reg,
        seqnames.field = "chr",
        start.field = "start",
        end.field = "end")

seqlengths(reg.gr) <- setNames(hg38_lengths$size, hg38_lengths$chrom)
genome(reg.gr) <- "hg38"

res.file = paste0("res_permTest_", nam, ".rds")
if (!file.exists(res.file)) {
# Test
res <- permTest(
	A = sv.gr, 
	B = reg.gr,
	mask = mask.gr,
	ntimes = 1000, 
	alternative = 'auto',
	randomize.function = myResample,
	randomize.function.name = "_SV_interval",
	universe = allowed_universe,
	evaluate.function = numOverlaps,
	count.once = T)

saveRDS(res, res.file)
} else {
res <- readRDS(res.file)
}

obs_hits <- res$numOverlaps$observed
exp_hits <- mean(res$numOverlaps$permuted)
total_n  <- length(sv.gr)

obs_pct <- (obs_hits / total_n) * 100
exp_pct <- (exp_hits / total_n) * 100
fold_enrichment <- obs_hits / exp_hits

cat(paste0("Percentage of SVs overlapping ", nam, " (Observed): ", round(obs_pct, 2), "%\n",
           "Percentage of SVs overlapping ", nam, " (Expected): ", round(exp_pct, 2), "%\n",
           "Fold Enrichment: ", round(fold_enrichment, 2), "x\n",
           "Z-score: ", round(res$numOverlaps$zscore, 2), "\n",
	   "P-value: ", round(res$numOverlaps$pval, 6), "\n"))

write.table(t(c(obs_pct, exp_pct, fold_enrichment, res$numOverlaps$zscore)), paste0("table_", nam, ".tsv"), sep='\t', quote=F, 
	row.names=F, col.names=F)

pdf(paste0("plot_perm_pval_", nam, ".pdf"), height=5, width=7)
par(mar = c(5, 6, 4, 2) + 0.1)
print(plot(res))
print(title(main=nam, adj=0))
dev.off()
}

