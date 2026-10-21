# SV-eQTL power simulation: inject a cis-eQTL effect into negative-binomial per-cell counts, aggregate to pseudobulk, re-run the identical linear eQTL test, permute genotypes for the null, and cross-check against powerEQTL; writes power_results.tsv and fig_power.pdf (Figure S20, Table S14).
library(data.table)
library(MASS)
library(ggplot2)
library(patchwork)
set.seed(1)
outdir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1]))
if (is.na(outdir) || outdir=="") outdir <- "."

# ---- model params (representative snRNA values; calibratable from counts_with_covariates.qs) ----
n_donors <- 100
disp     <- 0.5                                   # NB dispersion (theta = 1/disp)
sigma_bio<- 0.30                                  # donor-level biological SD (log scale); representative, calibratable
mu_expr  <- c(low=0.2, med=1.0, high=5.0)         # mean counts/cell for low/med/high eGenes
maf_grid <- c(0.05, 0.10, 0.20, 0.50)
beta_grid<- c(0, 0.10, 0.25, 0.50, 0.75, 1.00)    # log fold-change per alt allele (0 = null/FPR)
# round, log-spaced cells/donor (spanning DA neurons ~5 to abundant glia ~2000).
# The classification cutoffs (0.1 mean UMI/cell; >=9 cells/donor) are drawn on Panel A as dashed
# reference lines (computed separately below), NOT added to the simulation grid / Table S14.
ncell_grid <- c(5, 15, 30, 100, 300, 1000, 2000)
reps  <- 300
alpha <- 0.05

# one simulated eQTL test. The per-donor pseudobulk is the mean of `ncell` NB counts, so its
# sampling variance is (per-cell NB variance)/ncell exactly; we draw it directly (gamma matched
# to that mean/variance) instead of looping over cells, then use the closed-form simple-regression
# p-value. This is the identical linear eQTL model, just fast.
sim_p <- function(maf, beta, mu, ncell, permute=FALSE) {
  g  <- rbinom(n_donors, 2, maf)
  u  <- rnorm(n_donors, 0, sigma_bio)        # donor biological variation (does NOT shrink with ncell)
  m  <- mu*exp(beta*g + u)                    # per-cell mean, per donor
  v  <- (m + m^2*disp)/ncell                 # technical (sampling) variance of the pseudobulk mean
  pb <- rgamma(n_donors, shape=m^2/v, rate=m/v)   # mean m, var v, positive
  y  <- log(pmax(pb, 1e-6)); gt <- if (permute) sample(g) else g
  if (length(unique(gt)) < 2 || sd(y)==0) return(NA_real_)
  r <- suppressWarnings(cor(y, gt)); if (!is.finite(r) || abs(r)>=1) return(NA_real_)
  tstat <- r*sqrt((n_donors-2)/(1-r^2)); 2*pt(-abs(tstat), n_donors-2)
}
power_of <- function(maf,beta,mu,ncell,permute=FALSE)
  mean(replicate(reps, sim_p(maf,beta,mu,ncell,permute)) < alpha, na.rm=TRUE)

grid <- CJ(maf=maf_grid, beta=beta_grid, expr=names(mu_expr), ncell=ncell_grid)
grid[, mu := mu_expr[expr]]
cat(sprintf("simulating %d grid cells x %d reps ...\n", nrow(grid), reps))
grid[, power := power_of(maf,beta,mu,ncell), by=seq_len(nrow(grid))]
# permutation-null calibration (a few representative cells)
nullchk <- CJ(maf=maf_grid, expr="med", ncell=c(5,100,2000))[, mu:=mu_expr["med"]]
nullchk[, fpr_perm := power_of(maf,0.5,mu,ncell,permute=TRUE), by=seq_len(nrow(nullchk))]

# classification-cutoff additions: expression cutoff (0.1 mean UMI/cell) and cell cutoff (9 cells/donor),
# simulated as regular grid points and appended so Table S14 includes them (other grid rows unchanged).
mu_all    <- c(low=0.2, med=1.0, high=5.0, vlow=0.1)
ncell_all <- c(5, 9, 15, 30, 100, 300, 1000, 2000)
addgrid <- CJ(maf=maf_grid, beta=beta_grid, expr=names(mu_all), ncell=ncell_all)
addgrid[, mu := mu_all[expr]]
addgrid <- addgrid[expr=="vlow" | ncell==9]        # only the new (0.1 UMI or 9-cell) combinations
addgrid[, power := power_of(maf,beta,mu,ncell), by=seq_len(nrow(addgrid))]
allgrid <- rbind(grid, addgrid)
allgrid <- allgrid[ncell != 5]   # drop 5 cells/donor: below the >=9 classification cutoff, not shown/tabulated

fwrite(allgrid, file.path(outdir,"power_results.tsv"), sep="\t")
fwrite(nullchk, file.path(outdir,"power_null_calibration.tsv"), sep="\t")
cat("empirical FPR (beta=0), should be ~", alpha, ":", round(grid[beta==0, mean(power)],3), "\n")
cat("permutation-null FPR (median):", round(median(nullchk$fpr_perm),3), "\n")

# ---- powerEQTL cross-check (donor-level, if installed) ----
pc <- NULL
if (requireNamespace("powerEQTL", quietly=TRUE)) {
  pe <- CJ(maf=maf_grid, slope=beta_grid[beta_grid>0])
  pe[, power_powerEQTL := mapply(function(m,s) tryCatch(powerEQTL::powerEQTL.SLR(
        MAF=m, slope=s, n=n_donors, sigma.y=sigma_bio, FWER=alpha, nTests=1), error=function(e) NA), maf, slope)]
  fwrite(pe, file.path(outdir,"power_powerEQTL_crosscheck.tsv"), sep="\t")
  cat("wrote powerEQTL cross-check\n")

  # ---- powerEQTL fine grid for the simulation-vs-analytic cross-check (panel C; expr-independent) ----
  fine <- CJ(maf=maf_grid, slope=seq(0.05,1,by=0.05))
  fine[, powerEQTL := mapply(function(m,s) tryCatch(powerEQTL::powerEQTL.SLR(
        MAF=m, slope=s, n=n_donors, sigma.y=sigma_bio, FWER=alpha, nTests=1), error=function(e) NA_real_), maf, slope)]
  setorder(fine, maf, slope)
  fine[, powerEQTL := nafill(powerEQTL, type="locf"), by=maf]   # analytic returns NA once power saturates ~1
  sim <- grid[expr=="med" & ncell==2000 & beta>0, .(maf, slope=beta, power_sim=power)]
  m <- merge(sim, fine, by=c("maf","slope"))
  cat(sprintf("powerEQTL vs simulation: max abs power difference = %.3f (n=%d matched points)\n",
              max(abs(m$power_sim-m$powerEQTL), na.rm=TRUE), nrow(m)))
} else cat("powerEQTL not installed - skipping analytic cross-check\n")

# ---- figure: A (power vs effect size, gene-expression x MAF grid) + B (simulation vs powerEQTL) ----
TEAL<-"#0f857a"
lvl <- c("9","30","100","300","1000","2000")   # plotted cells/donor (9 = cell cutoff & floor; 5 below cutoff, 15 omitted)
explevels <- c("5.0 UMI/cell","1.0 UMI/cell","0.2 UMI/cell","0.1 UMI/cell")  # rows top->bottom (0.1 = expr cutoff)
gA <- allgrid[ncell %in% as.integer(lvl)]
gA[, ncf := factor(as.character(ncell), levels=lvl)]
gA[, expr_lab := factor(sprintf("%.1f UMI/cell", mu), levels=explevels)]
# Panel A: expression (rows) x MAF (cols), colored by cells/donor. The 0.1 UMI/cell row and the
# 9 cells/donor line are the two classification cutoffs (noted in the figure legend).
pa <- ggplot(gA, aes(beta, power, color=ncf)) +
  geom_line() + geom_point(size=1) +
  facet_grid(expr_lab ~ paste0("MAF ", maf)) +
  geom_hline(yintercept=0.8, linetype="dotted") +
  scale_color_viridis_d(end=0.92, direction=-1) +
  labs(x=expression("effect size (log"[e]*" FC per allele)"), y="power", color="cells/donor",
       title="Power vs. effect size") + theme_bw(base_size=10)
# Panel B: simulation vs powerEQTL cross-check (medium expression, 2,000 cells/donor -- the
# donor-limited regime where the powerEQTL analytic model applies; it is expression- and
# cell-count-independent). This is a slice of the grid, not all simulated conditions.
pb <- NULL
if (exists("fine")) {
  sim_e <- grid[expr=="med" & ncell==2000 & beta>0, .(maf, slope=beta, power_sim=power)]
  pb <- ggplot() +
    geom_line(data=fine, aes(slope, powerEQTL, color="powerEQTL"), linewidth=0.9) +
    geom_point(data=sim_e, aes(slope, power_sim, color="simulated"), size=1.8) +
    facet_wrap(~paste0("MAF ", maf), nrow=1) + geom_hline(yintercept=0.8, linetype="dotted") +
    scale_color_manual(values=c("powerEQTL"="#c0392b","simulated"="#0f857a")) +
    labs(x=expression("effect size (log"[e]*" FC per allele)"), y="power", color=NULL,
         title="Simulation vs. powerEQTL power (1 UMI/cell, 2,000 cells/donor)") +
    theme_bw(base_size=10) + theme(legend.position="top")
}
combined <- if (!is.null(pb)) (pa / pb + patchwork::plot_layout(heights=c(3,1.4))) else pa
combined <- combined + plot_annotation(tag_levels="A") &
  theme(plot.tag=element_text(face="bold", size=16))
ggsave(file.path(outdir,"fig_power.pdf"), combined, width=11, height=13)
ggsave(file.path(outdir,"fig_power.png"), combined, width=11, height=13, dpi=160)
cat("wrote fig_power.{pdf,png} (A: expression x MAF grid; B: powerEQTL cross-check)\n")
