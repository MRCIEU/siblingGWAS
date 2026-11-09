#!/usr/bin/env Rscript
# Simulates inputs for within-family.R and siblinggwas (Rust), including edge-case variants.
#
# Usage: Rscript simulate.R <out dir> <n unrelated> <n sibling families> <n random variants> <n PCs> [scenario]
#   scenario "standard" (default): numeric, logical and factor covariates with missing values
#   scenario "degenerate": additionally a collinear covariate and a factor with one level among siblings
#   scenario "clean": as "standard" but without the edge-case variants (whose NA results
#                     within-family-plots.R cannot handle)

suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
n_unrel <- as.integer(args[2])
n_fam <- as.integer(args[3])
n_rand <- as.integer(args[4])
n_pc <- as.integer(args[5])
scenario <- if (length(args) >= 6) args[6] else "standard"
set.seed(20261014)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Individuals: unrelated (FID == IID) followed by sibling families of 2-3.
fam_size <- sample(2:3, n_fam, replace = TRUE, prob = c(0.8, 0.2))
fid_fs <- rep(seq_len(n_fam), fam_size)
n_sib <- length(fid_fs)
n <- n_unrel + n_sib
iid <- 1000000L + seq_len(n)
is_sib <- c(rep(FALSE, n_unrel), rep(TRUE, n_sib))
fam <- c(seq_len(n_unrel) + n_fam, fid_fs) # family index for simulation

# Genotypes: siblings inherit from shared parental allele frequencies.
sim_snp <- function(maf) {
  g <- rbinom(n, 2, maf)
  par <- matrix(rbinom(2 * n_fam * 2, 1, maf), ncol = 4) # two parents x two alleles per family
  g[is_sib] <- par[cbind(fid_fs, sample(1:2, n_sib, TRUE))] + par[cbind(fid_fs, sample(3:4, n_sib, TRUE))]
  g
}
G <- sapply(runif(n_rand, 0.05, 0.5), sim_snp)
colnames(G) <- paste0("rs", seq_len(n_rand))
G[sample(length(G), round(0.01 * length(G)))] <- NA

# Edge cases.
edge <- list()
edge$mono <- rep(0, n)                                   # aliased in both models
edge$sib_mono <- ifelse(is_sib, 1, rbinom(n, 2, 0.3))     # FAM_MEAN and CENTRED aliased in model 2
nowithin <- rbinom(n, 2, 0.4)
nowithin[is_sib] <- rbinom(n_fam, 2, 0.4)[fid_fs]          # identical within family: WF aliased
edge$no_within <- nowithin
edge$all_na <- rep(NA_real_, n)                           # model 1 has no observations
sibs_na <- rbinom(n, 2, 0.2); sibs_na[is_sib] <- NA       # model 2 has no observations
edge$sibs_na <- sibs_na
edge$dosage <- round(pmin(2, pmax(0, sim_snp(0.3) + rnorm(n, 0, 0.15))), 3)
heavy <- sim_snp(0.25); heavy[sample(n, round(0.3 * n))] <- NA
edge$heavy_na <- heavy
E <- do.call(cbind, edge)
geno <- if (scenario == "clean") G else cbind(G, E)
nv <- ncol(geno)

# Covariates.
pcs <- matrix(rnorm(n * n_pc, sd = 0.02), ncol = n_pc, dimnames = list(NULL, paste0("PC", seq_len(n_pc))))
age <- round(runif(n, 40, 70), 1)
sex <- rbinom(n, 1, 0.5)
batch <- sample(c("Axiom_b1", "Axiom_b2", "BiLEVE_b1"), n, replace = TRUE)
batch[sample(which(!is_sib), 50)] <- "Axiom_b9_unrel" # a level with no siblings
smoker <- sample(c("TRUE", "FALSE"), n, replace = TRUE)
cov <- data.table(FID = ifelse(is_sib, iid, iid), IID = iid, age = age, sex = sex, pcs, batch = batch, smoker = smoker)
if (scenario == "degenerate") {
  cov$age_months <- cov$age * 12                       # collinear with age: aliased
  cov$centre <- ifelse(is_sib, "Bristol", sample(c("Bristol", "Leeds"), n, replace = TRUE)) # 1 level among siblings
}
cov$age[sample(n, 20)] <- NA

# Phenotype: direct effect of rs1 and rs2, family effect, covariates.
famfx <- rnorm(max(fam), 0, 0.7)[fam]
g1 <- ifelse(is.na(G[, 1]), 0, G[, 1]); g2 <- ifelse(is.na(G[, 2]), 0, G[, 2])
y <- 0.3 * g1 + 0.15 * g2 + 0.02 * age + 0.4 * sex + 30 * pcs[, 1] + famfx + rnorm(n)
y[sample(n, 15)] <- NA
phen <- data.table(FID = iid, IID = iid, bmi = round(y, 5))

# Plink .raw (with A1 suffix) and .bim
alleles <- matrix(sample(c("A", "C", "G", "T"), 2 * nv, replace = TRUE), ncol = 2)
raw <- data.table(FID = iid, IID = iid, PAT = 0, MAT = 0, SEX = sex + 1, PHENOTYPE = -9)
raw <- cbind(raw, setnames(as.data.table(geno), paste0(colnames(geno), "_", alleles[, 1])))
bim <- data.table(V1 = sample(1:22, nv, replace = TRUE), V2 = colnames(geno), V3 = 0,
                  V4 = sample(1e5:1e8, nv), V5 = alleles[, 1], V6 = alleles[, 2])
sibs <- data.table(IID = iid[is_sib], FID_FS = fid_fs)

# Shuffle row order between files so that the merge matters.
fwrite(raw[sample(.N)], file.path(outdir, "geno.raw"), sep = "\t", na = "NA")
fwrite(bim, file.path(outdir, "geno.bim"), sep = "\t", col.names = FALSE)
fwrite(phen[sample(.N)], file.path(outdir, "pheno.txt"), sep = "\t", na = "NA")
fwrite(cov[sample(.N)], file.path(outdir, "covariates.txt"), sep = "\t", na = "NA")
fwrite(sibs[sample(.N)], file.path(outdir, "siblings.txt"), sep = "\t")
cat(sprintf("Wrote %d individuals (%d siblings in %d families), %d variants, %d covariate columns to %s\n",
            n, n_sib, n_fam, nv, ncol(cov) - 2, outdir))
