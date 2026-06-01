# ============================================================
# LFMM Pipeline — Full Genome-Environment Association Workflow
# Author : Lidong
# Version: 1.0
#
# Description:
#   End-to-end pipeline for LFMM-based genotype-environment
#   association analysis. Four sequential steps:
#
#   Step 0 — Split merged PLINK bed into per-chromosome files
#             (calls PLINK via system())
#   Step 1 — Install dependencies; convert PED to LFMM format
#   Step 2 — Fit LFMM models; compute GC-corrected p-values
#   Step 3 — Extract per-env z-scores; merge genome-wide tables
#
# Usage:
#   Rscript lfmm_pipeline.R >> output.log 2>&1
#
# Requirements:
#   - PLINK 1.9 available on $PATH
#   - R packages: BiocManager, LEA, CMplot, tidyverse,
#                 data.table, foreach, doParallel
#
# Input files (must be in working directory before running):
#   - chr919qcldbeagle.bed / .bim / .fam
#   - env.env   (rows = individuals, cols = env variables)
#
# Output files:
#   - chr{N}.bed/.bim/.fam/.ped/.map/.lfmm/.geno  (Steps 0-1)
#   - results2/chr{N}/chr{N}_pvalues.txt           (Step 2)
#   - chr{N}_zs.median.d{D}.csv                    (Step 3A)
#   - d{D}zs.median.csv                            (Step 3B)
# ============================================================

cat("========================================================\n")
cat("LFMM pipeline started:", format(Sys.time()), "\n")
cat("========================================================\n")


# ============================================================
# Step 0: Split merged PLINK file into per-chromosome subsets
# ============================================================
# Calls PLINK via system(). PLINK 1.9 must be on $PATH.
#
# --make-bed   : write binary .bed/.bim/.fam per chromosome
# --recode     : also write .ped/.map, required by ped2geno()
#                in Step 1. Remove this flag if .ped files
#                already exist or disk space is limited.
#
# NOTE: PLINK expects numeric chromosome codes (1–18) by default.
#       If your BIM file uses "chr1"-style prefixes, append
#       --allow-extra-chr to the plink command below.

cat("\n[Step 0] Splitting merged PLINK file by chromosome...\n")

for (i in 1:18) {
  cmd <- paste(
    "plink",
    "--bfile chr919qcldbeagle",
    "--chr",   i,
    "--make-bed",
    "--recode",
    "--out",   paste0("chr", i)
  )
  exit_code <- system(cmd)
  if (exit_code != 0) {
    stop(paste("PLINK failed for chromosome", i,
               "— check that plink is on $PATH and the input files exist."))
  }
}

cat("[Step 0] Done.\n")


# ============================================================
# Step 1: Install dependencies and convert PED to LFMM format
# ============================================================
# ped2geno() reads chr{N}.ped + chr{N}.map and writes chr{N}.geno.
# geno2lfmm() reads chr{N}.geno and writes chr{N}.lfmm.
#
# Both functions write output files to the same directory as
# the input. Set setwd() at the top of the script if needed.
#
# The install blocks below are skipped if packages are already
# present; safe to leave in for reproducibility on new systems.

cat("\n[Step 1] Installing dependencies and converting to LFMM format...\n")

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!require("LEA", quietly = TRUE)) {
  BiocManager::install("LEA")
}
if (!require("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
if (!require("CMplot", quietly = TRUE)) {
  install.packages("CMplot")
}
if (!require("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
if (!require("foreach", quietly = TRUE)) {
  install.packages("foreach")
}
if (!require("doParallel", quietly = TRUE)) {
  install.packages("doParallel")
}

library(LEA)
library(tidyverse)
library(CMplot)
library(data.table)
library(foreach)
library(doParallel)

for (chr in 1:18) {
  ped_file  <- paste0("chr", chr, ".ped")
  geno_file <- paste0("chr", chr, ".geno")

  if (!file.exists(ped_file)) {
    warning(paste(
      "PED file not found for chromosome", chr,
      "— skipping. Ensure --recode was included in Step 0."
    ))
    next
  }

  ped2geno(ped_file)
  geno2lfmm(geno_file)

  message("Chromosome ", chr, " converted to LFMM format.")
}

cat("[Step 1] Done.\n")


# ============================================================
# Step 2: Fit LFMM models and compute GC-corrected p-values
# ============================================================
# One LFMM model is fit per chromosome using all env variables
# in env.env simultaneously (multivariate mode).
#
# K = 9         : number of latent factors; determine this
#                 empirically using sNMF or the Tracy-Widom
#                 test before running this pipeline.
# repetitions   : number of independent MCMC runs; the median
#                 z-score across runs is used for inference.
# CPU           : set to the number of physical cores available.
#
# project = "new" overwrites any existing project with the same
# name in the working directory without prompting. Archive or
# rename prior results before re-running.
#
# Output: results2/chr{N}/chr{N}_pvalues.txt

cat("\n[Step 2] Fitting LFMM models and computing p-values...\n")

for (chr in 1:18) {
  tryCatch({

    lfmm_file  <- paste0("chr", chr, ".lfmm")
    bim_file   <- paste0("chr", chr, ".bim")
    output_dir <- paste0("results2/chr", chr, "/")
    pval_file  <- paste0(output_dir, "chr", chr, "_pvalues.txt")

    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    if (!file.exists(lfmm_file) || !file.exists(bim_file)) {
      stop(paste("Input file missing for chromosome", chr))
    }

    obj.lfmm <- lfmm(
      input.file       = lfmm_file,
      environment.file = "env.env",
      K                = 9,
      project          = "new",
      CPU              = 50,
      repetitions      = 5
    )

    # K in lfmm.test() must equal K used in lfmm() above.
    # A mismatch does not raise an error but silently biases results.
    lfmm.result <- lfmm.test(obj.lfmm, K = 9)
    zs          <- lfmm.result$z.scores

    zs.median <- apply(zs, MARGIN = 1, median)

    # Lambda is estimated from the full z-score matrix (all SNPs × reps).
    # The divisor 0.456 is the median of a chi-squared(1) distribution,
    # which is standard practice for genomic control correction.
    lambda <- median(zs^2) / 0.456
    print(paste("Chromosome", chr, "- Lambda:", round(lambda, 4)))

    adj.p.values <- pchisq(zs.median^2 / lambda, df = 1, lower.tail = FALSE)

    bim    <- read.table(bim_file, header = FALSE)
    snp.id <- bim$V2

    pvalue.df <- data.frame(SNP = snp.id, P = adj.p.values)
    write.table(pvalue.df,
                file      = pval_file,
                sep       = "\t",
                row.names = FALSE,
                col.names = TRUE,
                quote     = FALSE)

    print(paste("Chromosome", chr, "complete. Results saved to", output_dir))

  }, error = function(e) {
    print(paste("Error on chromosome", chr, ":", e$message))
  })
}

cat("[Step 2] Done.\n")


# ============================================================
# Step 3A: Extract per-env z-scores for all chromosomes
# ============================================================
# Loads the saved .lfmmProject from Step 2 for each chromosome,
# then extracts z-scores for each of the 18 environmental
# variables and writes one CSV per chromosome × env combination.
#
# Parallelized over chromosomes. Each worker processes all 18
# env variables for its assigned chromosome.
#
# Raw_P is NOT GC-corrected here. Apply lambda correction before
# using these p-values in downstream candidate SNP detection.
#
# Output: chr{N}_zs.median.d{D}.csv

cat("\n[Step 3A] Extracting per-env z-scores (parallelized)...\n")

chromosomes <- 1:18
env_vars    <- 1:18

# Adjust worker count to match available physical cores.
# More workers than cores does not improve performance.
cl <- makeCluster(6)
registerDoParallel(cl)

foreach(chr = chromosomes, .packages = c("LEA", "data.table")) %dopar% {

  # The .lfmmProject file is created by Step 2.
  # Errors inside parallel workers are silently dropped —
  # if output CSVs are absent after this step, verify that
  # the project file exists and is not zero-byte.
  project_file <- paste0("chr", chr, "_env.lfmmProject")
  project      <- load.lfmmProject(project_file)

  bim_file <- paste0("chr", chr, ".bim")
  bim <- fread(bim_file,
               select    = c(1:4),
               col.names = c("CHR", "SNPID", "POS_CM", "POS_BP"))

  for (d in env_vars) {
    zs        <- z.scores(project, K = 9, d = d)
    zs_median <- apply(zs, 1, median)

    chr_d_data <- data.table(
      CHR      = bim$CHR,
      SNPID    = bim$SNPID,
      POS_BP   = bim$POS_BP,
      Z_median = zs_median,
      Raw_P    = pchisq(zs_median^2, df = 1, lower.tail = FALSE)
    )

    fwrite(chr_d_data, paste0("chr", chr, "_zs.median.d", d, ".csv"))
  }

  return(paste("Chromosome", chr, "processed"))
}

stopCluster(cl)
cat("[Step 3A] Done.\n")


# ============================================================
# Step 3B: Merge per-chromosome files into genome-wide tables
# ============================================================
# For each of the 18 env variables, collect all per-chromosome
# CSVs produced by Step 3A, sort by chromosome and position,
# and write one genome-wide file.
#
# list.files() searches the working directory only. If Step 3A
# output was written elsewhere, update the path= argument below.
#
# Output: d{D}zs.median.csv  (one per env variable)

cat("\n[Step 3B] Merging per-chromosome files genome-wide...\n")

for (d in 1:18) {

  chr_files <- list.files(pattern = paste0("_zs\\.median\\.d", d, "\\.csv$"))

  if (length(chr_files) == 0) {
    warning(paste("No files found for env variable d =", d, "— skipping."))
    next
  }

  merged_data <- rbindlist(lapply(chr_files, fread))
  setorder(merged_data, CHR, POS_BP)

  fwrite(merged_data, paste0("d", d, "zs.median.csv"))

  message("Env variable d=", d, " merged — ", nrow(merged_data), " SNPs total.")
}

cat("[Step 3B] Done.\n")

cat("\n========================================================\n")
cat("Pipeline complete:", format(Sys.time()), "\n")
cat("========================================================\n")
