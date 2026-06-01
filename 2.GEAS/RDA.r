# ============================================================
# Redundancy Analysis (RDA) for Genotype-Environment Association
# Author: Lidong
# Version: 1.0
# Description:
#   Identifies candidate SNPs associated with environmental
#   variables using RDA-based genome-environment association.
#
# Input files (place in working directory):
#   - snp.csv  : genotype matrix (rows = individuals, cols = SNPs)
#                encoded as 0/1/2 (ref homozygous/het/alt homozygous)
#   - env.csv  : environmental variables table
#                (first col = individual ID, remaining cols = predictors)
#
# Output files:
#   - rda_results.csv           : adjusted R2 and eigenvalue summary
#   - vif_values.csv            : variance inflation factors per predictor
#   - load_rda_species_scores.csv : SNP loadings on RDA axes 1-3
#   - cand.csv                  : candidate SNPs with top correlated predictor
#   - screeplot.pdf
#   - loadings_histograms.pdf
#   - rda_group_sites.pdf / .png
#   - rda_plot1.pdf / .png      : RDA1 x RDA2 candidate SNP plot
#   - rda_plot2.pdf / .png      : RDA1 x RDA3 candidate SNP plot
# ============================================================


# ---- 0. Initialization -----------------------------------------------

sink("log.txt", split = TRUE)
options(echo = TRUE)

library(adegenet)
library(psych)
library(vegan)
library(data.table)


# ---- 1. Load genotype matrix -----------------------------------------
# NOTE: fread() does not natively support rownames; the rownames=1 argument
# is passed to as.matrix() but has no effect here — row names must be
# assigned manually after reading.

gen <- as.matrix(fread("snp.csv"), rownames = 1)

# Check dimensions and missing data
dim(gen)
sum(is.na(gen))


# ---- 2. Load environmental data --------------------------------------

env <- read.csv("env.csv")
str(env)


# ---- 3. Harmonize individual IDs -------------------------------------
# fread() lowercases column names; apply the same transformation to
# individual IDs in env to ensure a consistent join key.

env$individual <- as.character(env$individual)
identical(rownames(gen), env[, 1])   # diagnostic check

rownames(gen) <- tolower(rownames(gen))
env[, 1]      <- tolower(as.character(env[, 1]))

identical(rownames(gen), env[, 1])   # should now be TRUE


# ---- 4. Select predictor variables -----------------------------------
# Variables below are excluded due to multicollinearity or low
# biological relevance for this study system.
# The retained columns are assumed to occupy positions 6-23 in the
# subsetted data frame — verify this against your own env.csv structure.

pred <- subset(env, select = -c(
  BIO1, BIO4, BIO5, BIO6, BIO7, BIO8, BIO9,
  BIO12, BIO15, BIO16, BIO17, BIO18,
  bdod15, bdod30, bdod60,
  cfvo15, cfvo30, cfvo60,
  phh3o15, phh4o30, phh5o60,
  clay30, clay60,
  silt15, silt30, silt60,
  soc5, soc15, soc30,
  nitrogen15, nitrogen30, nitrogen60,
  ocd5, ocd15,
  sand5, sand15, sand30, sand60,
  UVB2, UVB3, UVB4, UVB5, UVB6
))

pred1 <- pred[, 6:23]
colnames(pred1) <- paste0("eco", 1:18)


# ---- 5. Fit RDA model ------------------------------------------------
# scale = TRUE standardizes the response matrix (SNPs).
# All 18 eco variables are used as constraints.

rda_model <- rda(gen ~ ., data = pred1, scale = TRUE)
print(rda_model)

r_squared_adj <- RsquareAdj(rda_model)$adj.r.squared
print(paste("Adjusted R-squared:", r_squared_adj))

eigen_summary <- summary(eigenvals(rda_model, model = "constrained"))
print(eigen_summary)

output_data <- data.frame(
  Adj_R_Squared           = r_squared_adj,
  Eigenvalues_Constrained = eigen_summary,
  Proportion_Explained    = NA,
  Cumulative_Proportion   = NA
)
write.csv(output_data, file = "rda_results.csv", row.names = FALSE)


# ---- 6. Screeplot ----------------------------------------------------

pdf("screeplot.pdf")
screeplot(rda_model)
dev.off()


# ---- 7. Permutation significance test --------------------------------
# Default: 999 permutations. parallel uses all available cores.

signif.full <- anova.cca(rda_model, parallel = getOption("mc.cores"))
signif.full


# ---- 8. Variance inflation factors (VIF) ----------------------------
# VIF > 10 indicates strong collinearity; consider removing that variable
# and re-fitting before using the model for candidate SNP detection.

vif_values <- vif.cca(rda_model)
vif_data   <- data.frame(Variable = names(vif_values), VIF = vif_values)
write.csv(vif_data, "vif_values.csv", row.names = FALSE)


# ---- 9. RDA biplot — individuals colored by population group --------

env$group <- as.factor(env$group)
bg <- c("#ff7f00", "#1f78b4", "#ffff33", "#33a02c", "#6a3d9a")
bg <- rep(bg, length.out = length(env$group))

pdf("rda_group_sites.pdf", width = 8, height = 8)
plot(rda_model, type = "n", scaling = 3)
points(rda_model, display = "species", pch = 20, cex = 0.7, col = "gray32", scaling = 3)
points(rda_model, display = "sites",   pch = 21, cex = 1.3, col = "gray32", scaling = 3, bg = bg)
text(rda_model,   scaling = 3, display = "bp", col = "#0868ac", cex = 1)
legend("bottomright", legend = levels(env$group),
       bty = "n", col = "gray32", pch = 21, cex = 1, pt.bg = bg)
dev.off()


# ---- 10. High-resolution PNG output ----------------------------------
# Update file_path to match your local directory structure before running.

library(png)
file_path <- "rda_group_sites.png"
png(file_path, width = 8200, height = 8200, units = "px", res = 300)
plot(rda_model, type = "n", scaling = 3)
points(rda_model, display = "species", pch = 20, cex = 0.7, col = "gray32", scaling = 3)
points(rda_model, display = "sites",   pch = 21, cex = 1.3, col = "gray32", scaling = 3, bg = bg)
text(rda_model,   scaling = 3, display = "bp", col = "#0868ac", cex = 1)
legend("bottomright", legend = levels(env$group),
       bty = "n", col = "gray32", pch = 21, cex = 1, pt.bg = bg)
dev.off()


# ---- 11. Extract SNP loadings and identify outlier candidates --------
# Outliers are defined as SNPs whose loading on a given RDA axis exceeds
# ±3 standard deviations from the mean loading on that axis.

load.rda <- scores(rda_model, choices = c(1:3), display = "species")
write.csv(load.rda, file = "load_rda_species_scores.csv")

pdf("loadings_histograms.pdf")
par(mfrow = c(1, 3))
hist(load.rda[, 1], main = "Loadings on RDA1")
hist(load.rda[, 2], main = "Loadings on RDA2")
hist(load.rda[, 3], main = "Loadings on RDA3")
dev.off()

outliers <- function(x, z) {
  lims <- mean(x) + c(-1, 1) * z * sd(x)
  x[x < lims[1] | x > lims[2]]
}

cand1 <- outliers(load.rda[, 1], 3)
cand2 <- outliers(load.rda[, 2], 3)
cand3 <- outliers(load.rda[, 3], 3)


# ---- 12. Correlate candidate SNPs with environmental predictors ------

ncand <- length(cand1) + length(cand2) + length(cand3)
ncand

cand1 <- cbind.data.frame(rep(1, length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2, length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3, length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis", "snp", "loading")

cand     <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

foo <- matrix(nrow = ncand, ncol = 18)
colnames(foo) <- paste0("eco", 1:18)

for (i in 1:length(cand$snp)) {
  nam       <- cand[i, 2]
  snp.gen   <- gen[, nam]
  foo[i, ]  <- apply(pred1, 2, function(x) cor(x, snp.gen))
}

cand <- cbind.data.frame(cand, foo)
head(cand)

# Remove SNPs that appear on multiple axes — keep first occurrence only
length(cand$snp[duplicated(cand$snp)])
foo <- cbind(cand$axis, duplicated(cand$snp))
table(foo[foo[, 1] == 1, 2])
table(foo[foo[, 1] == 2, 2])
table(foo[foo[, 1] == 3, 2])
cand <- cand[!duplicated(cand$snp), ]

# Assign each candidate SNP to its most correlated predictor
for (i in 1:length(cand$snp)) {
  bar         <- cand[i, ]
  cand[i, 22] <- names(which.max(abs(bar[4:21])))
  cand[i, 23] <- max(abs(bar[4:21]))
}
colnames(cand)[22] <- "predictor"
colnames(cand)[23] <- "correlation"
table(cand$predictor)


# ---- 13. Candidate SNP visualization --------------------------------
# Each color represents a different ecological predictor variable.
# SNPs not identified as candidates are rendered transparent.
# NOTE: col.pred relies on grep("chr", ...) to detect non-candidate SNPs.
# If your SNP IDs do not contain "chr", all points will be invisible —
# update the pattern to match your naming convention.

sel     <- cand$snp
env_col <- cand$predictor

color_map <- c(
  eco1  = '#1f78b4', eco2  = '#a6cee3', eco3  = '#6a3d9a', eco4  = '#e31a1c',
  eco5  = '#ffff33', eco6  = '#BCA1F5', eco7  = '#fb9a99', eco8  = '#33a02c',
  eco9  = '#ff7f00', eco10 = '#fdbf6f', eco11 = '#cab2d6', eco12 = '#6b3d9a',
  eco13 = '#8dd3c7', eco14 = '#fb8072', eco15 = '#80b1d3', eco16 = '#fdb462',
  eco17 = '#1DD9F5', eco18 = '#b2df8a'
)
env_col <- color_map[env_col]

col.pred <- rownames(rda_model$CCA$v)
for (i in 1:length(sel)) {
  foo          <- match(sel[i], col.pred)
  col.pred[foo] <- env_col[i]
}
col.pred[grep("chr", col.pred)] <- '#f1eef6'
empty         <- col.pred
empty[grep("#f1eef6", empty)] <- rgb(0, 1, 0, alpha = 0)
empty.outline <- ifelse(empty == "#00FF0000", "#00FF0000", "gray32")
bg            <- unname(color_map)

# RDA1 x RDA2
pdf("rda_plot1.pdf", width = 8, height = 8)
plot(rda_model,  type = "n",       scaling = 3, xlim = c(-1, 1), ylim = c(-1, 1))
points(rda_model, display = "species", pch = 21, cex = 1, col = "gray32",      bg = col.pred,    scaling = 3)
points(rda_model, display = "species", pch = 21, cex = 1, col = empty.outline, bg = empty,        scaling = 3)
text(rda_model,  scaling = 3, display = "bp", col = "#0868ac", cex = 1)
legend("bottomright", legend = paste0("eco", 1:18),
       bty = "n", col = "gray32", pch = 21, cex = 0.5, pt.bg = bg)
dev.off()

file_path <- "rda_plot1.png"
png(file_path, width = 8200, height = 8200, units = "px", res = 300)
plot(rda_model,  type = "n",       scaling = 3, xlim = c(-1, 1), ylim = c(-1, 1))
points(rda_model, display = "species", pch = 21, cex = 1, col = "gray32",      bg = col.pred,    scaling = 3)
points(rda_model, display = "species", pch = 21, cex = 1, col = empty.outline, bg = empty,        scaling = 3)
text(rda_model,  scaling = 3, display = "bp", col = "#0868ac", cex = 1)
legend("bottomright", legend = paste0("eco", 1:18),
       bty = "n", col = "gray32", pch = 21, cex = 0.5, pt.bg = bg)
dev.off()

# RDA1 x RDA3
pdf("rda_plot2.pdf", width = 8, height = 8)
plot(rda_model,  type = "n",       scaling = 3, xlim = c(-1, 1), ylim = c(-1, 1), choices = c(1, 3))
points(rda_model, display = "species", pch = 21, cex = 1, col = "gray32",      bg = col.pred,    scaling = 3, choices = c(1, 3))
points(rda_model, display = "species", pch = 21, cex = 1, col = empty.outline, bg = empty,        scaling = 3, choices = c(1, 3))
text(rda_model,  scaling = 3, display = "bp", col = "#0868ac", cex = 1, choices = c(1, 3))
legend("bottomright", legend = paste0("eco", 1:18),
       bty = "n", col = "gray32", pch = 21, cex = 0.5, pt.bg = bg)
dev.off()

file_path <- "rda_plot2.png"
png(file_path, width = 8200, height = 8200, units = "px", res = 300)
plot(rda_model,  type = "n",       scaling = 3, xlim = c(-1, 1), ylim = c(-1, 1), choices = c(1, 3))
points(rda_model, display = "species", pch = 21, cex = 1, col = "gray32",      bg = col.pred,    scaling = 3, choices = c(1, 3))
points(rda_model, display = "species", pch = 21, cex = 1, col = empty.outline, bg = empty,        scaling = 3, choices = c(1, 3))
text(rda_model,  scaling = 3, display = "bp", col = "#0868ac", cex = 1, choices = c(1, 3))
legend("bottomright", legend = paste0("eco", 1:18),
       bty = "n", col = "gray32", pch = 21, cex = 0.5, pt.bg = bg)
dev.off()


# ---- 14. Export final candidate SNP table ---------------------------

write.csv(cand, file = "cand.csv")
sink()
