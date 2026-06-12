#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
GENE <- ifelse(length(args) >= 1, args[1], "BMPR2")

TPM_FILE   <- "TPM_matrix.csv"
GROUP_FILE <- "sample_group.tsv"
GTF_FILE   <- "00.index/Sus_scrofa.Sscrofa11.1.111.chr.gtf"  
n_perm <- 50
n_each <- 5
set.seed(42)

tpm <- read.csv(TPM_FILE, row.names = 1, check.names = FALSE)
grp <- read.delim(GROUP_FILE, stringsAsFactors = FALSE)

row_id <- NA
if (GENE %in% rownames(tpm)) {
    row_id <- GENE
} else if (file.exists(GTF_FILE)) {
    hit <- suppressWarnings(system(
        sprintf("grep -m1 'gene_name \"%s\"' %s", GENE, GTF_FILE), intern = TRUE))
    if (length(hit) > 0) row_id <- sub('.*gene_id "([^"]+)".*', "\\1", hit[1])
}
if (is.na(row_id) || !(row_id %in% rownames(tpm)))
    stop("not found: ", GENE)

vals  <- as.numeric(tpm[row_id, ])
g_of  <- grp$Group[match(colnames(tpm), grp$Sample)]
L_vals <- vals[which(g_of == "L")]
T_vals <- vals[which(g_of == "T")]
n_each <- min(n_each, length(L_vals), length(T_vals))
full_test <- wilcox.test(L_vals, T_vals)

all_samples <- vector("list", n_perm)
perm <- data.frame(iter = integer(n_perm), p_value = numeric(n_perm),
                   median_L = numeric(n_perm), median_T = numeric(n_perm),
                   diff_median = numeric(n_perm))
for (i in seq_len(n_perm)) {
    sub_L <- sample(L_vals, n_each)
    sub_T <- sample(T_vals, n_each)
    perm$iter[i]        <- i
    perm$p_value[i]     <- wilcox.test(sub_L, sub_T)$p.value
    perm$median_L[i]    <- median(sub_L)
    perm$median_T[i]    <- median(sub_T)
    perm$diff_median[i] <- median(sub_L) - median(sub_T)
    all_samples[[i]] <- data.frame(iter = i,
                                   Group = rep(c("L", "T"), each = n_each),
                                   value = c(sub_L, sub_T))
}
perm$sig  <- ifelse(perm$p_value < 0.05, "Sig (p<0.05)", "Not sig")
sig_count <- sum(perm$p_value < 0.05)
all_df <- merge(do.call(rbind, all_samples), perm[, c("iter", "sig")], by = "iter")
col_sig <- c("Not sig" = "#4E79A7", "Sig (p<0.05)" = "#E15759")

pair_df <- data.frame(Group = rep(c("L", "T"), each = n_perm),
                      median = c(perm$median_L, perm$median_T),
                      sig = rep(perm$sig, 2))
p1 <- ggplot() +
    geom_segment(data = perm, aes(x = "L", xend = "T", y = median_L,
                 yend = median_T, color = sig), alpha = 0.5, linewidth = 0.6) +
    geom_point(data = pair_df, aes(x = Group, y = median, color = sig),
               size = 2, alpha = 0.7) +
    scale_color_manual(values = c("Not sig" = "grey60",
                       "Sig (p<0.05)" = "#E15759"), name = "") +
    labs(title = paste0(n_perm, " resamples: median of L vs T"),
         subtitle = "Each line = one resample; slope = direction of difference",
         y = paste(GENE, "median (TPM)")) +
    theme_bw(base_size = 13) + theme(legend.position = "bottom")

pick <- sort(c(head(perm$iter[perm$sig == "Sig (p<0.05)"], 2),
               head(perm$iter[perm$sig == "Not sig"], 2)))
sub_df <- all_df[all_df$iter %in% pick, ]
sub_df$label <- paste0("Iter ", sub_df$iter, " (", sub_df$sig, ")")
sub_df$label <- factor(sub_df$label, levels = unique(sub_df$label[order(sub_df$iter)]))
p2 <- ggplot(sub_df, aes(x = Group, y = value, fill = Group)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.1, size = 2.5, alpha = 0.8) +
    facet_wrap(~label, nrow = 1) +
    scale_fill_manual(values = c("L" = "#4E79A7", "T" = "#E15759")) +
    labs(title = "Different samples, different conclusions",
         subtitle = "Same gene, same groups -- result changes with sample selection",
         y = paste(GENE, "(TPM)")) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none", strip.text = element_text(size = 10))

perm_s <- perm[order(perm$diff_median), ]; perm_s$rank <- seq_len(n_perm)
p3 <- ggplot(perm_s, aes(x = rank, y = diff_median, fill = sig)) +
    geom_col(width = 0.8, alpha = 0.85) +
    geom_hline(yintercept = 0, linewidth = 0.8) +
    scale_fill_manual(values = col_sig, name = "") +
    labs(title = "Resamples sorted by effect size (median L - T)",
         subtitle = paste0("Bars on both sides = unstable direction | Sig: ",
                           sig_count, "/", n_perm),
         x = "Resample (sorted)", y = "Median difference (L - T)") +
    theme_bw(base_size = 13) + theme(legend.position = "bottom")

perm_p <- perm[order(perm$p_value), ]; perm_p$cum_frac <- seq_len(n_perm) / n_perm
p4 <- ggplot(perm_p, aes(x = p_value, y = cum_frac)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_step(color = "#4E79A7", linewidth = 1) +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.6) +
    annotate("text", x = 0.5, y = 0.3, size = 3.5, color = "grey40",
             label = "Dashed = uniform distribution\n(expected if no real difference)") +
    labs(title = "Cumulative distribution of resample p-values",
         subtitle = "Curve near the diagonal = no consistent difference",
         x = "p-value", y = "Cumulative fraction") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) + theme_bw(base_size = 13)

combined <- (p1 | p4) / p2 / p3 +
    plot_annotation(
        title = paste0(GENE, ": is the group difference real or sampling-dependent?"),
        subtitle = paste0("Full-sample Wilcoxon p = ", signif(full_test$p.value, 3),
                          " | ", sig_count, "/", n_perm, " resamples reached p < 0.05"),
        theme = theme(plot.title = element_text(size = 16, face = "bold"),
                      plot.subtitle = element_text(size = 12, color = "grey30")))

out_pdf <- paste0(GENE, "_resampling.pdf")
ggsave(out_pdf, combined, width = 14, height = 16, device = cairo_pdf)
cat("save:", out_pdf, "\n")
