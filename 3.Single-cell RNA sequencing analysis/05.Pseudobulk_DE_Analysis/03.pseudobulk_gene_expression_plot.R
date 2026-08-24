# ============================================================
# Pseudobulk Gene Expression Plot (Violin + Boxplot + Jitter)
# Description: Generates publication-ready expression plots for
#              specific genes with significance brackets.
# Prerequisite: Objects from the main pseudobulk pipeline must exist
#               (deg_all, dge, sample_meta)
# ============================================================

# 0. Load packages
library(edgeR); library(dplyr); library(ggplot2); library(patchwork)

# ============================================================
# 1. Set DEG cutoff criteria (FDR < 0.05, |log2FC| > 0.5)
# ============================================================
fc_cutoff <- 0.5; fdr_cutoff <- 0.05

deg_all <- deg_all %>% dplyr::mutate(
  change = dplyr::case_when(
    FDR < fdr_cutoff & logFC >  fc_cutoff ~ "Up",
    FDR < fdr_cutoff & logFC < -fc_cutoff ~ "Down",
    TRUE ~ "NS"
  )
)
table(deg_all$change)

# ============================================================
# 2. Extract significant DEGs
# ============================================================
deg_sig <- deg_all %>%
  dplyr::filter(FDR < fdr_cutoff, abs(logFC) > fc_cutoff) %>%
  dplyr::arrange(FDR)
dim(deg_sig); head(deg_sig)

# ============================================================
# 3. Check statistics of target gene (DDIT4)
# ============================================================
target_genes <- c("DDIT4")
target_stats <- deg_all %>%
  dplyr::filter(gene %in% target_genes) %>%
  dplyr::select(gene, logFC, logCPM, F, PValue, FDR, change) %>%
  dplyr::arrange(match(gene, target_genes))
target_stats

write.csv(target_stats, "DDIT4_pseudobulk_statistics.csv", row.names = FALSE)

# ============================================================
# 4. Animal-level pseudobulk normalized expression
# ============================================================
logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
dim(logCPM)

# ============================================================
# 5. Check if target genes exist in logCPM
# ============================================================
target_genes %in% rownames(logCPM)

# ============================================================
# 6. Color scheme (TB: teal / LY: orange-red)
# ============================================================
group_cols <- c("TB" = "#22BFC4", "LY" = "#F8766D")

# ============================================================
# 7. Significance stars (|logFC| > 0.5 + FDR < 0.05, then by FDR)
# ============================================================
get_significance <- function(logfc, fdr) {
  if (is.na(logfc) || is.na(fdr) || abs(logfc) <= fc_cutoff || fdr >= fdr_cutoff) return("ns")
  if (fdr < 0.0001)     return("****")
  else if (fdr < 0.001) return("***")
  else if (fdr < 0.01)  return("**")
  else                  return("*")
}

# ============================================================
# 8. General plotting function (trim=TRUE, violin stays within axis)
# ============================================================
plot_pseudobulk_gene <- function(gene_name) {
  if (!gene_name %in% rownames(logCPM)) {
    stop(paste0(gene_name, " is not present in logCPM."))
  }

  df_gene <- data.frame(
    animal_id = colnames(logCPM),
    expression = as.numeric(logCPM[gene_name, ]),
    stringsAsFactors = FALSE
  ) %>% dplyr::left_join(sample_meta %>% dplyr::select(animal_id, group), by = "animal_id")

  df_gene$group <- factor(df_gene$group, levels = c("TB", "LY"))

  stat_gene <- deg_all %>% dplyr::filter(gene == gene_name)
  if (nrow(stat_gene) == 0) stop(paste0(gene_name, " is not present in deg_all."))
  logfc <- stat_gene$logFC[1]; fdr <- stat_gene$FDR[1]
  sig_label <- get_significance(logfc, fdr)

  ymax <- max(df_gene$expression, na.rm = TRUE)
  ymin <- min(df_gene$expression, na.rm = TRUE)
  yrange <- ymax - ymin
  if (yrange == 0) yrange <- 1

  bracket_y <- ymax + 0.15 * yrange
  star_y    <- ymax + 0.23 * yrange
  info_y    <- ymax + 0.42 * yrange
  upper_y   <- ymax + 0.57 * yrange

  p <- ggplot(df_gene, aes(x = group, y = expression, fill = group, color = group)) +
    geom_violin(trim = TRUE, alpha = 0.75, linewidth = 0.45) +
    geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.45) +
    geom_jitter(width = 0.06, size = 2.5, alpha = 0.9) +
    scale_fill_manual(values = group_cols) +
    scale_color_manual(values = group_cols) +
    annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y,
             linewidth = 0.5, color = "black") +
    annotate("segment", x = 1, xend = 1, y = bracket_y, yend = bracket_y - 0.05 * yrange,
             linewidth = 0.5, color = "black") +
    annotate("segment", x = 2, xend = 2, y = bracket_y, yend = bracket_y - 0.05 * yrange,
             linewidth = 0.5, color = "black") +
    annotate("text", x = 1.5, y = star_y, label = sig_label, size = 5, fontface = "bold") +
    annotate("text", x = 1.5, y = info_y,
             label = paste0("log2FC = ", sprintf("%.2f", logfc),
                            "\nFDR = ", formatC(fdr, format = "e", digits = 2)),
             size = 3.4) +
    labs(x = NULL, y = "Expression level (logCPM)", title = gene_name) +
    coord_cartesian(ylim = c(ymin - 0.05 * yrange, upper_y)) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "italic", hjust = 0.5, size = 12),
      axis.title.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 9),
      axis.line = element_line(linewidth = 0.5),
      plot.margin = margin(10, 10, 5, 5)
    )

  return(p)
}

# ============================================================
# 9. Plot DDIT4
# ============================================================
p_DDIT4 <- plot_pseudobulk_gene("DDIT4")
p_DDIT4
