# ============================================================
# Diagonal Volcano Plot for scRNA-seq DEGs
# Description: Generates a diagonal volcano plot (x = signed -log10(FDR), 
#              y = logFC) with manual gene labeling.
# ============================================================

# 0. Load packages
library(ggplot2); library(ggrepel); library(dplyr)

# 1. Read data
deg <- read.csv("ScRNA_DEGS.regulation.csv", stringsAsFactors = FALSE)

# 2. Grouping and diagonal axis transformation
deg$change <- factor(tolower(deg$Regulation), levels = c("down", "ns", "up"))
deg$change <- droplevels(deg$change)
deg$Difference <- sign(deg$logFC) * -log10(deg$FDR)
print(table(deg$change))

# Legend labels and colors
lab_map <- c(down = paste0("down(", sum(deg$change == "down"), ")"),
             ns = "ns", up = paste0("up(", sum(deg$change == "up"), ")"))
col_map <- c(down = "#0d1b46", ns = "grey", up = "tomato")

# 3. Gene labels (15 genes, unified black style) + manual blank positions
label_pos <- data.frame(
  gene = c("ARPC5L", "ATP8A1", "DNM3", "DOCK10", "DDIT4", "DDIT4", "NRCAM",
           "PLD1", "PRKN", "PSMD5", "RPL35", "SERGEF", "SLC25A13", "TEP1", "WDR38"),
  label_x = c(1.5, -9.3, -3.8, 3.0, 6.0, 6.5, -1.2, -2.6, -8.0, -2.0, -0.5, -0.6, -5.8, 0.8, 5.0),
  label_y = c(2.5, -3.8, -6.2, 3.5, 12.0, -2.0, -8.6, -4.4, -5.4, 2.0, 1.5, -5.2, -4.6, -7.6, 7.5)
)

label_data <- deg %>%
  filter(gene %in% label_pos$gene) %>%
  distinct(gene, .keep_all = TRUE) %>%
  left_join(label_pos, by = "gene")

print(label_data$gene) # Expected: 15 genes

# 4. Diagonal volcano plot
volcano_plot <- ggplot(deg, aes(x = Difference, y = logFC, color = change)) +
  geom_point(size = 1.2) +
  scale_color_manual(NULL, values = col_map[levels(deg$change)], labels = lab_map[levels(deg$change)]) +
  # Indicator lines (unified black)
  geom_segment(data = label_data, aes(x = label_x, y = label_y, xend = Difference, yend = logFC),
               color = "black", linewidth = 0.3) +
  # Gene names (unified black, same size)
  geom_text(data = label_data, aes(x = label_x, y = label_y, label = gene),
            color = "black", size = 3.5) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Signed -log10(FDR)", y = "Log-Fold change") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = c(0.72, 0.25), legend.justification = c(0.5, 0.5),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 10), aspect.ratio = 1
  )

print(volcano_plot)

# 5. Save as square PDF
ggsave("ScRNA_DEGS.regulation.pdf", volcano_plot, width = 8, height = 8, bg = "white")
