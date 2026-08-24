# ============================================================
# 0. Load packages
# ============================================================
library(Seurat); library(Matrix); library(edgeR)
library(dplyr); library(tibble); library(ggplot2); library(ggrepel)

# ============================================================
# 1. Check metadata
# ============================================================
colnames(seu_obj@meta.data)
table(seu_obj$orig.ident, seu_obj$group)

# Number of animals per group (Expected: TB = 6, LY = 5)
length(unique(seu_obj$orig.ident[seu_obj$group == "TB"]))
length(unique(seu_obj$orig.ident[seu_obj$group == "LY"]))

# ============================================================
# 2. Extract Epithelium cells
# ============================================================
seu_sub <- subset(seu_obj, subset = celltype == "Epithelium")
table(seu_sub$orig.ident)
table(seu_sub$orig.ident, seu_sub$group)

# ============================================================
# 3. Extract RNA raw counts
# ============================================================
DefaultAssay(seu_sub) <- "RNA"
counts <- GetAssayData(seu_sub, assay = "RNA", layer = "counts") # For Seurat v5
dim(counts)

# Check if counts are integers (sum should be 0)
summary(counts@x)
sum(counts@x %% 1 != 0) 

# ============================================================
# 4. Create cell-level metadata
# ============================================================
meta <- seu_sub@meta.data %>% tibble::rownames_to_column("cell")
meta$animal_id <- meta$orig.ident
table(meta$animal_id, meta$group)

# ============================================================
# 5. Construct pseudobulk counts
# Note: Directly use all Epithelium raw counts. No pre-filtering.
# ============================================================
animal_ids <- unique(meta$animal_id)
pb_counts <- sapply(animal_ids, function(id) {
  cells_i <- meta$cell[meta$animal_id == id]
  Matrix::rowSums(counts[, cells_i, drop = FALSE])
})
pb_counts <- as.matrix(pb_counts)
dim(pb_counts) # Expected: genes x 11 animals
colnames(pb_counts)

# ============================================================
# 6. Create animal-level metadata
# ============================================================
sample_meta <- meta %>% dplyr::select(animal_id, group) %>% dplyr::distinct()
sample_meta <- sample_meta[match(colnames(pb_counts), sample_meta$animal_id), ]
rownames(sample_meta) <- sample_meta$animal_id
sample_meta

# ============================================================
# 7. Check pseudobulk
# ============================================================
ncol(pb_counts) # Expected = 11
table(sample_meta$group) # Expected: LY = 5, TB = 6
identical(colnames(pb_counts), rownames(sample_meta)) # Must be TRUE

# ============================================================
# 8. Epithelium cell number per animal
# ============================================================
cell_number <- meta %>% dplyr::count(animal_id, group, name = "n_cells")
cell_number
write.csv(cell_number, "Epithelium_cell_number_per_animal.csv", row.names = FALSE)

# ============================================================
# 9. Set factor levels for groups
# ============================================================
sample_meta$group <- factor(sample_meta$group, levels = c("LY", "TB")) # LY = reference
table(sample_meta$group)

# ============================================================
# 10. Create edgeR DGEList
# ============================================================
dge <- edgeR::DGEList(counts = pb_counts, group = sample_meta$group)

# ============================================================
# 11. filterByExpr (Low-expression filtering)
# ============================================================
keep <- edgeR::filterByExpr(dge, group = sample_meta$group)
table(keep)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dim(dge)

# ============================================================
# 12. TMM normalization
# ============================================================
dge <- edgeR::calcNormFactors(dge)
dge$samples

# ============================================================
# 13. Design matrix
# ============================================================
design <- model.matrix(~ group, data = sample_meta)
design
# Since LY is reference: groupTB > 0 (TB > LY), groupTB < 0 (TB < LY)

# ============================================================
# 14. Estimate dispersion
# ============================================================
dge <- edgeR::estimateDisp(dge, design)
plotBCV(dge)

# ============================================================
# 15. Quasi-likelihood model
# ============================================================
fit <- edgeR::glmQLFit(dge, design)
qlf <- edgeR::glmQLFTest(fit, coef = "groupTB")

# ============================================================
# 16. Extract all results
# ============================================================
deg_all <- edgeR::topTags(qlf, n = Inf, adjust.method = "BH")$table %>% tibble::rownames_to_column("gene")
head(deg_all) # Includes: gene, logFC, logCPM, F, PValue, FDR

# ============================================================
# 17. DEG cutoff criteria
# ============================================================
fc_cutoff <- log2(1.5)
deg_all <- deg_all %>% dplyr::mutate(
  change = dplyr::case_when(
    FDR < 0.05 & logFC >= fc_cutoff ~ "Up",
    FDR < 0.05 & logFC <= -fc_cutoff ~ "Down",
    TRUE ~ "NS"
  )
)
table(deg_all$change)

# ============================================================
# 18. Extract significant DEGs
# ============================================================
deg_sig <- deg_all %>% dplyr::filter(FDR < 0.05, abs(logFC) >= fc_cutoff)
dim(deg_sig)
head(deg_sig)

# ============================================================
# 19. Save results
# ============================================================
write.csv(deg_all, "Epithelium_TB_vs_LY_pseudobulk_edgeR_all_genes.csv", row.names = FALSE)
write.csv(deg_sig, "Epithelium_TB_vs_LY_pseudobulk_edgeR_DEGs.csv", row.names = FALSE)

# ============================================================
# 20. Check ASXL1
# ============================================================
deg_all %>% dplyr::filter(gene == "ASXL1")

# ============================================================
# 21. Check all target genes
# ============================================================
target_genes <- c("GLUL", "ASXL1", "EGLN1", "DDIT4", "LIMD1", "ROCK2", 
                  "DNM3", "ATP8A", "PRKN", "PLD1", "SLC25A13", "EPAS1")
target_result <- deg_all %>% 
  dplyr::filter(gene %in% target_genes) %>% 
  dplyr::arrange(match(gene, target_genes))
target_result

# ============================================================
# 22. Save target genes
# ============================================================
write.csv(target_result, "Epithelium_target_genes_pseudobulk.csv", row.names = FALSE)

# ============================================================
# 23. Volcano plot
# ============================================================
deg_all <- deg_all %>% dplyr::mutate(
  plot_FDR = pmax(FDR, 1e-300),
  minus_log10_FDR = -log10(plot_FDR)
)

top_up <- deg_all %>% dplyr::filter(change == "Up") %>% dplyr::arrange(FDR, desc(logFC)) %>% head(10)
top_down <- deg_all %>% dplyr::filter(change == "Down") %>% dplyr::arrange(FDR, logFC) %>% head(10)
top_label <- dplyr::bind_rows(top_up, top_down)

p_volcano <- ggplot(deg_all, aes(x = logFC, y = minus_log10_FDR)) +
  geom_point(aes(shape = change), alpha = 0.7, size = 1.5) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  ggrepel::geom_text_repel(data = top_label, aes(label = gene), size = 3, max.overlaps = 20) +
  labs(x = "log2 fold change (TB vs LY)", y = "-log10(FDR)", shape = NULL) +
  theme_classic()

p_volcano
ggsave("Epithelium_TB_vs_LY_pseudobulk_volcano.pdf", p_volcano, width = 6, height = 5)

# ============================================================
# 24. Pseudobulk normalized expression
# ============================================================
logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
dim(logCPM)

# ============================================================
# 25. Check ASXL1 expression per animal
# ============================================================
logCPM["ASXL1", ]

# ============================================================
# 26. ASXL1 animal-level expression
# ============================================================
gene_to_plot <- "ASXL1"
df_gene <- data.frame(
  animal_id = colnames(logCPM),
  expression = as.numeric(logCPM[gene_to_plot, ])
) %>% dplyr::left_join(sample_meta %>% dplyr::select(animal_id, group), by = "animal_id")
df_gene

p_gene <- ggplot(df_gene, aes(x = group, y = expression)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 3) +
  labs(x = NULL, y = "logCPM", title = gene_to_plot) +
  theme_classic()

p_gene
