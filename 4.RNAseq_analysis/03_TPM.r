sample_list_path <- "list"  
if (!file.exists(sample_list_path)) {
  stop("error：", sample_list_path)
}


samples <- readLines(sample_list_path)
if (length(samples) == 0) {
  stop("error")
}


counts_dir <- "."  
full_filenames <- file.path(counts_dir, paste0(samples, ".counts.txt"))


missing_files <- full_filenames[!file.exists(full_filenames)]
if (length(missing_files) > 0) {
  stop("error：\n", paste(missing_files, collapse="\n"))
}


first_file <- full_filenames[1]
first_df <- read.table(
  first_file, 
  header=TRUE, 
  row.names=1, 
  comment.char="#",  
  sep="\t",         
  check.names=FALSE  
)


if (ncol(first_df) < 6) {
  stop("error：", first_file, "\n：", ncol(first_df), 
       "\n")
}


gene_length <- first_df$Length


count_list <- lapply(1:length(full_filenames), function(i) {
  file <- full_filenames[i]
  sample_name <- samples[i]  
  
  
  df <- read.table(
    file, 
    header=TRUE, 
    row.names=1, 
    comment.char="#", 
    sep="\t", 
    check.names=FALSE
  )
  
  
  sample_data <- df[, 6, drop=FALSE]  
  colnames(sample_data) <- sample_name  
  
  return(sample_data)
})

expr_matrix <- do.call(cbind, count_list)


keep_genes <- rowSums(expr_matrix) > 0
expr_matrix_filtered <- expr_matrix[keep_genes, , drop=FALSE]
gene_length_filtered <- gene_length[keep_genes]


cat("rawgene：", nrow(expr_matrix), "\n")
cat("filter：", nrow(expr_matrix_filtered), "\n")

output_file <- "expression_matrix_filtered.csv"  
write.csv(expr_matrix_filtered, output_file, quote=FALSE)
cat("success：", output_file, "\n")


length_file <- "gene_length_filtered.txt"
write.table(data.frame(Geneid=names(gene_length_filtered), Length=gene_length_filtered),
            length_file, sep="\t", quote=FALSE, row.names=FALSE)
cat("save：", length_file, "\n")
