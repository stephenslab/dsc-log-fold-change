
#' @param per_gene_zero_cell For each gene, the mininum number of samples with zero count.
#'        Default is non-zero count in at least 5 samples.

# filter_genes <- function(Y1, Y2, min_cell_detected) {
#   counts <- cbind(Y1, Y2)
#   genes_to_keep <- which(rowSums(counts > 0) >= per_gene_zero_cell)
#   return(genes_to_keep)
# }

filter_genes <- function(Y, min_cell_detected=5) {
#  counts <- cbind(Y1, Y2)
  keep_genes <- which(rowSums(Y > 0) >= min_cell_detected)
  return(keep_genes)
}

# get_data_celltypes = function(counts, meta, args) {
#   counts = as.matrix(counts)
#   group1 = counts[ , which(meta$celltype == args$groups[1])]
#   group2 = counts[ , which(meta$celltype == args$groups[2])]
#   return(list(x=sample_data(group1, args$n1, args$p),
#               y=sample_data(group2, args$n2, args$p)))
# }
