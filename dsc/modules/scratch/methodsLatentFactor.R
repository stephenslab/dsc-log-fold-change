#' Estimate number of latent factors
#'
#' Use SVA package to estimate the number of latent factors.
#' This number is usually required for differential expression analysis
#' that accounts for latent factors. Most methods will estimate this
#' number for the user if not provided, in which case the computation will
#' take more time.

run_sva <- function(Y_log_normed, X) {
  library(sva)
  num_sv <- num.sv(Y_log_normed,X,method="leek")
  X_0 <- model.matrix(~1)
  res_sv <- sva(Y_log_normed, X, X_0, n.sv=num_sv)
  return(res_sv)
}


# estimate latent variables ----------


# estimate_latent_confound <- function() {
# }
# run_sva <- function(Y1, Y2) {
#
#   library(sva)
# #  library(BiocParallel)
#
#   Y <- cbind(as.matrix(Y1), as.matrix(Y2))
#
#   n1 <- dim(Y1)[2]
#   n2 <- dim(Y2)[2]
#   x <- rep(c(1,2), times = c(n1, n2))
#   x <- factor(x)
#
#   if (sum(duplicated(colnames(Y))) > 0) {
#     colnames(Y) <- paste0("cell.", c(1:ncol(Y))) }
#
#
#   dds <- DESeqDataSetFromMatrix(countData = round(Y),
#                                 colData = data.frame(condition = x),
#                                 design = ~condition)
#
#   dds <- estimateSizeFactors(dds,type="poscounts")
#   dds <- estimateDispersions(dds, minmu = 1e-3)
#   dds <- nbinomLRT(dds, minmu=1e-3, reduced=~1)
#   res <- results(dds, name="condition_2_vs_1")
#
#     # res <- results(dds, contrast = c("condition", levels(factor(x))[1],
#     #                                levels(factor(x))[2]), alpha = 0.05);
#   return(list(pval = res$pvalue,
#               est = res$log2FoldChange,
#               se = res$lfcSE))
# }
#
#
#
