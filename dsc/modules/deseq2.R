run_deseq2 <- function(Y, X) {

  library(DESeq2)
  #  library(BiocParallel)

  #  Y <- cbind(as.matrix(Y1), as.matrix(Y2))

  # n1 <- dim(Y1)[2]
  # n2 <- dim(Y2)[2]
  #  x <- rep(c(1,2), times = c(n1, n2))


  if (sum(duplicated(colnames(Y))) > 0) {
    colnames(Y) <- paste0("cell.", c(1:ncol(Y))) }

  rownames(X) <- colnames(Y)
  x <- factor(X[,2])

  dds <- DESeqDataSetFromMatrix(countData = round(Y),
                                colData = data.frame(condition = x),
                                design = ~condition)

  dds <- estimateSizeFactors(dds,type="poscounts")
  dds <- estimateDispersions(dds, minmu = 1e-3)
  dds <- nbinomLRT(dds, minmu=1e-3, reduced=~1)
  res <- results(dds, name="condition_1_vs_0")

  # res <- results(dds, contrast = c("condition", levels(factor(x))[1],
  #                                levels(factor(x))[2]), alpha = 0.05);
  return(list(pval = res$pvalue,
              est = res$log2FoldChange,
              se = res$lfcSE))
}
