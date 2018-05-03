#' @param Y1 count matrix of group 1; p by n1
#' @param Y2 count matrix of group 2; p by n2

run_DESeq2 <- function(Y1, Y2) {
  Y <- cbind($(Y1), $(Y2));
  x <- rep(c(1,2), each = c(dim($(Y1))[2],dim($(Y2))[2]));
  dds <- DESeqDataSetFromMatrix(countData = round(Y),
                                colData = data.frame(condition = x),
                                design = ~x);
  dds <- DESeq(dds);
  res <- results(dds, contrast = c("condition", levels(factor(x))[1],
                                   levels(factor(x))[2]), alpha = 0.05);

  return(list(pval = res$pval,
              est = res$log2FoldChange,
              se = res$lfcSE))
}
