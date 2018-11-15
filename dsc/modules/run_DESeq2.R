run_DESeq2 <- function(Y1, Y2) {
  Y <- cbind(as.matrix(Y1), as.matrix(Y2))
  n1 <- dim(Y1)[2]
  n2 <- dim(Y2)[2]
  x <- rep(c(1,2), times = c(n1, n2))
  x <- factor(x)

  data.is.simulated <- sum(duplicated(colnames(Y))) > 0

  if (data.is.simulated) { colnames(Y) <- paste0("cell.", c(1:ncol(Y)))}

  dds <- DESeqDataSetFromMatrix(countData = round(Y),
                                colData = data.frame(condition = x),
                                design = ~condition)

  ncores_default <- detectCores()
  if (ncores_default == 8) {
    register(MulticoreParam(4))
    } else {
    register(MulticoreParam(ncores_default))
    }
  dds <- DESeq(dds, parallel = TRUE, sftype="iterate")
  res <- results(dds, contrast = c("condition", levels(factor(x))[1],
                                   levels(factor(x))[2]), alpha = 0.05);

  return(list(pval = res$pval,
              est = res$log2FoldChange,
              se = res$lfcSE))
}
