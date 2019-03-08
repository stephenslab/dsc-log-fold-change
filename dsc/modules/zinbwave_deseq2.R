run_zinbwave_deseq2 <- function(Y1, Y2){
  library(DESeq2)
  library(zinbwave)

  Y <- cbind(as.matrix(Y1), as.matrix(Y2))

  n1 <- dim(Y1)[2]
  n2 <- dim(Y2)[2]
  x <- rep(c(1,2), times = c(n1, n2))
  x <- factor(x)

  if (sum(duplicated(colnames(Y))) > 0) {
    colnames(Y) <- paste0("cell.", c(1:ncol(Y))) }
  if (is.null(rownames(Y))) {
    rownames(Y) <- paste0("gene.", c(1:nrow(Y))) }

  condition <- x
  design <- model.matrix(~ condition)

  # compute zinbwave weights
  zinb <- zinbFit(Y, X = design, epsilon = 1e12)
  weights <- computeObservationalWeights(zinb, Y)

  # use DESeq2
  d <- DESeqDataSetFromMatrix(Y,
                              colData= DataFrame(data.frame(condition=condition)),
                              design= ~ condition)
  d <- estimateSizeFactors(d, type="poscounts")
  dimnames(weights) = NULL
  assays(d)[["weights"]] = weights
  d <- estimateDispersions(d, minmu = 1e-3)
  #dse = nbinomWaldTest(dse, betaPrior=TRUE, useT=TRUE, df=rowSums(weights)-2, minmu = 1e-3)
  d <- nbinomLRT(d, minmu=1e-3, reduced=~1)
  res <- results(d, name="condition_2_vs_1")

  return(list(pval = res$pvalue,
              betahat = res$log2FoldChange) )
}
