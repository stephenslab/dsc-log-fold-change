zinbwave_edger <- function(Y1, Y2){
  library(edgeR)
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
  # use -edgeR
  d <- DGEList(Y)
  d <- edgeR::calcNormFactors(d)
  d$weights <- weights
  d <- estimateDisp(d, design)
  fit <- glmFit(d,design)
  lrt <- glmWeightedF(fit,coef=2, independentFiltering = TRUE)
  pvalues <- lrt$table$PValue

  return(list(pval = pvalues,
              betahat = lrt$table$logFC) )
}
