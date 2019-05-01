edger <- function(Y, X) {

  library(edgeR)

  #  Y <- cbind(as.matrix(Y1), as.matrix(Y2))
  #  n1 <- sum(X[,2]==1)
  #  n2 <- sum(X[,2]==0)
  #  x <- rep(c(1,2), times = c(n1, n2))
  #  x <- factor(x)

  if (sum(duplicated(colnames(Y))) > 0) {
    colnames(Y) <- paste0("cell.", c(1:ncol(Y))) }
  if (is.null(rownames(Y))) {
    rownames(Y) <- paste0("gene.", c(1:nrow(Y))) }

  #<--------------------------------------
  # Make "DGEList" object
  dge <- edgeR::DGEList(counts = Y,
                        group = X[,2],
                        genes = rownames(Y))

  dge <- edgeR::calcNormFactors(dge)

  # estimate dispersion
  dge <- edgeR::estimateDisp(dge, design = X)

  # Run DE analysis; dispersion = NULL will extract tagwise (genewise) dispersion estimates
  # for DE analysis
  fit <- edgeR::glmFit(dge, dispersion = NULL)

  # Run LRT test
  lrt <- edgeR::glmLRT(fit, coef = 2)

  betahat <- lrt$coefficients[,2]
  pvalue <- lrt$table$PValue

  res <- list(log2FoldChange=betahat,
              pvalue=pvalue)

  return(list(pval = res$pvalue,
              est = res$log2FoldChange))
}

