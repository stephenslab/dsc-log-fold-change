limma_voom <- function(Y, X, normalize_counts=F){

  library(limma)
  library(edgeR)
  #  Y <- as.matrix(cbind(Y1, Y2))
  #  condition <- c(rep(1, ncol(Y1)), rep(2, ncol(Y2)))

  if (normalize_counts) {
    design <- X
    Y <- DGEList(Y)
    dge <- edgeR::calcNormFactors(Y)
    v <- voom(dge,design,plot=FALSE)
    fit <- lmFit(v,design)
    fit.ebayes <- eBayes(fit)
  }

  # input is log expression, already normalized
  if (normalize_counts==F) {
    design <- X
#    Y <- DGEList(Y)
#    dge <- edgeR::calcNormFactors(Y)
#    v <- voom(dge,design,plot=FALSE)
    fit <- lmFit(Y,design, weights=NULL)
    fit.ebayes <- eBayes(fit)
  }

  # given that the condition is a binary vector
  # extract the coefficient corresponds to the difference between the two conditions
  betahat <- fit.ebayes$coefficients[,2]
  sebetahat <- with(fit.ebayes, stdev.unscaled[,2]*sigma)
  pvalue <- fit.ebayes$p.value[,2]
  df <- fit.ebayes$df.total

  return(list(betahat=betahat, sebetahat=sebetahat,
              df=df, pvalue = pvalue))
}
