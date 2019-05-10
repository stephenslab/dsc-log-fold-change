limma_voom <- function(Y, X, libnorm_factors=NULL){

  library(limma)
  library(edgeR)

  design <- X
  if (is.null(libnorm_factors)) {
    # divided by library size with no adjustment
    v <- voom(Y,design=design,plot=FALSE)
  } else {
    # multiple library size by normalizing factors
    v <- voom(Y,design=design, plot=F, lib.size=colSums(Y)*libnorm_factors)
  }
   fit <- lmFit(v,design)
   fit.ebayes <- eBayes(fit)

  # given that the condition is a binary vector
  # extract the coefficient corresponds to the difference between the two conditions
  betahat <- fit.ebayes$coefficients[,2]
  sebetahat <- with(fit.ebayes, stdev.unscaled[,2]*sigma)
  pvalue <- fit.ebayes$p.value[,2]
  df <- fit.ebayes$df.total

  return(list(betahat=betahat, sebetahat=sebetahat,
              df=df, pvalue = pvalue))
}
