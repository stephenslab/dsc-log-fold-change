limma_voom <- function(Y, X){

  library(limma)
  library(edgeR)

  # default is TMM normalization
   design <- X
   v <- voom(Y,design,plot=FALSE)
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
