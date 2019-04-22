# function directly copied from https://github.com/dsgerard/mouthwash_sims/Code/non_nc_methods.R
sva_voom <- function(Y, X, num_sv=NULL) {

  if (is.null(num_sv)) {
    num_sv <- sva::num.sv(dat = Y, mod = X)
    cat("num_sv = ", num_sv, "\n")
  }
  trash      <- capture.output(sva_out <- sva::sva(dat = Y, mod = X, n.sv = num_sv))
  X.sv       <- cbind(X, sva_out$sv)
  voom_out   <- limma::voom(counts = Y, design = X.sv)
  limma_out  <- limma::lmFit(object = voom_out)
  ebayes_out <- limma::ebayes(fit = limma_out)
  betahat    <- limma_out$coefficients[, 2]
  sebetahat  <- sqrt(ebayes_out$s2.post) * limma_out$stdev.unscaled[, 2]
  df         <- ebayes_out$df.total[1]
  pvalues    <- ebayes_out$p.value[, 2]
  return(list(betahat = betahat, sebetahat = sebetahat, df = df,
              pvalues = pvalues))
}
