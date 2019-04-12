# function directly copied from https://github.com/dsgerard/mouthwash_sims/Code/non_nc_methods.R
sva <- function(Y, X, num_sv) {
  trash     <- capture.output(sva_out <- sva::sva(dat = t(Y), mod = X, n.sv = num_sv))
  X.sv      <- cbind(X, sva_out$sv)
  limma_out <- limma::lmFit(object = t(Y), design = X.sv)
  betahat   <- limma_out$coefficients[, 2]
  sebetahat <- limma_out$stdev.unscaled[, 2] * limma_out$sigma
  df        <- limma_out$df.residual[1]
  tstats    <- betahat / sebetahat
  pvalues   <- 2 * pt(-abs(tstats), df = df)
  return(list(betahat = betahat, sebetahat = sebetahat, df = df,
              pvalues = pvalues))
}
