# function directly copied from https://github.com/dsgerard/mouthwash_sims/Code/non_nc_methods.R
cate_rr <- function(Y, X, num_sv, calibrate = FALSE) {
  calibrate <- as.logical(calibrate)
  cate_rr <- cate::cate.fit(Y = Y, X.primary = X[, 2, drop = FALSE],
                            X.nuis = X[, -2, drop = FALSE],
                            r = num_sv, adj.method = "rr",
                            calibrate = calibrate, fa.method = "pc")
  betahat     <- c(cate_rr$beta)
  sebetahat   <- c(sqrt(c(cate_rr$beta.cov.row) * c(cate_rr$beta.cov.col)) / sqrt(nrow(X)))
  pvalues     <- c(cate_rr$beta.p.value)
  df          <- Inf

  if (calibrate) {
    lambda <- stats::mad(x = betahat / sebetahat, center = 0)
    sebetahat <- sebetahat * lambda
  }

  return(list(betahat = betahat, sebetahat = sebetahat, df = df,
              pvalues = pvalues))
}
