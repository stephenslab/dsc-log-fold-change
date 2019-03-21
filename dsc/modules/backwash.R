# function directly copied from https://github.com/dsgerard/mouthwash_sims/Code/non_nc_methods.R
backwash <- function(Y, X, num_sv, alpha = 0, scale_var = FALSE, var_inflate_pen = 0) {
  mout <- vicar::backwash(Y = Y, X = X, k = num_sv, scale_var = scale_var, sprop = alpha,
                          var_inflate_pen = var_inflate_pen)
  return_list <- list()
  return_list$betahat <- mout$result$PosteriorMean
  return_list$lfdr    <- mout$result$lfdr
  return_list$pi0hat  <- mout$pi0
  return(return_list)
}
