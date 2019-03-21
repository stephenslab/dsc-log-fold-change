# function directly copied from https://github.com/dsgerard/mouthwash_sims/Code/non_nc_methods.R

mouthwash <- function(Y, X, num_sv, likelihood = c("normal", "t"), alpha = 0, scale_var = FALSE,
                      var_inflate_pen = 0, mixing_dist = NULL) {
  likelihood <- match.arg(likelihood)
  if (is.null(mixing_dist) & likelihood == "t") {
    mixing_dist <- "sym_uniform"
  } else if (is.null(mixing_dist)) {
    mixing_dist <- "normal"
  }
  mout <- vicar::mouthwash(Y = Y, X = X, k = num_sv, likelihood = likelihood,
                           scale_var = scale_var, sprop = alpha, mixing_dist = mixing_dist,
                           var_inflate_pen = var_inflate_pen)
  return_list <- list()
  return_list$betahat <- mout$result$PosteriorMean
  return_list$lfdr    <- mout$result$lfdr
  return_list$pi0hat  <- mout$pi0
  return(return_list)
}
