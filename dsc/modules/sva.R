#' Estimate number of latent factors
#'
#' Use SVA package to estimate the number of latent factors.
#' This number is usually required for differential expression analysis
#' that accounts for latent factors. Most methods will estimate this
#' number for the user if not provided, in which case the computation will
#' take more time.

run_sva <- function(Y_log_normed, X) {
  library(sva)
  num_sv <- num.sv(Y_log_normed,X,method="leek")
  X_0 <- model.matrix(~1)
  res_sv <- sva(Y_log_normed, X, X_0, n.sv=num_sv)
  return(res_sv)
}

