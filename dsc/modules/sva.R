# function directly copied from https://github.com/dsgerard/mouthwash_sims/Code/non_nc_methods.R
sva <- function(Y, X, num_sv=NULL) {
  trash     <- capture.output(sva_out <- sva::svaseq(dat = Y, mod = X, n.sv = NULL))
  X.sv      <- cbind(X, sva_out$sv)
  return(list(X.sv = X.sv))
}
