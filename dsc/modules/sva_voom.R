# function directly copied from https://github.com/dsgerard/mouthwash_sims/Code/non_nc_methods.R

# svaseq input is count matrix,
# within the function, count matrix is transformed to log(count+constant), where the constant defaults 1
#
sva_voom <- function(Y, X, num_sv=NULL) {

  if (is.null(num_sv)) {
    num_sv <- sva::num.sv(dat = Y, mod = X)
    cat("num_sv = ", num_sv, "\n")
  }
#  trash      <- capture.output(sva_out <- sva::sva(dat = Y, mod = X, n.sv = num_sv))
  try_sva      <- tryCatch(sva_out <- sva::svaseq(dat = Y, mod = X, n.sv = num_sv), error = function(e) "error")
  if (try_sva=="error")  {
    pvalues <- rep(NA, nrow(Y))
    betahat <- rep(NA, nrow(Y))
    sebetahat <- rep(NA, nrow(Y))
    df <- rep(NA, nrow(Y))
  } else {
    X.sv       <- cbind(X, sva_out$sv)
    dge <- edgeR::calcNormFactors(dge, method="none")

    voom_out   <- limma::voom(counts = Y, design = X.sv)
    limma_out  <- limma::lmFit(object = voom_out)
    ebayes_out <- limma::ebayes(fit = limma_out)
    betahat    <- limma_out$coefficients[, 2]
    sebetahat  <- sqrt(ebayes_out$s2.post) * limma_out$stdev.unscaled[, 2]
    df         <- ebayes_out$df.total[1]
    pvalues    <- ebayes_out$p.value[, 2]
  }
  return(list(betahat = betahat, sebetahat = sebetahat, df = df,
              pvalues = pvalues))
}
