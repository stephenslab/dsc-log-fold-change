# estimate number of latent factors ----------

nf_sva: methodsLatenFactor.R + \
      R(res <- nf_sva(Y1, Y2))
  Y1: $Y1
  Y2: $Y2
  $log_fold_change_est: res$est
  $s_hat: res$se
  $p: res$pval


# estimate latent variables ----------

factor_est: methodsLatenFactor.R + \
      R(res <- run_sva(Y, X))
  Y1: $Y1
  Y2: $Y2
  $log_fold_change_est: res$est
  $s_hat: res$se
  $p: res$pval

