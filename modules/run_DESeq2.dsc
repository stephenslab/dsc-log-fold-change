DESeq2: run_DESeq2.R + \
          R(res <- run_DESeq2(Y1, Y2))
  @CONF: R_libs = (DESeq2, BiocParallel)
  Y1: $Y1
  Y2: $Y2
  $log_fold_change_est: res$est
  $s_hat: res$se
  $p: res$pval
