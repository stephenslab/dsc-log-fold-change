tpr: scores.R + \
      R(out <- getTPR(response, predictor = args$, fdr_cutoff=.05))
  Y1: $Y1
  Y2: $Y2
  $log_fold_change_est: res$est
  $p: res$pval
