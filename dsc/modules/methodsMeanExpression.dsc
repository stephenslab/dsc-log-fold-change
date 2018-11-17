deseq2: methodsMeanExpression.R + \
      R(res <- run_DESeq2(Y1, Y2))
  @CONF: R_libs = (DESeq2, BiocParallel)
  Y1: $Y1
  Y2: $Y2
  $log_fold_change_est: res$est
  $s_hat: res$se
  $p: res$pval


#genes are rows
#input is raw counts
glm_pois: methodsMeanExpression.R + \
          R(res <- run_glm(Y1, Y2);
            est <- res$Estimate;
            se <- res$Std..Error;
            p <- res$Pr...z..)
  $log_fold_change_est: est
  $s_hat: se
  $p: p


ttest: methodMeanExpressin.R + \
       R(res <- run_ttest(Y2, Y2))
   x: $Y1
   y: $Y2
   $p: res[2,]
   $log_fold_change_est: res[1,]


wilcoxon: methodMeanExpressin.R + \
        R(res <- run_wilcoxon(Y1, Y2))
   x: $Y1
   y: $Y2
   $p: res[2,]
   $log_fold_change_est: res[1,]

