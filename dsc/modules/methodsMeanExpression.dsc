deseq2: methodsMeanExpression.R + \
      R(res <- run_deseq2(Y1, Y2))
  Y1: $Y1
  Y2: $Y2
  $log_fold_change_est: res$est
  $s_hat: res$se
  $p: res$pval


#genes are rows
#input is raw counts
glm_pois: methodsMeanExpression.R + \
          R(res <- run_glm(Y1, Y2, family);
            est <- res$Estimate;
            se <- res$Std..Error;
            p <- res$Pr)
  Y1: $Y1
  Y2: $Y2
  family: "poisson"
  $log_fold_change_est: est
  $s_hat: se
  $p: p

glm_quasipois: methodsMeanExpression.R + \
          R(res <- run_glm(Y1, Y2, family);
            est <- res$Estimate;
            se <- res$Std..Error;
            p <- res$Pr)
  Y1: $Y1
  Y2: $Y2
  family: "quasipoisson"
  $log_fold_change_est: est
  $s_hat: se
  $p: p


limma_voom: methodsMeanExpression.R + \
       R(res <- run_limma_voom(Y1, Y2))
#  @CONF: R_libs = (limma, edgeR)
   Y1: $Y1
   Y2: $Y2
   $p: res$pvalue
   $log_fold_change_est: res$betahat
   $s_hat: res$sebetahat
   $df: res$df


t_test: methodsMeanExpression.R + \
       R(res <- run_t_test(Y1, Y2))
   Y1: $Y1
   Y2: $Y2
   $p: res[2,]
   $log_fold_change_est: res[1,]


wilcoxon: methodsMeanExpression.R + \
        R(res <- run_wilcoxon(Y1, Y2))
   Y1: $Y1
   Y2: $Y2
   $p: res[2,]
   $log_fold_change_est: res[1,]

