#!/usr/bin/env dsc

# pipeline variables  --------------------------------------------------
# $Y1: 'ngene' by `nsamp/2` matrix of counts for samples in group 1
# $Y2: 'ngene' by `nsamp/2` matrix of counts for samples in group 0
# $beta: an 'ngene' vector of simulated true values beta (used 'poisthin' function)
# $log_fold_change_est: an `ngene` vector of estimated values beta
# $s_hat: an 'ngene' vector of estimated values standard error
# $pval: an 'ngene' vector of p-values
# $df: an 'ngene' vector of degrees of freedom
# $type_one_error: an 'ngene' vector of degrees of freedom
# $qval: an 'ngene' vector of adjusted p-values, currently we use 'qvalue' from the 'qvalue'
# $fdr_est: an 'ngene' vector of estimated vaules for false discover rate (depend on 'fdr_thres' level)
# $auc_est: an 'ngene' vector of estimated values for area under the curve (using pROC package)


# module groups --------------------------------------------------------
# data:
#   input: "data/pbmc_counts.rds"
#   output:  $Y1, $Y2, $beta
# method:
#   input: $Y1, $Y2
#   output:  $log_fold_change_est, $s_hat, $pval, $df
# pval_rank:
#   input: $pval
#   output: $qval
# score:
#   input: $pval, $qval, $beta
#   output: $type_one_error, $fdr_est, $auc_est


# Define DSC modules and pipelines --------------------------------------

DSC:
  define:
    data: data_poisthin
    method: edger, deseq2, glm_pois, glm_quasipois, limma_voom, mast, t_test, wilcoxon
    pval_rank: qvalue
    score: type_one_error, fdr, auc
  run:
    data * method * pval_rank
  exec_path: modules
#  output:
#    /scratch/midway2/joycehsiao/dsc-log-fold-change/benchmark


# simulate modules ----------------------------------------------------------

data_poisthin: R(counts = readRDS(dataFile)) + \
       dataSimulate.R + \
       R(set.seed(seed=seed); out = poisthin(mat=t(counts), nsamp=nsamp, ngene=ngene, gselect=gselect, shuffle_sample=shuffle_sample, signal_dist=signal_dist, prop_null = prop_null)) + \
       R(groupInd = out$X[,2]; Y1 = t(out$Y[groupInd==1,]); Y2 = t(out$Y[groupInd==0,]))
  dataFile: "data/pbmc_counts.rds"
  seed: R{2:101}
  nsamp: 90
  ngene: 1000
  prop_null: .5, .9, 1
  shuffle_sample: T, F
  gselect: "random"
  signal_dist: "bignormal"
  $Y1: Y1
  $Y2: Y2
  $beta: out$beta





# method modules ------------------------------------------------------------------

deseq2: methodsMeanExpression.R + \
      R(res <- run_deseq2(Y1, Y2))
  Y1: $Y1
  Y2: $Y2
  $log_fold_change_est: res$est
  $s_hat: res$se
  $pval: res$pval

edger: methodsMeanExpression.R + \
      R(res <- run_edger(Y1, Y2))
  Y1: $Y1
  Y2: $Y2
  $log_fold_change_est: res$est
  $pval: res$pval

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
  $pval: p

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
  $pval: p


limma_voom: methodsMeanExpression.R + \
       R(res <- run_limma_voom(Y1, Y2))
   Y1: $Y1
   Y2: $Y2
   $pval: res$pvalue
   $log_fold_change_est: res$betahat
   $s_hat: res$sebetahat
   $df: res$df

mast: methodsMeanExpression.R + \
       R(res <- run_mast(Y1, Y2))
   Y1: $Y1
   Y2: $Y2
   $pval: res$pval
   $log_fold_change_est: res$betahat
   $s_hat: res$sebetahat
   $df: res$df

t_test: methodsMeanExpression.R + \
       R(res <- run_t_test(Y1, Y2))
   Y1: $Y1
   Y2: $Y2
   $pval: res[2,]
   $log_fold_change_est: res[1,]


wilcoxon: methodsMeanExpression.R + \
        R(res <- run_wilcoxon(Y1, Y2))
   Y1: $Y1
   Y2: $Y2
   $pval: res[2,]
   $log_fold_change_est: res[1,]


zinbwave_deseq: methodsMeanExpression.R + \
      R(res <- run_zinbwave_deseq2(Y1,Y2))
   Y1: $Y1
   Y2: $Y2
   $pval: res$pval
   $log_fold_change_est: res$betahat


zinbwave_edger: methodsMeanExpression.R + \
      R(res <- run_zinbwave_edger(Y1,Y2))
   Y1: $Y1
   Y2: $Y2
   $pval: res$pval
   $log_fold_change_est: res$betahat


# Scoring modules --------------------------------------------------------------

type_one_error: R(truth_vec <- beta !=0; pval_null <- pval[which(truth_vec ==F)]) +\
              R(out <- mean(pval_null < pval_thres, na.rm = TRUE) )
  pval_thres: .01
  pval: $pval
  beta: $beta
  $type_one_error: out

qvalue: R(library(qvalue); qval <- qvalue(p=pval)$qvalues)
  pval: $pval
  $qval: qval

fdr: R(truth_vec <- beta !=0) + \
    R(fdr_est <- sum(qval < fdr_thres & !truth_vec, na.rm=TRUE)/sum(qval < fdr_thres, na.rm=TRUE))
  fdr_thres: .05
  beta: $beta
  qval: $qval
  $fdr_est: fdr_est

auc: R(truth_vec <- beta !=0) + \
      R(library(pROC); if(sum(truth_vec==F) == length(truth_vec)) {auc_est <- NA} else {auc_est <- roc(response=truth_vec, predictor=qval)$auc})
    beta: $beta
    qval: $qval
    $auc_est: auc_est






