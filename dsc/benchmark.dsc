#!/usr/bin/env dsc

# pipeline variables  --------------------------------------------------
# $Y1: `ngene` by `nsamp/2` matrix of counts for samples in group 1
# $Y2: `ngene` by `nsamp/2` matrix of counts for samples in group 0
# $beta: an `ngene` vector of simulated true values beta (used `poisthin` function)
# $log_fold_change_est: an `ngene` vector of estimated values beta
# $s_hat: an `ngene` vector of estimated values standard error
# $pval: an `ngene` vector of p-values
# $df: an `ngene` vector of degrees of freedom
# $type_one_error: an 'ngene' vector of degrees of freedom
# $pval_adj: an 'ngene' vector of adjusted p-values, currently we use 'qvalue' from the 'qvalue'
# $fdr_est: an 'ngene' vector of estimated vaules for false discover rate (depend on 'fdr_thres' level)
# $auc_est: an 'ngene' vector of estimated values for area under the curve (using pROC package)


# module groups --------------------------------------------------------
# data:
#   input: "data/pbmc_counts.rds"
#   output:  $Y1, $Y2, $beta
# method:
#   input: $Y1, $Y2
#   output:  $log_fold_change_est, $s_hat, $pval, $df
# score:
#   input: $pval, $pval_adj, $beta
#   output: $type_one_error, $pval_adj, $fdr_est, $auc_est


# Define DSC modules and pipelines --------------------------------------

DSC:
  define:
    data: data_poisthin_null, data_poisthin_signal
    method: edger, deseq2, glm_pois, glm_quasipois, limma_voom, mast, t_test, wilcoxon
    score: type_one_error, pval_adj, fdr, auc
  run:
    pipe_typeone: data_poisthin_null * method * type_one_error
    pipe_power: data_poisthin_signal * method * pval_adj * (fdr, auc)
  exec_path: modules
  global:
    dataFile: "data/pbmc_counts.rds"
#  output:
#    /scratch/midway2/joycehsiao/dsc-log-fold-change/benchmark


# simulate modules ----------------------------------------------------------

data_poisthin_signal: R(counts = readRDS(dataFile)) + \
       dataSimulate.R + \
       R(out = poisthin(mat=t(counts), nsamp=nsamp, ngene=ngene, gselect=gselect, signal_dist=signal_dist, prop_null = prop_null)) + \
       R(groupInd = out$X[,2]; Y1 = t(out$Y[groupInd==1,]); Y2 = t(out$Y[groupInd==0,]))
  dataFile: "data/pbmc_counts.rds"
  seed: R{2:101}
  nsamp: 90
  ngene: 1000
  prop_null: .9
  gselect: "random"
  signal_dist: "bignormal"
  $Y1: Y1
  $Y2: Y2
  $beta: out$beta


data_poisthin_null (data_poisthin_signal): R(counts = readRDS(dataFile)) + \
     dataSimulate.R + \
     R(out = poisthin(mat=t(counts), nsamp=nsamp, ngene=ngene, gselect=gselect, signal_dist=signal_dist, prop_null = prop_null)) + \
     R(groupInd = out$X[,2]; Y1 = t(out$Y[groupInd==1,]); Y2 = t(out$Y[groupInd==0,]))
  prop_null: 1
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

type_one_error: R(out <- mean(pval < pval_thres, na.rm = TRUE) )
  pval_thres: .05
  pval: $pval
  $type_one_error: out

pval_adj: R(library(qvalue); pval_adj <- qvalue(p=pval)$qvalues)
  pval: $pval
  $pval_adj: pval_adj

fdr: R(truth_vec <- beta !=0) + \
    R(fdr_est <- sum(pval_adj < fdr_thres & !truth_vec, na.rm=TRUE)/sum(pval_adj < fdr_thres, na.rm=TRUE))
  fdr_thres: .05
  beta: $beta
  pval_adj: $pval_adj
  $fdr_est: fdr_est

auc: R(truth_vec <- beta !=0) + \
      R(library(pROC); auc_est <- roc(response=truth_vec, predictor=pval_adj)$auc)
    beta: $beta
    pval_adj: $pval_adj
    $auc_est: auc_est






