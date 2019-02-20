#!/usr/bin/env dsc


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


# Simulate data ----------------------------------------------------------

data_poisthin_signal: R(counts = readRDS(${dataFile})) + \
       dataSimulate.R + \
       R(out = poisthin(mat=t(counts), nsamp = args$nsamp, ngene = args$ngene, gselect = args$gselect, signal_dist = "bignormal", prop_null = args$prop_null)) + \
       R(groupInd = out$X[,2]; Y1 = t(out$Y[groupInd==1,]); Y2 = t(out$Y[groupInd==0,]))
  seed: R{2:101}
  nsamp: 90
  ngene: 1000
  prop_null: .9
  gselect: "random"
  signal_dist: "bignormal"
  @ALIAS: args = List(!dataFile)
  $Y1: Y1
  $Y2: Y2
  $beta: out$beta


data_poisthin_null: R(counts = readRDS(${dataFile})) + \
     dataSimulate.R + \
     R(out = poisthin(mat=t(counts), nsamp = args$nsamp, ngene = args$ngene, gselect = args$gselect, signal_dist = args$signal_dist, prop_null = args$prop_null)) + \
     R(groupInd = out$X[,2]; Y1 = t(out$Y[groupInd==1,]); Y2 = t(out$Y[groupInd==0,]))
  seed: R{2:101}
  nsamp: 90
  ngene: 1000
  prop_null: 1
  gselect: "random"
  signal_dist: "big_normal"
  @ALIAS: args = List(!dataFile)
  $Y1: Y1
  $Y2: Y2
  $beta: out$beta


# DE methods -------------------------------------------------------------

deseq2: methodsMeanExpression.R + \
      R(res <- run_deseq2(Y1, Y2))
  Y1: $Y1
  Y2: $Y2
  $log_fold_change_est: res$est
  $s_hat: res$se
  $p: res$pval

edger: methodsMeanExpression.R + \
      R(res <- run_edger(Y1, Y2))
  Y1: $Y1
  Y2: $Y2
  $log_fold_change_est: res$est
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
   Y1: $Y1
   Y2: $Y2
   $p: res$pvalue
   $log_fold_change_est: res$betahat
   $s_hat: res$sebetahat
   $df: res$df

mast: methodsMeanExpression.R + \
       R(res <- run_mast(Y1, Y2))
   Y1: $Y1
   Y2: $Y2
   $p: res$pval
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


zinbwave_deseq: methodsMeanExpression.R + \
      R(res <- run_zinbwave_deseq2(Y1,Y2))
   Y1: $Y1
   Y2: $Y2
   $p: res$pval
   $log_fold_change_est: res$betahat


zinbwave_edger: methodsMeanExpression.R + \
      R(res <- run_zinbwave_edger(Y1,Y2))
   Y1: $Y1
   Y2: $Y2
   $p: res$pval
   $log_fold_change_est: res$betahat


# Score methods --------------------------------------------------------------

type_one_error: R(out <- mean(pval < .05, na.rm = TRUE) )
  pval: $p
  $type_one_error: out

pval_adj: R(library(qvalue); pval_adj <- qvalue(p=pval)$qvalues)
  pval: $p
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






