#!/usr/bin/env dsc

# pipeline variables  --------------------------------------------------
# $Y: 'ngene' by `nsamp` matrix of counts for samples
# $X: model design matrix of 'nsamp' rows and 2 columns: the second column encodes to the group effect (0,1), and the first column corresponds to the intercept
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
#   output:  $Y, $X, $beta
# method:
#   input: $Y, $X
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
    # generate simulated data
    data: data_poisthin, data_poisthin_null

    # differential expression analysis
    method: edger, deseq2, limma_voom, t_test_log2cpm_quant, wilcoxon

    # scoring methods
    pval_rank: qvalue
    score: type_one_error, fdr, auc
  run:
    pipe_power: data_poisthin * method * pval_rank
    pipe_type1: data_poisthin_null * method * pval_rank

  exec_path: modules

#  output:
#    /scratch/midway2/joycehsiao/dsc-log-fold-change/benchmark


# data modules ----------------------------------------------------------

data_poisthin: R(counts = readRDS(dataFile)) + \
       poisthin.R + \
       R(set.seed(seed=seed); out = poisthin(mat=t(counts), nsamp=n1+n2, ngene=ngene, gselect=gselect, shuffle_sample=shuffle_sample, signal_fun=signal_fun, signal_params=list(betapi=1, betamu=0, betasd=betasd), prop_null = prop_null)) + \
       R(X <- out$X; Y <- t(out$Y); beta <- out$beta)
  dataFile: "data/pbmc_counts.rds"
  seed: R{2:11}
  n1, n2: (45, 45), (250, 250)
  ngene: 1000
  prop_null: .9
  shuffle_sample: T, F
  gselect: "random"
  signal_fun: "bignormal"
  betasd: .5, 1, 2, 4
  $Y: Y
  $X: X
  $beta: beta



data_poisthin_null(data_poisthin): R(counts = readRDS(dataFile)) + \
       poisthin.R + \
       R(set.seed(seed=seed); out = poisthin(mat=t(counts), nsamp=n1+n2, ngene=ngene, gselect=gselect, shuffle_sample=shuffle_sample, signal_fun=signal_fun, signal_params=list(betapi=1, betamu=0, betasd=betasd), prop_null = prop_null)) + \
       R(X <- out$X; Y <- t(out$Y); beta <- out$beta)
  seed: R{2:11}
  n1, n2: (45, 45), (250, 250)
  ngene: 1000, 10000
  prop_null: 1
  shuffle_sample: T, F
  gselect: "random"
  signal_fun: "bignormal"
  betasd: 1,



# method modules ------------------------------------------------------------------

deseq2: deseq2.R + \
      R(res <- deseq2(Y, X))
  Y: $Y
  X: $X
  $log_fold_change_est: res$est
  $s_hat: res$se
  $pval: res$pval

edger: edger.R + \
      R(res <- edger(Y, X))
  Y: $Y
  X: $X
  $log_fold_change_est: res$est
  $pval: res$pval

#genes are rows
#input is raw counts
glm_pois: glm.R + \
          R(res <- glm(Y, X, family);
            est <- res$Estimate;
            se <- res$Std..Error;
            p <- res$Pr)
  Y: $Y
  X: $X
  family: "poisson"
  $log_fold_change_est: est
  $s_hat: se
  $pval: p

glm_quasipois: glm.R + \
          R(res <- glm(Y, X, family);
            est <- res$Estimate;
            se <- res$Std..Error;
            p <- res$Pr)
  Y: $Y
  X: $X
  family: "quasipoisson"
  $log_fold_change_est: est
  $s_hat: se
  $pval: p


limma_voom: limma_voom.R + \
       R(res <- limma_voom(Y, X))
   Y: $Y
   X: $X
   $pval: res$pvalue
   $log_fold_change_est: res$betahat
   $s_hat: res$sebetahat
   $df: res$df

mast: mast.R + \
       R(res <- mast(Y, X))
   Y: $Y
   X: $X
   $pval: res$pval
   $log_fold_change_est: res$betahat
   $s_hat: res$sebetahat
   $df: res$df

sva: sva_voom.R + \
    R(sva(Y, X, num_sv = num_sv))
    Y: $Y
    X: $X

t_test_log2cpm_quant: t_test.R + normalize_counts.R + \
     R(log2cpm_qqnormed <- normalize_log2cpm(Y)) + \
     R(res <- t_test(log2cpm_qqnormed, X))
   Y: $Y
   X: $X
   $pval: res[2,]
   $log_fold_change_est: res[1,]

wilcoxon: wilcoxon.R + \
        R(res <- wilcoxon(Y, X))
   Y: $Y
   X: $X
   $pval: res[2,]
   $log_fold_change_est: res[1,]


#zinbwave_deseq2: zinbwave_deseq2.R + \
#      R(res <- zinbwave_deseq2(Y1,Y2))
#   Y1: $Y1
#   Y2: $Y2
#   $pval: res$pval
#   $log_fold_change_est: res$betahat


#zinbwave_edger: zinbwave_edger.R + \
#      R(res <- zinbwave_edger(Y1,Y2))
#   Y1: $Y1
#   Y2: $Y2
#   $pval: res$pval
#   $log_fold_change_est: res$betahat





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

#try(qvalue(p=a$pval))

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






