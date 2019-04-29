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
    data: data_poisthin

#    # normalization
#    count_normalize: log2_cpm

    # differential expression analysis
    method: edger, deseq2, limma_voom, t_test_log2cpm, t_test_log2cpm_quant, wilcoxon

    # estimate latent confounding variables
#    method_latent_factor: sva

    # scoring methods
    pval_rank: qvalue
    score: type_one_error, fdr, auc
  run:
    pipe_type1: data * method * pval_rank

  exec_path: modules

#  output:
#    /scratch/midway2/joycehsiao/dsc-log-fold-change/benchmark


# data modules ----------------------------------------------------------

data_poisthin: R(counts = readRDS(dataFile)) + \
       dataSimulate.R + \
       R(set.seed(seed=seed); out = poisthin(mat=t(counts), nsamp=nsamp, ngene=ngene, gselect=gselect, shuffle_sample=shuffle_sample, signal_fun=signal_fun, signal_params=list(betapi=1, betamu=0, betasd=betasd), prop_null = prop_null)) +\
       R(X <- out$X; Y <- t(out$Y))
#       R(groupInd = out$X[,2]; Y1 = t(out$Y[groupInd==1,]); Y2 = t(out$Y[groupInd==0,]))
  dataFile: "data/pbmc_counts.rds"
  seed: R{2:11}
  nsamp: 90
  ngene: 1000, 10000
  prop_null: 1
  shuffle_sample: T, F
  gselect: "random"
  signal_fun: "bignormal"
  betasd: 1
  $Y: Y
  $X: X
  $beta: out$beta


#log2_cpm: R(counts = cbind(Y1, Y2) + \
#      R(libsize=colSums(count); log2cpm=log2(t(10^6*(t(counts)/libsize)+1)))
#  Y1: $Y1
#  Y2: $Y2


# method modules ------------------------------------------------------------------

deseq2: methodsMeanExpression.R + \
      R(res <- run_deseq2_ed(Y, X))
  Y: $Y
  X: $X
  $log_fold_change_est: res$est
  $s_hat: res$se
  $pval: res$pval

edger: methodsMeanExpression.R + \
      R(res <- run_edger_ed(Y, X))
  Y: $Y
  X: $X
  $log_fold_change_est: res$est
  $pval: res$pval

#genes are rows
#input is raw counts
glm_pois: methodsMeanExpression.R + \
          R(res <- run_glm_ed(Y, X, family);
            est <- res$Estimate;
            se <- res$Std..Error;
            p <- res$Pr)
  Y: $Y
  X: $X
  family: "poisson"
  $log_fold_change_est: est
  $s_hat: se
  $pval: p

glm_quasipois: methodsMeanExpression.R + \
          R(res <- run_glm_ed(Y, X, family);
            est <- res$Estimate;
            se <- res$Std..Error;
            p <- res$Pr)
  Y: $Y
  X: $X
  family: "quasipoisson"
  $log_fold_change_est: est
  $s_hat: se
  $pval: p


limma_voom: methodsMeanExpression.R + \
       R(res <- run_limma_voom_ed(Y, X))
   Y: $Y
   X: $X
   $pval: res$pvalue
   $log_fold_change_est: res$betahat
   $s_hat: res$sebetahat
   $df: res$df

mast: methodsMeanExpression.R + \
       R(res <- run_mast_ed(Y, X))
   Y: $Y
   X: $X
   $pval: res$pval
   $log_fold_change_est: res$betahat
   $s_hat: res$sebetahat
   $df: res$df

t_test_counts: methodsMeanExpression.R + \
       R(res <- run_t_test_ed(Y, X))
   Y: $Y
   X: $X
   $pval: res[2,]
   $log_fold_change_est: res[1,]

t_test_counts: methodsMeanExpression.R + \
       R(res <- run_t_test_ed(Y, X))
   Y: $Y
   X: $X
   $pval: res[2,]
   $log_fold_change_est: res[1,]

t_test_log2cpm: methodsMeanExpression.R + \
       R(counts=Y; libsize=colSums(counts); log2cpm=log2(t(10^6*(t(counts)/libsize)+1))) + \
       R(res <- run_t_test_ed(Y=log2cpm, X))
   Y: $Y
   X: $X
   $pval: res[2,]
   $log_fold_change_est: res[1,]

t_test_log2cpm_quant: methodsMeanExpression.R + \
     R(counts=Y; libsize=colSums(counts); log2cpm=log2(t(10^6*(t(counts)/libsize)+1))) + \
     R(log2cpm_qqnormed <- do.call(rbind, lapply(1:nrow(log2cpm), function(g) {qqnorm(log2cpm[g,])$x}))) + \
   R(res <- run_t_test_ed(Y=log2cpm_qqnormed, X))
   Y: $Y
   X: $X
   $pval: res[2,]
   $log_fold_change_est: res[1,]


wilcoxon: methodsMeanExpression.R + \
        R(res <- run_wilcoxon_ed(Y, X))
   Y: $Y
   X: $X
   $pval: res[2,]
   $log_fold_change_est: res[1,]


#zinbwave_deseq: methodsMeanExpression.R + \
#      R(res <- run_zinbwave_deseq2(Y1,Y2))
#   Y1: $Y1
#   Y2: $Y2
#   $pval: res$pval
#   $log_fold_change_est: res$betahat


#zinbwave_edger: methodsMeanExpression.R + \
#      R(res <- run_zinbwave_edger(Y1,Y2))
#   Y1: $Y1
#   Y2: $Y2
#   $pval: res$pval
#   $log_fold_change_est: res$betahat


#sva: methodsLatentFactor.R + \
#    R(sva(Y, X, method = ))
#    Y: $Y
#    X: $X



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

