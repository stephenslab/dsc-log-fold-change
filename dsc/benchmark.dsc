#!/usr/bin/env dsc

# pipeline variables  --------------------------------------------------
# $Y: count matrix; 'ngene' by 'nsamp' matrix
# $X: model design matrix; samples in rows and variables in columns
# $beta: simulated effect sizes, i.e., log2 fold-change; a length 'ngene' vector
# $log_fold_change_est: estimated log2 fold-changes (i.e., estimated betas); a length `ngene` vector
# $s_hat: estimated standard errors; a length 'ngene' vector
# $df: degrees of freedom; a 'ngene' vector
# $pval: p-values; a length 'ngene' vector
# $type1error: type I error in one dataset at a given alpha level
# $qval: adjusted p-values, computed using 'qvalue' package; a 'ngene' vector
# $fdr: false discovery rate in one dataset at a given qvalue threshold ('fdr_thres'); a length 'ngene' vector
# $auc: area under the curve of sensitivty vs false discovery rate, computed using 'pROC' package; a length 'ngene' vector


# Outline module groups -------------------------------------------------
# make_data:
#   input: "data/pbmc_counts.rds"
#   output:  $Y, $X, $beta
# run_methods:
#   input: $Y, $X
#   output:  $log_fold_change_est, $s_hat, $pval, $df
# rank_pvals:
#   input: $pval
#   output: $qval
# score_methods:
#   input: $pval, $qval, $beta
#   output: $type1error, $fdr, $auc


# Define DSC modules and pipelines --------------------------------------

DSC:
  define:
    # generate simulated data
    make_data:
      data_poisthin_choose_betasd, data_poisthin, data_poisthin_libsize, data_poisthin_gtex

    # differential expression analysis
    run_methods: edger, deseq2, limma_voom, t_test_log2p1, t_test_log2p1_cpm_quant, wilcoxon

    # scoring methods
    rank_pvals: qvalue
    score_methods: type1error, fdr, auc

  run:
    pipe_choose_betasd: data_poisthin_choose_betasd * (edger, limma_voom) * rank_pvals
    pipe_pbmc: data_poisthin * run_methods * rank_pvals
    pipe_gtex: data_poisthin_gtex * run_methods * rank_pvals
    pipe_null_libsize: data_poisthin_libsize * run_methods * rank_pvals

  exec_path: modules

#  output:
#    /scratch/midway2/joycehsiao/dsc-log-fold-change


# data modules ----------------------------------------------------------

# null data with different library sizes
data_poisthin_libsize(data_poisthin_choose_betasd): R(counts = readRDS(dataFile)) + \
       poisthin.R + filter_genes.R + \
       R(out = poisthin(mat=t(counts), nsamp=n1+n2, ngene=ngene, gselect=gselect, signal_fun=function(n) rep(libsize_factor, n), signal_params=list(), prop_null = prop_null)) + \
       R(X <- out$X; Y <- t(out$Y); beta <- out$beta) + \
       R(keep_genes <- filter_genes(Y, min_cell_detected=1)) + \
       R(Y <- Y[keep_genes,]; beta <- beta[keep_genes])
  n1, n2: (50, 50)
  ngene: 1000
  prop_null: 0
  gselect: "random"
  libsize_factor: 0, 1, 2, 4
  $Y: Y
  $X: X
  $beta: beta


# Description:
#   90% null genes, equal sample sizes across conditions
# for evaluating power in relation to sample size
data_poisthin_choose_betasd: R(counts = readRDS(dataFile)) + \
       poisthin.R + \
       R(out = poisthin(mat=t(counts), nsamp=n1+n2, ngene=ngene, gselect=gselect, signal_params=list(mean=0, sd=betasd), prop_null = prop_null)) + \
       R(X <- out$X; Y <- t(out$Y); beta <- out$beta)
  dataFile: "data/pbmc_counts.rds"
  n1, n2: (50, 50), (100, 100), (250, 250)
  ngene: 1000
  prop_null: .9
  gselect: "random"
  betasd: .5, 1, 2, 4
  $Y: Y
  $X: X
  $beta: beta


# Description:
#   90% null genes, scRNA-seq
#   equal library sizes and sample sizes (across conditions)
data_poisthin: R(counts = readRDS(dataFile)) + \
       poisthin.R + \
       R(out = poisthin(mat=t(counts), nsamp=n1+n2, ngene=ngene, gselect=gselect, signal_params=list(mean=0, sd=betasd), prop_null = prop_null)) + \
       R(X <- out$X; Y <- t(out$Y); beta <- out$beta)
  dataFile: "data/pbmc_counts.rds"
  n1, n2: (50, 50), (250, 250)
  ngene: 1000
  prop_null: .9, 1
  gselect: "random"
  betasd: 1
  $Y: Y
  $X: X
  $beta: beta



# Description:
#   90% null genes, bulk RNA-seq
#   equal library sizes and sample sizes (across conditions)
data_poisthin_gtex: R(counts = readRDS(dataFile)) + \
       poisthin.R + \
       R(out = poisthin(mat=t(counts), nsamp=n1+n2, ngene=ngene, gselect=gselect, signal_params=list(mean=0, sd=betasd), prop_null = prop_null)) + \
       R(X <- out$X; Y <- t(out$Y); beta <- out$beta)
  dataFile: "data/gtex_lung.rds"
  n1, n2: (5, 5), (10, 10), (50, 50), (150, 150)
  ngene: 1000
  prop_null: .9, 1
  gselect: "random"
  betasd: 1
  $Y: Y
  $X: X
  $beta: beta







# Method modules ------------------------------------------------------------------
# For example, one can define method module groups as follows,
#    run_methods_1: edger, deseq2, limma_voom, t_test_log2p1, t_test_log2p1cpm_quant, wilcoxon
#    run_methods_2: edger, deseq2, limma_voom, t_test_log2p1, wilcoxon


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

sva_limma_voom: sva.R + limma_voom.R + \
    R(output_sva <- sva(Y, X)) +\
    R(res <- limma_voom(Y, X=output_sva$X.sv))
    Y: $Y
    X: $X
    $pval: res$pvalues
    $log_fold_change_est: res$betahat
    $s_hat: res$sebetahat
    $df: res$df

#sva_ttest: t_test.R + sva.R +\
#    R(output_sva <- sva(Y, X)) +\
#    R(log2_counts <- log2(Y+1)) + \
#    R(res <- t_test(Y=log2_counts, X=cbind(X,output_sva$X.sv)))
#    Y: $Y
#    X: $X
#   $pval: res[2,]
#   $log_fold_change_est: res[1,]

t_test_log2p1: t_test.R + \
     R(log2_counts <- log2(Y+1)) + \
     R(res <- t_test(log2_counts, X))
   Y: $Y
   X: $X
   $pval: res[2,]
   $log_fold_change_est: res[1,]

t_test_log2p1_cpm_quant: t_test.R + normalize_counts.R + \
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
# rank_pvals: qvalue
# score_methods: type1error, qvalue, fdr, auc

qvalue: qvalue.R + \
  R(qval <- get_qvalue(pval) )
  pval: $pval
  $qval: qval

type1error: R(truth_vec <- beta !=0; pval_null <- pval[which(truth_vec ==F)]) +\
              R(out <- mean(pval_null < pval_thres, na.rm = TRUE) )
  pval_thres: .01
  pval: $pval
  beta: $beta
  $type_one_error: out

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






