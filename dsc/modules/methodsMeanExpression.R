



run_deseq2 <- function(Y1, Y2) {

  library(DESeq2)
#  library(BiocParallel)

  Y <- cbind(as.matrix(Y1), as.matrix(Y2))

  n1 <- dim(Y1)[2]
  n2 <- dim(Y2)[2]
  x <- rep(c(1,2), times = c(n1, n2))
  x <- factor(x)

  if (sum(duplicated(colnames(Y))) > 0) {
    colnames(Y) <- paste0("cell.", c(1:ncol(Y))) }


  dds <- DESeqDataSetFromMatrix(countData = round(Y),
                                colData = data.frame(condition = x),
                                design = ~condition)

  dds <- estimateSizeFactors(dds,type="poscounts")
  dds <- estimateDispersions(dds, minmu = 1e-3)
  dds <- nbinomLRT(dds, minmu=1e-3, reduced=~1)
  res <- results(dds, name="condition_2_vs_1")

    # res <- results(dds, contrast = c("condition", levels(factor(x))[1],
    #                                levels(factor(x))[2]), alpha = 0.05);
  return(list(pval = res$pvalue,
              est = res$log2FoldChange,
              se = res$lfcSE))
}



run_edger <- function(Y1, Y2) {

  library(edgeR)

  Y <- cbind(as.matrix(Y1), as.matrix(Y2))
  n1 <- dim(Y1)[2]
  n2 <- dim(Y2)[2]
  x <- rep(c(1,2), times = c(n1, n2))
  x <- factor(x)


  if (sum(duplicated(colnames(Y))) > 0) {
    colnames(Y) <- paste0("cell.", c(1:ncol(Y))) }
  if (is.null(rownames(Y))) {
    rownames(Y) <- paste0("gene.", c(1:nrow(Y))) }

  #<--------------------------------------
  # Make "DGEList" object
  dge <- edgeR::DGEList(counts = Y,
                        group = x,
                        genes = rownames(Y))

  dge <- edgeR::calcNormFactors(dge)

  # estimate dispersion
  dge <- edgeR::estimateDisp(dge, design = model.matrix(~x))

  # Run DE analysis; dispersion = NULL will extract tagwise (genewise) dispersion estimates
  # for DE analysis
  fit <- edgeR::glmFit(dge, dispersion = NULL)

  # Run LRT test
  lrt <- edgeR::glmLRT(fit, coef = 2)

  betahat <- lrt$coefficients[,2]
  pvalue <- lrt$table$PValue

  res <- list(log2FoldChange=betahat,
              pvalue=pvalue)

  return(list(pval = res$pvalue,
              est = res$log2FoldChange))
}





run_glm <- function(Y1, Y2, family) {
  Y <- cbind(Y1, Y2)

  x <- rep(c(0, 1), c(ncol(Y1), ncol(Y2)))
  results <- apply(Y, 1, FUN=function(y){
    fit_try <- try(glm(y~x, family=family))
    if (any(class(fit_try) == "try-error")) {
      res <- rep(NA,4)
    } else {
      res <- summary(fit_try)$coefficients[2,]
    }
    return(res[c(1, 2, 4)])
  })
  results <- data.frame(t(results))
  return(results)
}


run_limma_voom <- function(Y1, Y2){

  library(limma)
  Y <- as.matrix(cbind(Y1, Y2))

  condition <- c(rep(1, ncol(Y1)), rep(2, ncol(Y2)))

  design <- model.matrix(~factor(condition))

  dge <- DGEList(Y)
  dge <- edgeR::calcNormFactors(dge)
  v <- voom(dge,design,plot=FALSE)
  fit <- lmFit(v,design)
  fit.ebayes <- eBayes(fit)


    # given that the condition is a binary vector
  # extract the coefficient corresponds to the difference between the two conditions
  betahat <- fit.ebayes$coefficients[,2]
  sebetahat <- with(fit.ebayes, stdev.unscaled[,2]*sigma)
  pvalue <- fit.ebayes$p.value[,2]
  df <- fit.ebayes$df.total

  return(list(betahat=betahat, sebetahat=sebetahat,
              df=df, pvalue = pvalue))
}



run_mast <- function(Y1, Y2,
                     pseudocount=1) {

  library(MAST)
  Y <- cbind(as.matrix(Y1), as.matrix(Y2))

  n1 <- dim(Y1)[2]
  n2 <- dim(Y2)[2]

  x <- rep(c(1,2), times = c(n1, n2))
  x <- factor(x)
  if (sum(duplicated(colnames(Y))) > 0) {
    colnames(Y) <- paste0("cell.", c(1:ncol(Y))) }
  if (is.null(rownames(Y))) {
    rownames(Y) <- paste0("gene.", c(1:nrow(Y))) }
  # compute log2CPM
  log2CPM <- log2(Y + pseudocount)

  # make data.frame into singleCellAssay object
  colData <- data.frame(condition = x, row.names = colnames(Y))
  rowData <- data.frame(gene = rownames(Y))
  sca <- suppressMessages(MAST::FromMatrix(log2CPM, colData, rowData))

  # calculate cellualr detection rate; normalized to mean 0 and sd 1
  colData(sca)$cdr.normed <- scale(colSums(assay(sca) > 0))
  # the default method for fitting is bayesGLM
  fit <- suppressMessages(MAST::zlm(~ condition + cdr.normed,
                                    sca=sca))

  # LRT test for the significance of the condition effect
  lrt <-  suppressMessages(MAST::lrTest(fit, "condition"))

  # extract p.value from their "hurdle" model
  pvalue <- lrt[,3,3]

  # extract effect size, standard error, and df from the non-zero component
  betahat <- fit@coefC[,2]
  setbetahat <- sqrt(sapply(1:dim(fit@vcovC)[3], function(i) {
    diag(fit@vcovC[,,i])[2]} ) )
  df <- fit@df.resid[,1]

  return(list(betahat=betahat,
              sebetahat=setbetahat,
              df=df,
              pval=pvalue))
}

run_t_test <- function(Y1, Y2)
{
    res <-  sapply(1:nrow(Y1), function(i) {
    t <- try(t.test(Y1[i,],Y2[i,]))
    if(class(t) == "try-error") {
      return(c(NA, NA))
      } else {
    c(t$estimate[1]-t$estimate[2],t$p.value) }})
  return(res)
}

run_wilcoxon <- function(Y1, Y2) {

  res <- sapply(1:nrow(Y1), function(i){
    w <- try(wilcox.test(Y1[i,],Y2[i,], conf.int = TRUE))
    if(class(w) == "try-error") {
      return(c(NA, NA))
      } else {
    c(w$estimate, w$p.value)}  })
  return(res)
}


run_zinbwave_deseq2 <- function(Y1, Y2){
  library(DESeq2)
  library(zinbwave)

  Y <- cbind(as.matrix(Y1), as.matrix(Y2))

  n1 <- dim(Y1)[2]
  n2 <- dim(Y2)[2]
  x <- rep(c(1,2), times = c(n1, n2))
  x <- factor(x)

  if (sum(duplicated(colnames(Y))) > 0) {
    colnames(Y) <- paste0("cell.", c(1:ncol(Y))) }
  if (is.null(rownames(Y))) {
    rownames(Y) <- paste0("gene.", c(1:nrow(Y))) }

  condition <- x
  design <- model.matrix(~ condition)

  # compute zinbwave weights
  zinb <- zinbFit(Y, X = design, epsilon = 1e12)
  weights <- computeObservationalWeights(zinb, Y)

  # use DESeq2
  d <- DESeqDataSetFromMatrix(Y,
                              colData= DataFrame(data.frame(condition=condition)),
                              design= ~ condition)
  d <- estimateSizeFactors(d, type="poscounts")
  dimnames(weights) = NULL
  assays(d)[["weights"]] = weights
  d <- estimateDispersions(d, minmu = 1e-3)
  #dse = nbinomWaldTest(dse, betaPrior=TRUE, useT=TRUE, df=rowSums(weights)-2, minmu = 1e-3)
  d <- nbinomLRT(d, minmu=1e-3, reduced=~1)
  res <- results(d, name="condition_2_vs_1")

  return(list(pval = res$pvalue,
              betahat = res$log2FoldChange) )
}


run_zinbwave_edger <- function(Y1, Y2){
  library(edgeR)
  library(zinbwave)

  Y <- cbind(as.matrix(Y1), as.matrix(Y2))

  n1 <- dim(Y1)[2]
  n2 <- dim(Y2)[2]
  x <- rep(c(1,2), times = c(n1, n2))
  x <- factor(x)

  if (sum(duplicated(colnames(Y))) > 0) {
    colnames(Y) <- paste0("cell.", c(1:ncol(Y))) }
  if (is.null(rownames(Y))) {
    rownames(Y) <- paste0("gene.", c(1:nrow(Y))) }

  condition <- x
  design <- model.matrix(~ condition)

  # compute zinbwave weights
  zinb <- zinbFit(Y, X = design, epsilon = 1e12)
  weights <- computeObservationalWeights(zinb, Y)
  # use -edgeR
  d <- DGEList(Y)
  d <- edgeR::calcNormFactors(d)
  d$weights <- weights
  d <- estimateDisp(d, design)
  fit <- glmFit(d,design)
  lrt <- glmWeightedF(fit,coef=2, independentFiltering = TRUE)
  pvalues <- lrt$table$PValue

  return(list(pval = pvalues,
              betahat = lrt$table$logFC) )
}


###############################################################
# Some utility functions

#' @title limma + voom
#'
#' @description Implement limma with voom and apply the empirical Bayes method
#'    in limma. This is a modified version of voom which allows to specify pseudocount and
#'    pseudo library size.
#'    Hence, voom: y <- t(log2(t(counts + pseudocount)/(lib.size + pseudo_libsizes) * 1e+06)).
#'
#'
#' @param counts gene by sample expression count matrix.
#' @param design design matrix, generated by R function model.matrix()
#' @param pseudocount default .5
#' @param pseudo_libsizes default 1.
#'
#' @return
#'    \code{w} Weights of dimension G by N.
#' @author Chiaowen Joyce Hsiao
#'
#' @export
voom.controlPseudocount <- function(counts, design,
                                    pseudocount = .5,
                                    pseudo_libsizes = 1,
                                    span = .5) {
  # this function wad adpated from the limma voom function
  # prior count and library size adjustment are set to be argument
  lib.size <- colSums(counts)

  y <- t(log2(t(counts + pseudocount)/(lib.size + pseudo_libsizes) * 1e+06))
  fit <- lmFit(y, design)
  if (is.null(fit$Amean))  fit$Amean <- rowMeans(y, na.rm = TRUE)

  # compute variance-mean dependency
  sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)
  sy <- sqrt(fit$sigma)
  allzero <- rowSums(counts) == 0
  if (any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  l <- lowess(sx, sy, f = span)
  f <- approxfun(l, rule = 2)
  if (fit$rank < ncol(design)) {
    j <- fit$pivot[1:fit$rank]
    fitted.values <- fit$coef[, j, drop = FALSE] %*% t(fit$design[,
                                                                  j, drop = FALSE])
  }
  else {
    fitted.values <- fit$coef %*% t(fit$design)
  }
  fitted.cpm <- 2^fitted.values
  fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + pseudo_libsizes))
  fitted.logcount <- log2(fitted.count)
  w <- 1/f(fitted.logcount)^4
  dim(w) <- dim(fitted.logcount)

  return(w)
}



