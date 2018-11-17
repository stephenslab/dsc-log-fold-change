
run_deseq2 <- function(Y1, Y2) {
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

  ncores_default <- detectCores()
  if (ncores_default == 8) {
    register(MulticoreParam(4))
  } else {
    register(MulticoreParam(ncores_default))
  }
  dds <- DESeq(dds, parallel = TRUE)
  res <- results(dds, contrast = c("condition", levels(factor(x))[1],
                                   levels(factor(x))[2]), alpha = 0.05);

  return(list(pval = res$pvalue,
              est = res$log2FoldChange,
              se = res$lfcSE))
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
