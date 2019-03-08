run_mast_ed <- function(Y, X,
                        pseudocount=1) {

  library(MAST)
  #  Y <- cbind(as.matrix(Y1), as.matrix(Y2))

  #  n1 <- dim(Y1)[2]
  #  n2 <- dim(Y2)[2]

  # x <- rep(c(1,2), times = c(n1, n2))
  # x <- factor(x)
  x <- X[,2]
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

