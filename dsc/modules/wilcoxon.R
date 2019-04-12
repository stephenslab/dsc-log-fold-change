wilcoxon <- function(Y, X) {
  groupInd <- X[,2]
  Y1 <- Y[,groupInd==1]
  Y2 <- Y[,groupInd==0]

  res <- sapply(1:nrow(Y1), function(i){
    w <- try(wilcox.test(Y1[i,],Y2[i,], conf.int = TRUE))
    if(class(w) == "try-error") {
      return(c(NA, NA))
    } else {
      c(w$estimate, w$p.value)}  })
  return(res)
}
