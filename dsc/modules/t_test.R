t_test <- function(Y, X)
{
  groupInd <- X[,2]
  Y1 <- Y[,groupInd==1]
  Y2 <- Y[,groupInd==0]

  res <-  sapply(1:nrow(Y1), function(i) {
    t <- try(t.test(Y1[i,],Y2[i,]))
    if(class(t) == "try-error") {
      return(c(NA, NA))
    } else {
      c(t$estimate[1]-t$estimate[2],t$p.value) }})
  return(res)
}
