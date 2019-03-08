run_glm <- function(Y, X, family) {
  #  Y <- cbind(Y1, Y2)

  #  x <- rep(c(0, 1), c(ncol(Y1), ncol(Y2)))
  x <- X[,2]

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
