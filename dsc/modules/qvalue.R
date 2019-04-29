library(qvalue)

get_qvalue <- function(pval) {
  try_qval <- tryCatch(qvalue(p=pval)$qvalues, error=function(e) "error")

  if (length(try_qval)==1) {
    ifelse(try_qval=="error", qval <- rep(NA, length(pval)), qval <- try_qval)
  } else {
    qval <- try_qval
  }
  return(qval)
}

