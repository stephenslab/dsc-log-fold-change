# This document contains codes for method performance evaluation

#' @title Compute sensitivity or true positive rate given fixed FDR
#'
#' @param true_positive_rate numeric vector of sensitivity.
#' @param false_positive_rate numeric vector of 1 - specificity.
#' @param fdr_cutoff default .05
#'
#' @return
#'    \code{tpr} true postive rate at the given FDR threshold.
#' @export
getTPR <- function(true_positive_rate, false_positive_rate,
                   fdr_cutoff = .05) {

  df <- data.frame(true_positive_rate = true_positive_rate,
                   false_positive_rate = false_positive_rate)
  df <- df[order(df$false_positive_rate, df$true_positive_rate), ]
  which_max_fdr <- max(which(df$false_positive_rate < .05))
  tpr <- with(df, true_positive_rate[which_max_fdr])
  return(tpr)
}
