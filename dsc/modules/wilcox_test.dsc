# Two sample Test: Wilcoxon Test, input: the log count, output: the estimated log fold change and the p value

wilcoxon_test: R(res= sapply(1:nrow(x), function(i){ w = try(wilcox.test(x[i,],y[i,], conf.int = TRUE)); if(class(w) == "try-error") { return(c(NA, NA))} else{ c(w$estimate, w$p.value)} }))
   x: $Y1
   y: $Y2
   $p: res[2,]
   $log_fold_change_est: res[1,]
