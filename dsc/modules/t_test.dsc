t_test: R(  res = sapply(1:nrow(x), function(i) { t = try(t.test(x[i,],y[i,])); if(class(t) == "try-error") {return(c(NA, NA))} else {c(t$estimate[1]-t$estimate[2],t$p.value) }}) )
   x: $Y1
   y: $Y2
   $p: res[2,]
   $log_fold_change_est: res[1,]
