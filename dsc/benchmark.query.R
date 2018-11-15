library(dscrutils)
setwd("~/Dropbox/GitHub/dsc-log-fold-change/")
out <- dscquery("benchmark", c("get_data", "get_data.n1", "get_data.n2", "method", "method.p"))

# random sampling same labels across genes

out.sub <- out[out$get_data=="random_sample" & out$get_data.n1==50,]

res <- vector("list",4)
for (i in 1:nrow(out.sub)) {
  print(i)
  fl <- readRDS(paste0("benchmark/",
                       as.character(out.sub$method.p[i]), ".rds"))
  res[[i]] <- data.frame(method = as.character(out.sub$method)[i],
             n1_n2 = paste0(out.sub$get_data.n1[i],".",
                            out.sub$get_data.n2[i]),
             pval = fl$p,
             stringsAsFactors = F)
}
names(res) <- as.character(out.sub$method)

cols <- c("gray50", "forestgreen", "blue", "orange")
par(mfrow=c(2,2))
for (i in 1:length(res)) {
  hist(res[[i]]$pval, main="",
       xlab = "p-values", ylab = "Frequency",
       nclass = 20, col=cols[i])
  title(main=names(res)[i])
}

qq <- lapply(1:length(res), function(i) {
  qqplot(x=runif(100,0,1), y=res[[i]]$pval, plot.it=F)
})
plot(qq[[1]]$x, qq[[1]]$y, col = "gray50", cex=.7, pch = 16,
     xlab = "Uniform(0,1)", ylab = "Empirical distribution",
     main = "QQ-plot")
points(qq[[2]]$x, qq[[2]]$y, col = "forestgreen", cex=.7, pch = 16)
points(qq[[2]]$x, qq[[2]]$y, col = "blue", cex=.7, pch = 16)
points(qq[[3]]$x, qq[[3]]$y, col = "orange", cex=.7, pch = 16)
abline(0,1, col = "black")
title("random sample; 50 & 50", outer=TRUE, line=-1)




##### random sampling per gene
out.sub <- out[out$get_data=="random_gene" & out$get_data.n1==50,]

res <- vector("list",4)
for (i in 1:nrow(out.sub)) {
  print(i)
  fl <- readRDS(paste0("benchmark/",
                       as.character(out.sub$method.p[i]), ".rds"))
  res[[i]] <- data.frame(method = as.character(out.sub$method)[i],
                         n1_n2 = paste0(out.sub$get_data.n1[i],".",
                                        out.sub$get_data.n2[i]),
                         pval = fl$p,
                         stringsAsFactors = F)
}
names(res) <- as.character(out.sub$method)

cols <- c("gray50", "forestgreen", "blue", "orange")
par(mfrow=c(2,2))
for (i in 1:length(res)) {
  hist(res[[i]]$pval, main="",
       xlab = "p-values", ylab = "Frequency",
       nclass = 20, col=cols[i])
  title(main=names(res)[i])
}

qq <- lapply(1:length(res), function(i) {
  qqplot(x=runif(100,0,1), y=res[[i]]$pval, plot.it=F)
})
plot(qq[[1]]$x, qq[[1]]$y, col = "gray50", cex=.7, pch = 16,
     xlab = "Uniform(0,1)", ylab = "Empirical distribution",
     main = "QQ-plot")
points(qq[[2]]$x, qq[[2]]$y, col = "forestgreen", cex=.7, pch = 16)
points(qq[[2]]$x, qq[[2]]$y, col = "blue", cex=.7, pch = 16)
points(qq[[3]]$x, qq[[3]]$y, col = "orange", cex=.7, pch = 16)
abline(0,1, col = "black")
title("random gene; 50 & 50", outer=TRUE, line=-1)

