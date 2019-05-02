#' @param Y gene by sample
#'
normalize_log2cpm <- function(Y) {
  libsize <- base::colSums(Y)
  log2cpm <- log2(t(10^6*(t(Y+1)/libsize)))
  log2cpm_qqnormed <- do.call(rbind,
                              lapply(1:nrow(log2cpm), function(g) {qqnorm(log2cpm[g,], plot.it=F)$x}))
  return(log2cpm_qqnormed)
}




# normalize <- function(matExpr=NULL,eSet=NULL,geneLen=NULL,log=FALSE,method=c("quantile","rpkm","Cquant","cpm")) {
#   #eSet = eSetNone
#
#   suppressMessages(library(limma))
#   suppressMessages(library(DESeq))
#   suppressMessages(library(edgeR))
#
#
#   if (method=="cpm") {
#     # 1. CPM normalzation
#     matExpr = exprs(eSet)
#
#     libSize = colSums(matExpr,na.rm=T)
#
#     mat.cpm = t(t(matExpr)/libSize)*(10^6)
#
#     mat = log2(mat.cpm)
#     return(mat)
#   }
#
#   if (method=="quantile") {
#     if (log==TRUE) { mat = normalizeQuantiles(log2(matExpr)) }
#     if (log==FALSE) { mat = normalizeQuantiles(matExpr) }
#
#     return(mat)
#   }
#
#   if (method=="Cquant") {
#     # 1. CPM normalzation
#     # 2. Convert to RPKM
#     # 3. quantile normalization
#     matExpr = exprs(eSet)
#
#     ii = match(featureNames(eSet),geneLen$ENSGID)
#     length = geneLen[ii,2]
#
#     libSize = colSums(matExpr)
#
#     mat.cpm = t(t(matExpr)/libSize)*(10^6)
#     mat.rpkm = (mat.cpm/length)*1000
#
#     mat = normalizeQuantiles(log2(mat.rpkm+1))
#     return(mat)
#   }
#
#   if (method=="rpkm") {
#     #     eSet=eSetNone[,eSetNone$seqData=="rna"]
#     #     geneLen=fData(eSet[which(fData(eSet)$ENSGID %in% fData(eSetNone)$ENSGID),])
#
#     matExpr = exprs(eSet)
#
#     ii = match(featureNames(eSet),geneLen$ENSGID)
#     length = geneLen[ii,2]
#
#     libSize = colSums(matExpr,na.rm=T)
#
#     mat.cpm = t(t(matExpr)/libSize)*(10^6)
#     mat.rpkm = (mat.cpm/length)*1000
#
#     #     if (log==TRUE) {mat = log2(mat.rpkm+1) }
#     #     if (log==FALSE) {mat = mat.rpkm}
#     return(mat.rpkm)
#   }
#
#   #   if (method=="DESeq") {
#   #     cds = estimateSizeFactorsForMatrix(matExpr)
#   #     dmat = t(t(matExpr)/cds)
#   #     dmat = log2(dmat+1)
#   #     return(dmat)
#   #   }
#
# }
