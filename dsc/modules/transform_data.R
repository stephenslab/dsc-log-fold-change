
#' @description transform expression count data
#' @param Y count data; sample by gene
#' @param X design matrix (sample by variables)
#'
transform_data <- function(Y, X=NULL, pseudo_count = 1, log=c("none", "log2", "log10"),
                           libscale_method = c("none", "sum", "TMM", "RLE",
                                               "pearsons_residual")) {
  libsize <- colSums(Y)

  if (log!="none") {
    if (libscale_method == "none") {
      transformed_Y <- log2(Y + pseudo_count)
    }
    if (libscale_method == "sum") {
      # assume all UMI data, and total library size ~ 100K
      transformed_Y <- log2(10e+05 * t(t(Y)/libsize) + pseudo_count)
    }
    if (libscale_method == "TMM") {
      libscale_factors <- edgeR::calcNormFactors(Y, method = libscale_method)
      libsize_normed <- libsize*libscale_factors
      transformed_Y <- log2(10e+05 * t(t(Y)/libsize_normed) + pseudo_count)
    }
    if (libscale_method == "RLE") {
      libscale_factors <- edgeR::calcNormFactors(Y, method = libscale_method)
      libsize_normed <- libsize*libscale_factors
      transformed_Y <- log2(10e+05 * t(t(Y)/libsize_normed) + pseudo_count)
    }
  }
  if (log=="none" & libscale_method == "pearsons_residual") {
#      library(sctransform)
      vst_out <- sctransform::vst(Y, latent_var = c('log_umi'),
                                  return_gene_attr = TRUE, return_cell_attr = TRUE,
                                  show_progress = FALSE)
      transformed_Y <- vst_out$y
  }
  # if (log=="none" & libscale_method == "TMM") {
  #   libscale_factors <- edgeR::calcNormFactors(Y, method = libscale_method)
  #   libsize_normed <- libsize*libscale_factors
  #   transformed_Y <- t(t(Y)/libsize_normed)
  # }
  # if (log=="none" & libscale_method == "RLE") {
  #   libscale_factors <- edgeR::calcNormFactors(Y, method = libscale_method)
  #   libsize_normed <- libsize*libscale_factors
  #   transformed_Y <- t(t(Y)/libsize_normed)
  # }

  return(transformed_Y)
}
