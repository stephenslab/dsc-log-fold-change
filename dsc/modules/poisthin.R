library(parallel)
#' Apply Poisson thinning to a matrix of count data.
#'
#' Given a matrix of RNA-seq counts, this function will randomly select two groups of
#' samples and add signal to a known proportion of the genes. This signal
#' is the log (base 2) effect size of the group indicator in a linear model.
#' The user may specify the distribution of the effects.
#'
#' The Poisson thinning approach first randomly assigns samples to be in one of two groups. Then,
#' given this assignment, will Binomially sample counts with a sample size of the gene expression
#' counts and a probability that is a function of the effect size. For details, see
#' Gerard and Stephens (2017).
#'
#' @param mat A matrix of count data. The rows index the individuals and
#'     the columns index the genes.
#' @param nsamp The number of samples to select from \code{mat}.
#' @param ngene The number of genes to select from \code{mat}.
#' @param gselect How should we select the subset of genes? Should we choose
#'     the \code{ngene} most median expressed genes (\code{"max"}), a random sample
#'     of the genes (\code{"random"}), a random sample of the most expressed
#'     genes (\code{"rand_max"}), a user-provided list (\code{"custom"}), or by maximum
#'     mean expression level (\code{"mean_max"})?
#'     If \code{"custom"}, then \code{gvec} should be specified. Expression levels
#'     of a gene are measured by median expression across individuals with ties broken
#'     by mean expression.
#' @param gvec A logical of length \code{ncol(mat)}. A \code{TRUE} in position \eqn{i}
#'     indicates inclusion into the smaller dataset. Hence, \code{sum(gvec)} should
#'     equal \code{ngene}.
#' @param shuffle_sample TRUE to shuffle sample per gene. Suppose samples are correlated
#'  across genes, shuffling sample labels at each gene would reduce sample correlation.
#' @param signal_dist What's the distribution of beta - the true effecrs? Options
#' 	include \code{"big_normal"} and \code{"near_normal"}.
#' @param signal_params Specify parameters for beta distribution. This is a list of
#' 	arguments to pass to \code{signal_dist}.
#' @param skip_gene The number of maximally expressed genes to skip.
#'     Not used if \code{gselect = "custom"}.
#' @param prop_null The proportion of genes that are null.
#' @param alpha If \eqn{b} is an effect and \eqn{s} is an empirical standard deviation, then
#'     we model \eqn{b/s^\alpha} as being exchangeable.
#' @param ncores Default 2 cores used in parallelized computing. Use mcapply in parallel package.
#'
#' @return A list with the following elements:
#' \itemize{
#'  \item{\code{Y}: }{A matrix of altered counts with \code{nsamp} rows
#'        and \code{ngene} columns.}
#'  \item{\code{X}: }{A design matrix. The first column contains a vector ones (for an
#'        intercept term) and the second column contains an indicator for group membership.}
#'  \item{\code{beta}: }{The approximately true effect sizes of \eqn{log(Y) ~ X\beta}.}
#' }
#'
#' @author David Gerard, Joyce Hsiao
#'
#' @export
poisthin <- function(mat, nsamp = nrow(mat), ngene = ncol(mat),
                     shuffle_sample = F,
                     gselect = c("max", "random", "rand_max", "custom", "mean_max"),
                     gvec = NULL,
                     skip_gene = 0,
                     signal_fun = c("big_normal","near_normal"),
                     signal_params = list(betapi=1, betamu=0, betasd=1),
                     prop_null = 1,
                     alpha = 0,
                     ncores=2) {
  ## Check Input -------------------------------------------------------------
#  assertthat::assert_that(is.matrix(mat))
  assertthat::assert_that(nsamp <= nrow(mat))
  assertthat::assert_that(ngene + skip_gene <= ncol(mat))
#  assertthat::assert_that(is.function(signal_fun))
  assertthat::assert_that(is.list(signal_params))
  assertthat::assert_that(prop_null >= 0, prop_null <= 1)

  gselect <- match.arg(gselect)

  if (gselect == "custom") {
    assertthat::assert_that(is.logical(gvec))
    assertthat::are_equal(length(gvec), ncol(mat))
    assertthat::are_equal(sum(gvec), ngene)
  } else {
    if (!is.null(gvec)) {
      warning('gvec is specified but being ignored since gselect is not "custom"')
    }
  }

  ## subset matrix -----------------------------------------------------------
  med_express <- apply(mat, 2, stats::median)
  mean_express <- colMeans(mat)
  order_vec <- order(med_express, mean_express, decreasing = TRUE)

  if (gselect == "max") {
    gindices <- order_vec[(skip_gene + 1):(skip_gene + ngene)]
  } else if (gselect == "rand_max") {
    first_zero <- match(0, med_express)
    max_gene <- min(c(first_zero, 2 * ngene + skip_gene, ncol(mat)), na.rm = TRUE)
    if (max_gene < ngene + skip_gene) {
      warning("including some low-expressed genes in the sample.")
      max_gene <- ngene + skip_gene
    }
    gindices <- sample(x = order_vec[(skip_gene + 1):max_gene], size = ngene)
  } else if (gselect == "random") {
    gindices <- sample(x = sample(order_vec[(skip_gene + 1):(skip_gene + ngene)]), size = ngene)
  } else if (gselect == "custom") {
    gindices <- (1:ncol(mat))[gvec]
    if (skip_gene > 0) {
      warning('ignoring skip_gene because gselect = "custom"')
    }
  } else if (gselect == "mean_max") {
    order_vec_means <- order(mean_express, decreasing = TRUE)
    gindices <- order_vec_means[(skip_gene + 1):(skip_gene + ngene)]
  }

  samp_indices <- sample(1:nrow(mat), size = nsamp)

  gindices <- sort(gindices)
  samp_indices <- sort(samp_indices)

  submat <- mat[samp_indices, gindices, drop = FALSE]
  group_indicator <- rep(FALSE, length = nsamp)
  group_indicator[sample(1:nsamp, size = floor(nsamp / 2))] <- TRUE

  ## Shuffle sample labels per gene or not
  if (shuffle_sample==TRUE) {
    submat <- do.call(rbind, mclapply(1:nrow(submat), function(g) {
      submat[g,sample(ncol(submat))]
    }, mc.cores = ncores))
  } else{
    submat <- submat
  }

  ## Draw signal -------------------------------------------------------------
#  if (signal_fun=="big_normal") {
#    signal_params <- list(betapi=1, betamu=0, betasd=1)
#  } else {
#    signal_params <- list(betapi=c(2/3,1/3), betamu=c(0,0),
#                          betasd=c(1,2)/sqrt(2*nsamp-2))
#  }

  make_normalmix <- function(nsignal, beta_args) {
    k <- length(beta_args$betapi) # number of components
    comp <- sample(1:k, nsignal, beta_args$betapi,replace=TRUE) #randomly draw a component
    beta <- rnorm(nsignal, beta_args$betamu[comp], beta_args$betasd[comp])
    return(beta)
  }

  nsignal <- round(ngene * (1 - prop_null))
  if (nsignal > 0) {
    signal_params$nsignal <- nsignal
#    signal_vec      <- do.call(what = signal_fun, args = signal_params) ## log2-fold change
    signal_vec <- make_normalmix(nsignal, signal_params)

    assertthat::are_equal(length(signal_vec), nsignal)

    which_signal <- sort(sample(1:ncol(submat), nsignal)) # location of signal

    ## Deal with alpha here ----------------------------
    if (abs(alpha) > 10 ^ -6) {
      sd_vec <- apply(log2(submat[, which_signal, drop = FALSE] + 1), 2, stats::sd) / sqrt(nrow(submat))
      assertthat::are_equal(length(sd_vec), length(signal_vec))
      signal_vec <- signal_vec * (sd_vec ^ alpha)
    }

    sign_vec  <- sign(signal_vec) # sign of signal
    bin_probs <- 2 ^ -abs(signal_vec) # binomial prob

    submat[group_indicator, which_signal[sign_vec > 0]] <-
      matrix(stats::rbinom(n = sum(sign_vec > 0) * nsamp / 2,
                           size = c(submat[group_indicator, which_signal[sign_vec > 0]]),
                           prob = rep(bin_probs[sign_vec > 0], each = nsamp / 2)),
             nrow = nsamp / 2)

    submat[!group_indicator, which_signal[sign_vec < 0]] <-
      matrix(stats::rbinom(n = sum(sign_vec < 0) * nsamp / 2,
                           size = c(submat[!group_indicator, which_signal[sign_vec < 0]]),
                           prob = rep(bin_probs[sign_vec < 0], each = nsamp / 2)),
             nrow = nsamp / 2)

    # for (gn in 1:length(signal_vec)) {
    #   if (sign_vec[gn] == 1) {
    #     current_count <- submat[group_indicator, which_signal[gn]]
    #     submat[group_indicator, which_signal[gn]] <-
    #       sapply(current_count, FUN = stats::rbinom, n = 1, prob = bin_probs[gn])
    #   } else if (sign_vec[gn] == -1) {
    #     current_count <- submat[!group_indicator, which_signal[gn]]
    #     submat[!group_indicator, which_signal[gn]] <-
    #       sapply(current_count, FUN = stats::rbinom, n = 1, prob = bin_probs[gn])
    #   }
    # }
    beta <- rep(0, ngene)
    beta[which_signal] <- -1 * signal_vec ## -1 because of way design matrix is created
  } else if (nsignal == 0 & abs(prop_null - 1) > 10 ^ -6) {
    warning('no genes were given signal since (1 - prop_null) * ngene was very close to zero')
    beta <- rep(0, ngene)
  } else {
    beta <- rep(0, ngene)
  }

  X <- stats::model.matrix(~group_indicator)
  return_list <- list(Y = submat, X = X, beta = beta)

  return(return_list)
}





#log2_cpm <- function() {}



# older code ----------


library(Matrix)

sample_data <- function(counts, n, p=NULL) {
    # each row is a gene
    # each column is a sample
    counts = as.matrix(counts)
    nn = ncol(counts)
    ngenes = nrow(counts)
    n = min(nn, n)

    if (!is.null(p)) {
      pp = min(ngenes, p)
      genes_idx = sample(1:ngenes, pp)
      }
    if (is.null(p)) {
      pp = ngenes
      genes_idx = 1:ngenes
      }

    return(counts[genes_idx, sample(1:nn, n)])
}

get_sample_correlated <- function(counts, args) {
  # if (is.null(args$seed)) {args$seed <- 99}
  # set.seed(args$seed)
  df <- sample_data(counts, args$n1+args$n2, args$p)
  group <- c(rep(1, args$n1), rep(2, args$n2))
  df.perm <- do.call(rbind, lapply(1:nrow(df), function(g) {
      df[g,sample(ncol(df))]
  }))
  x <- df.perm[,group==1]
  y <- df.perm[,group==2]
  return(list(x=x,
              y=y))
}

get_sample_uncorrelated <- function(counts, args) {
  # if (is.null(args$seed)) {args$seed <- 99}
  # set.seed(args$seed)
  df <- sample_data(counts, args$n1+args$n2, args$p)
  group <- c(rep(1, args$n1), rep(2, args$n2))
  x <- df[,group==1]
  y <- df[,group==2]
  return(list(x=x,
              y=y))
}


