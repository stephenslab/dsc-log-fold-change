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


