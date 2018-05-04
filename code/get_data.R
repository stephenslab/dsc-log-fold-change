sample_data = function(counts, n, p) {
    # each row is a gene
    # each column is a sample
    counts = as.matrix(counts)
    nn = ncol(counts)
    pp = nrow(counts)
    n = min(nn, n)
    p = min(pp, p)
    genes_idx = sample(1:pp, p)
    return(counts[genes_idx, sample(1:nn, n)])
}

get_data_random = function(counts, args) {
  return(list(x=sample_data(counts, args$n1, args$p),
              y=sample_data(counts, args$n2, args$p)))
}

get_data_celltypes = function(counts, meta, args) {
  counts = as.matrix(counts)
  group1 = counts[ , which(meta$celltype == args$groups[1])]
  group2 = counts[ , which(meta$celltype == args$groups[2])]
  return(list(x=sample_data(group1, args$n1, args$p),
              y=sample_data(group2, args$n2, args$p)))
}
