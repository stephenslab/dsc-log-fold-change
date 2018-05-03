get_data = function(counts, args) {
    # each row is a gene
    # each column is a sample
    counts = as.matrix(counts)
    n = ncol(counts)
    p = nrow(counts)
    args$n1 = min(n, args$n1)
    args$n2 = min(n, args$n2)
    args$p = min(args$p, p)
    genes_idx = sample(1:p, args$p)
    group1 = counts[genes_idx, sample(1:n, args$n1)]
    group2 = counts[genes_idx, sample(1:n, args$n2)]
    return(list(x=group1, y=group2, p=p))
}
