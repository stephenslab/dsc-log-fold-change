filtering_gene: process_data.R + \
               R(genes_to_keep = filter_genes(Y1, Y2, per_gene_zero_cell=5)) + \
               R(Y1_filt <- Y1[genes_to_keep,]; Y2_filt <- Y2[genes_to_keep,])
  Y1: $Y1
  Y2: $Y2
  $Y1: Y1_filt
  $Y2: Y2_filt



