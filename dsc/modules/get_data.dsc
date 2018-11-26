sample_correlated: R(counts = readRDS(${data_file})) + \
               get_data.R + \
               R(splitted = get_sample_correlated(counts, args))
  seed: R{1:1000}
  n1, n2: (50, 50), (300, 300)
  p: 1000
  @ALIAS: args = List(!data_file, !meta_file)
  $Y1: splitted$x
  $Y2: splitted$y

sample_uncorrelated(sample_correlated): R(counts = readRDS(${data_file})) + \
               get_data.R + \
               R(splitted = get_sample_uncorrelated(counts, args))

#celltype_sample(random_sample):  R(counts = readRDS(${data_file}); \
#                                 meta = readRDS(${meta_file})) + \
#                                 get_data.R + \
#                                 R(splitted = get_data_celltypes(counts, meta, args))
#  groups: ("CD8 T cells", "CD14+ Monocytes")
