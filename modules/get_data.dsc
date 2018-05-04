random_sample: R(counts = readRDS(${data_file})) + \
               get_data.R + \
               R(splitted = get_data_random(counts, args))
  n1, n2: (50, 50)
  p: 100
  @ALIAS: args = List(!data_file, !meta_file)
  $Y1: splitted$x
  $Y2: splitted$y

celltype_sample(random_sample):  R(counts = readRDS(${data_file}); \
                                 meta = readRDS(${meta_file})) + \
                                 get_data.R + \
                                 R(splitted = get_data_celltypes(counts, meta, args))
  groups: ("CD8 T cells", "CD14+ Monocytes")
