data_poisthin_signal: R(counts = readRDS(${dataFile})) + \
       dataSimulate.R + \
       R(out = poisthin(mat=t(counts), nsamp = args$nsamp, ngene = args$ngene, gselect = args$gselect, signal_dist = "bignormal", prop_null = args$prop_null)) + \
       R(groupInd = out$X[,2]; Y1 = t(out$Y[groupInd==1,]); Y2 = t(out$Y[groupInd==0,]))
  seed: R{2:101}
  nsamp: 90
  ngene: 1000
  prop_null: .9
  gselect: "random"
  signal_dist: "bignormal"
  @ALIAS: args = List(!dataFile)
  $Y1: Y1
  $Y2: Y2
  $beta: out$beta


data_poisthin_null: R(counts = readRDS(${dataFile})) + \
     dataSimulate.R + \
     R(out = poisthin(mat=t(counts), nsamp = args$nsamp, ngene = args$ngene, gselect = args$gselect, signal_dist = args$signal_dist, prop_null = args$prop_null)) + \
     R(groupInd = out$X[,2]; Y1 = t(out$Y[groupInd==1,]); Y2 = t(out$Y[groupInd==0,]))
  seed: R{2:101}
  nsamp: 90
  ngene: 1000
  prop_null: 1
  gselect: "random"
  signal_dist: "big_normal"
  @ALIAS: args = List(!dataFile)
  $Y1: Y1
  $Y2: Y2
  $beta: out$beta


sample_correlated: R(counts = readRDS(${dataFile})) + \
               get_data.R + \
               R(set.seed(args$seed); splitted = get_sample_correlated(counts, args))
  seed: R{2:101}
  n1, n2: (45, 45), (300, 300)
  p: 1000
  @ALIAS: args = List(!dataFile)
  $Y1: splitted$x
  $Y2: splitted$y

sample_uncorrelated(sample_correlated): R(counts = readRDS(${dataFile})) + \
               get_data.R + \
               R(set.seed(args$seed); splitted = get_sample_uncorrelated(counts, args))


make_null_berge: R(counts = readRDS(${dataFile})) + \
               get_data.R + \
               R(set.seed(args$seed); splitted = get_sample_correlated(counts, args))
  seed: R{2:101}
  n1, n2: (45, 45) #, (300, 300)
  p: NULL
  @ALIAS: args = List(!dataFile)
  $Y1: splitted$x
  $Y2: splitted$y


#get_sim_data: R(counts = readRDS(${dataFile})) + \
#              R(celltype = readRDS(${metaFile})) + \
#              R(Y1 <- counts[, celltype==1]; Y1 <- counts[, celltype==2])
#  @ALIAS: args = List(!dataFile, !metaFile)
#  $Y1: Y1
#  $Y2: Y2

