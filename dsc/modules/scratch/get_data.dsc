sim_data: R(counts = readRDS(${dataFile})) + \
       get_data.R + \
       R(sim_counts = poisthin(counts, args))
  seed: R{2:3}
  nsamp: 90
  ngene: 1000
  prop_null: 1
  gselect: "random"
  @ALIAS: args = List(!dataFile)
  $Y: out$Y
  $X: out$X
  $beta: out$beta

library(seqgendiff)
out <- poisthin(t(a), nsamp=90, ngene=100,
gselect="random", gvec=NULL, signal_fun=stats::rnorm,
signal_params=list(mean=0,sd=1),
prop_null=1)

#  n1, n2: (45, 45), (300, 300)
#  p: 1000
#  @ALIAS: args = List(!dataFile)
#  $Y1: splitted$x
#  $Y2: splitted$y





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


get_sim_data: R(counts = readRDS(${dataFile})) + \
              R(celltype = readRDS(${metaFile})) + \
              R(Y1 <- counts[, celltype==1]; Y1 <- counts[, celltype==2])
  @ALIAS: args = List(!dataFile, !metaFile)
  $Y1: Y1
  $Y2: Y2

