#' @description Simulation expression based on Van den Berge et al., 2018

source("/project2/gilad/joycehsiao/dsc-log-fold-change/dsc/code/zinbwaveZinger/zingeRsimulationFunctions/simulationHelpFunctions_v7_diffInZero.R")

sim_data_berge_params <- function(counts) {
  # if (is.null(args$seed)) {args$seed <- 99}
  # set.seed(args$seed)
  params <- getDatasetMoMPositive(counts = counts)

  return(params)
}


sim_data_berge <- function(counts, args) {
  # if (is.null(args$seed)) {args$seed <- 99}
  # set.seed(args$seed)
#  params <- getDatasetMoMPositive(counts = counts)

  nSamples <- ncol(counts)
  grp <- as.factor(rep(1:2, each = nSamples/2)) #two-group comparison
  nTags <- args$p #nr of features
#  set.seed(args$seed)
  DEind <- sample(1:nTags,floor(nTags*args$pi1),replace=FALSE) #10% DE
  fcSim <- (2 + rexp(length(DEind), rate = 1/2)) #fold changes
  libSizes <- sample(colSums(counts),nSamples,replace=TRUE) #library sizes
  simData <- NBsimSingleCell(foldDiff = fcSim, ind = DEind,
                             dataset = counts, nTags = nTags,
                             group = grp,
                             verbose = TRUE, params = args$params,
                             lib.size = libSizes, cpm="AveLogCPM", normalizeLambda=TRUE,
                             min.dispersion=1e-3)

  return(simData)
}
