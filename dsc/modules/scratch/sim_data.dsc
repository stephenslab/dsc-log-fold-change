sim_data: R(counts <- readRDS(${dataFile})) +\
          R(sim_counts <- poisthin(counts, nsamp=)) + \


#sim_data_berge_params: R(counts=readRDS(${dataFile})) + \
#              sim_data.R + \
#              R(params = sim_data_berge_params(counts=counts))
#  @ALIAS: args = List(!dataFile)
#  $params: params

#sim_data_berge: R(counts = readRDS(${dataFile})) + \
#               sim_data.R + \
#               R(set.seed(args$seed); output = sim_data_berge(counts, args)) + \
#               R(Y1 <- output$Y[,output$group==1]; Y2 <- output$Y[,output$group==2])
#  seed: R{2:3}
#  pi1: .1
#  param: $params
#  p: 1000
#  @ALIAS: args = List(!dataFile)
#  $Y1: Y1
#  $Y2: Y2
#  $DEind: output$indDE

