get_sim_data: R(df = readRDS(${dataFile})) + \
              R(counts <- df$counts; group <- df$group; indDE <- df$indDE) + \
              R(out <- list(Y1= counts[, group==1], Y2 = counts[, group==2]))
  @ALIAS: args = List(!dataFile)
  $Y1: out$Y1
  $Y2: out$Y2

