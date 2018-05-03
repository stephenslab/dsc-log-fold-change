get_data: R(counts = readRDS(${data_file})) + \
          get_data.R + \
          R(splitted = get_data(counts, args))
  n1, n2: (50, 50)
  p: 100
  @ALIAS: args = List(!data_file)
  $Y1: splitted$x
  $Y2: splitted$y
