# testing scripts ------

dsc benchmark.dsc --target pipe_power_choose_betasd --truncate --host config.yaml  -o /scratch/midway2/joycehsiao/dsc-log-fold-change/pipe_power_choose_betasd

dsc benchmark.dsc --target pipe_type1 --truncate --host config.yaml  -o /scratch/midway2/joycehsiao/dsc-log-fold-change/pipe_type1_test

dsc benchmark.dsc --target pipe_power --truncate --host config.yaml  -o /scratch/midway2/joycehsiao/dsc-log-fold-change/pipe_power_test

# running scripts

dsc benchmark.dsc --target pipe_power_choose_betasd --host config.yaml  -o /scratch/midway2/joycehsiao/dsc-log-fold-change/pipe_power_choose_betasd

dsc benchmark.dsc --target pipe_type1 --host config.yaml  -o /scratch/midway2/joycehsiao/dsc-log-fold-change/pipe_type1

dsc benchmark.dsc --target pipe_power --host config.yaml  -o /scratch/midway2/joycehsiao/dsc-log-fold-change/pipe_power




# else

#dsc benchmark.dsc --target "make_null_berge*process_data*t_test" --host config.yaml -o /scratch/midway2/joycehsiao/dsc-log-fold-change/type1_berge

#dsc benchmark.dsc --target "sample_correlated*t_test" --host config.yaml -o /scratch/midway2/joycehsiao/dsc-log-fold-change/type1_simple

#dsc type1_berge.dsc --host config.yaml #-o /scratch/midway2/joycehsiao/dsc-log-fold-change/type1_berge


dsc benchmark.dsc --target pipe_power --host config.yaml -o /scratch/midway2/joycehsiao/dsc-log-fold-change/benchmark

dsc benchmark.dsc --target pipe_null --host config.yaml -o /scratch/midway2/joycehsiao/dsc-log-fold-change/benchmark
