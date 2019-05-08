# testing scripts ------

dsc benchmark.dsc --target pipe_power_choose_betasd --truncate --host config.yaml  -o /scratch/midway2/joycehsiao/dsc-log-fold-change/pipe_power_choose_betasd

dsc benchmark.dsc --target pipe_type1 --truncate --host config.yaml  -o /scratch/midway2/joycehsiao/dsc-log-fold-change/pipe_type1

dsc benchmark.dsc --target pipe_power --truncate --host config.yaml  -o /scratch/midway2/joycehsiao/dsc-log-fold-change/pipe_power

dsc benchmark.dsc --target pipe_type1_libsize --truncate --host config.yaml  -o /scratch/midway2/joycehsiao/dsc-log-fold-change/pipe_type1_libsize

dsc benchmark.dsc --target pipe_gtex --truncate --host config.yaml  -o /scratch/midway2/joycehsiao/dsc-log-fold-change/pipe_gtex

# running scripts

dsc benchmark.dsc --target pipe_choose_betasd --host config.yaml  -o /scratch/midway2/joycehsiao/dsc-log-fold-change/pipe_power_choose_betasd

dsc benchmark.dsc --target pipe_type1 --host config.yaml  -o /scratch/midway2/joycehsiao/dsc-log-fold-change/pipe_type1

dsc benchmark.dsc --target pipe_power --host config.yaml  -o /scratch/midway2/joycehsiao/dsc-log-fold-change/pipe_power

dsc benchmark.dsc --target pipe_type1_libsize --host config.yaml  -o /scratch/midway2/joycehsiao/dsc-log-fold-change/pipe_type1_libsize

dsc benchmark.dsc --target pipe_gtex --host config.yaml  -o /scratch/midway2/joycehsiao/dsc-log-fold-change/pipe_gtex



#
dsc benchmark.dsc --target pipe_choose_betasd --truncate --host config.yaml  -o /scratch/midway2/joycehsiao/dsc-log-fold-change/pipe_choose_betasd

dsc benchmark.dsc --target pipe_pbmc --replicate 2 --host config.yaml  -o /scratch/midway2/joycehsiao/dsc-log-fold-change/pipe_pbmc


dsc benchmark.dsc --target pipe_gtex --replicate 10 --host config.yaml  -o /scratch/midway2/joycehsiao/dsc-log-fold-change/pipe_gtex

dsc benchmark.dsc --target pipe_null_libsize --replicate 10 --host config.yaml  -o /scratch/midway2/joycehsiao/dsc-log-fold-change/pipe_null_libsize


# else

#dsc benchmark.dsc --target "make_null_berge*process_data*t_test" --host config.yaml -o /scratch/midway2/joycehsiao/dsc-log-fold-change/type1_berge

#dsc benchmark.dsc --target "sample_correlated*t_test" --host config.yaml -o /scratch/midway2/joycehsiao/dsc-log-fold-change/type1_simple

#dsc type1_berge.dsc --host config.yaml #-o /scratch/midway2/joycehsiao/dsc-log-fold-change/type1_berge


dsc benchmark.dsc --target pipe_power --host config.yaml -o /scratch/midway2/joycehsiao/dsc-log-fold-change/benchmark

dsc benchmark.dsc --target pipe_null --host config.yaml -o /scratch/midway2/joycehsiao/dsc-log-fold-change/benchmark
