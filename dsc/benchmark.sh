#dsc benchmark.dsc --target "make_null_berge*process_data*t_test" --host config.yaml -o /scratch/midway2/joycehsiao/dsc-log-fold-change/type1_berge

#dsc benchmark.dsc --target "sample_correlated*t_test" --host config.yaml -o /scratch/midway2/joycehsiao/dsc-log-fold-change/type1_simple

#dsc type1_berge.dsc --host config.yaml #-o /scratch/midway2/joycehsiao/dsc-log-fold-change/type1_berge


dsc benchmark.dsc --target pipe_power --host config.yaml -o /scratch/midway2/joycehsiao/dsc-log-fold-change/benchmark

dsc benchmark.dsc --target pipe_null --host config.yaml -o /scratch/midway2/joycehsiao/dsc-log-fold-change/benchmark
