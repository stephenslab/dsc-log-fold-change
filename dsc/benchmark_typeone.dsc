#!/usr/bin/env dsc

%include modules/dataSimulate
%include modules/methodsMeanExpression

DSC:
  define:
    data: data_poisthin_null, data_poisthin_signal
    method: edger, deseq2, glm_pois, glm_quasipois, limma_voom, mast, t_test, wilcoxon
#    score: typeone, fdr
  run:
    data_poisthin_null * method
  exec_path: modules
  global:
    dataFile: "data/pbmc_counts.rds"
  output:
    /scratch/midway2/joycehsiao/dsc-log-fold-change/benchmark_typeone


