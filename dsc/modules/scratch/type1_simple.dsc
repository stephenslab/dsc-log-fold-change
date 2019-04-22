#!/usr/bin/env dsc

%include modules/get_data
%include modules/process_data
%include modules/methodsMeanExpression

DSC:
  define:
    get_data: sample_correlated, sample_uncorrelated
    method: edger, deseq2, glm_pois, glm_quasipois, limma_voom, mast, t_test, wilcoxon
    #zinbwave_edger, zinbwave_deseq2
  run:
    get_data*method
  exec_path: modules
  global:
    dataFile: "data/rawcounts.rds"
  output:
    /scratch/midway2/joycehsiao/dsc-log-fold-change/type1_simple


