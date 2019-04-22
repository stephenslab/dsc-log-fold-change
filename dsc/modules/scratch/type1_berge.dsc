#!/usr/bin/env dsc

%include modules/get_data
%include modules/process_data
%include modules/methodsMeanExpression

DSC:
  define:
    get_data: make_null_berge
    process_data: filtering_gene
    method: edger, deseq2, glm_pois, glm_quasipois, limma_voom, mast, t_test, wilcoxon
    #zinbwave_edger, zinbwave_deseq2
  run:
    first_pass: get_data*process_data*method
  exec_path: modules
  global:
    dataFile: "data/pbmc_counts.rds"
  output:
    /scratch/midway2/joycehsiao/dsc-log-fold-change/type1_berge


