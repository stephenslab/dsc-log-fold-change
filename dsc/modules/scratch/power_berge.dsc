#!/usr/bin/env dsc

%include modules/get_sim_data
%include modules/methodsMeanExpression

DSC:
  define:
   get_data: get_sim_data
   method: edger, deseq2, glm_pois, glm_quasipois, limma_voom, mast, t_test, wilcoxon
#    method: edger, deseq2, glm_pois, glm_quasipois, limma_voom, mast, t_test, wilcoxon, zinbwave_edger, zinbwave_deseq2
  run:
    get_data*method
  exec_path: modules
  global:
    dataFile: "data/pbmc_simdata_berge.rds"
  output:
    /scratch/midway2/joycehsiao/dsc-log-fold-change/power_berge


