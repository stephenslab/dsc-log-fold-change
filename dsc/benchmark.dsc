#!/usr/bin/env dsc

%include modules/get_data
%include modules/methodsMeanExpression

DSC:
  define:
    get_data: sample_correlated, sample_uncorrelated
#    method: glm_pois, glm_quasipois, t_test, wilcoxon, deseq2, limma_voom
    method: t_test, limma_voom
  run:
    first_pass: get_data * method
  exec_path: modules
  global:
    data_file: "data/rawcounts.rds"
    meta_file: "data/metadata.rds"
