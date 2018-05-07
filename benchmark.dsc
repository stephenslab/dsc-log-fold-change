#!/usr/bin/env dsc

%include modules/get_data
%include modules/glm
%include modules/t_test
%include modules/wilcox_test
%include modules/run_DESeq2

DSC:
  define:
    get_data: random_sample, random_gene #, celltype_sample
    method: glm_pois, t_test, wilcoxon_test #DESeq2
  run:
    first_pass: get_data * method
  exec_path: modules
  global:
    data_file: data/rawcounts.rds
    meta_file: data/metadata.rds
