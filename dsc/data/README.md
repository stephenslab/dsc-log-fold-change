## Experimental data 

* PBMC data of 8+ cell types, 13,713 genes from 2,638 samples
    + Location: 'data/pbmc_counts.rds'
    + Experiment: single cell RNA-seq from 10X genomics  
    + Source:
      1. Downloaded from Seurat tutrial: http://satijalab.org/seurat/pbmc3k_tutorial.html
      2. This is a subset of data from Zheng, G. X. Y. et al. Massively parallel digital transcriptional profiling of single cells. Nat. Commun. 8, 14049 doi: 10.1038/ncomms14049 (2017)
    + Code: 
      1. Filter the downloaed samples following Van den Berge et al., 2018. See code [here]((dsc-log-fold-change/docs/pbmc_berge_null.html)). 
    + Others:
      1. `data/pbmc.rds` (an ExpressionSet object) for associated meta data

* PBMC data of 2 cell types, 13,713 genes from 1,153 samples
    + Location: `data/pbmc_counts_sub.rds`
    + Experiment: single cell RNA-seq from 10X genomics
    + Source: A subset of `data/pbmc_counts_sub.rds`

* PBMC data of 2 cell types, 946 genes from 787 samples
    + Location: `data/rawcounts.rds`, `data/metadata.rds` 
    + Experiment: single cell RNA-seq from 10X genomics  
    + CD8 T cells (308) and CD14+ Monocytes (479) 
    + Description: The single cell samples include CD14+ Monocytes and CD8 T cells, a total of 787 samples and 946 genes. The original data include 2,700 single cells and 32,738 genes from 8 immune cell types: B cells, CD4 T cells, CD8 T cells, CD14+ Monocytes, Dendritic cells, FCGR3A+ Monocytes, Megakaryocytes, NK cells. We filtered the original data to include samples that are detected in > 200 genes and genes that are detected in > 20% of cells.




