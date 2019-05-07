## Experimental data 

* GTEx V6 bulk expression data
    + gtex_lung.rds: 320 samples and 16,069 genes.

* PBMC data of 8+ cell types, 13,713 genes from 2,638 samples
    + Location: 'data/pbmc_counts.rds'
    + Experiment: frozen human PBMCs, single cell RNA-seq from 10X genomics  
    + Summary (provided on 10x website): 
      1. Single Cell Gene Expression Dataset by Cell Ranger 1.1.0
      2. Frozen PBMCs from Donor A (Human)
      3. ~2,900 cells detected (raw data before filtering)
      4. Sequenced on Illumina NextSeq 500 High Output with ~25,000 reads per cell
      5. 98bp read1 (transcript), 8bp I5 sample barcode, 14bp I7 GemCode barcode and 5bp read2 (UMI)
      6. Analysis run with --expected-cells=3000 (barcode cut-off)
    + Reference: Zheng, G. X. Y. et al. Massively parallel digital transcriptional profiling of single cells. Nat. Commun. 8, 14049 doi: 10.1038/ncomms14049 (2017)
    + Downloaded from: Seurat tutrial: http://satijalab.org/seurat/pbmc3k_tutorial.html. 
      1. This dataset contains a filtered set of gene expression data. The original dataset contains 2,900 single-cell samples and can be found here https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/frozen_pbmc_donor_a (in cellranger format).
      2. See `data/pbmc.rds` (an ExpressionSet object) for associated meta data, such as cell type labels assigned in clustering analysis. 
          + Clustering analysis following Van den Berge et al., 2018. Reproduce results [here]((dsc-log-fold-change/docs/pbmc_berge_null.html)). 


* PBMC data of 2 cell types, 13,713 genes from 1,153 samples
    + Location: `data/pbmc_counts_sub.rds`
    + Experiment: single cell RNA-seq from 10X genomics
    + Source: A subset of `data/pbmc_counts_sub.rds`

* PBMC data of 2 cell types, 946 genes from 787 samples
    + Location: `data/rawcounts.rds`, `data/metadata.rds` 
    + Experiment: single cell RNA-seq from 10X genomics  
    + CD8 T cells (308) and CD14+ Monocytes (479) 
    + Description: The single cell samples include CD14+ Monocytes and CD8 T cells, a total of 787 samples and 946 genes. The original data include 2,700 single cells and 32,738 genes from 8 immune cell types: B cells, CD4 T cells, CD8 T cells, CD14+ Monocytes, Dendritic cells, FCGR3A+ Monocytes, Megakaryocytes, NK cells. We filtered the original data to include samples that are detected in > 200 genes and genes that are detected in > 20% of cells.




