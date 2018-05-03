## PBMC example dataset

This is an example dataset of single cell RNA-seq data from 10x genomics technology. The single cell samples include CD14+ Monocytes and CD8 T cells, a total of 787 samples and 946 genes. The original data include 2,700 single cells and 32,738 genes from 8 immune cell types: B cells, CD4 T cells, CD8 T cells, CD14+ Monocytes, Dendritic cells, FCGR3A+ Monocytes, Megakaryocytes, NK cells. We filtered the original data to include samples that are detected in > 200 genes and genes that are detected in > 20% of cells.

How to use the example data:

Import raw counts

```counts <- readRDS("~/data/rawcounts.rds")```

Import metadata

```metadata <- readRDS("~/data/metadata.rds")```

Check sample size of each group

```table(metadata$celltype)```


How to cite the data:

Zheng, G. X. Y. et al. Massively parallel digital transcriptional profiling of single cells. Nat. Commun. 8, 14049 doi: 10.1038/ncomms14049 (2017)


How to obtain the full dataset:

See the Seurat tutrial: http://satijalab.org/seurat/pbmc3k_tutorial.html
