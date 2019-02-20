## PBMC example dataset

This is an example dataset of single cell RNA-seq data from 10x genomics technology. The single cell samples include CD14+ Monocytes and CD8 T cells, a total of 787 samples and 946 genes. The original data include 2,700 single cells and 32,738 genes from 8 immune cell types: B cells, CD4 T cells, CD8 T cells, CD14+ Monocytes, Dendritic cells, FCGR3A+ Monocytes, Megakaryocytes, NK cells. 

We filtered the original data to include samples that are detected in > 200 genes and genes that are detected in > 20% of cells.

This folder contains

* A small dataset for evaluating type I error
    + This is a subset of two cell types from the PBMC data
    + CD8 T cells (308) and CD14+ Monocytes (479) 
    + `data/rawcounts.rds` 
    + `data/metadata.rds` 

* A large dataset for evaluating type I error
    + expressionSet object in `data/pbmc.rds` and count matrix in `data/pbmc_counts.rds`
    + Follow Van den Berge et al., 2018 for data curation steps. See [here]((dsc-log-fold-change/docs/pbmc_berge_null.html)) for the codes.
    
* A small dataset for evaluating power and false positive rate control
    + count matrix in `data/pbmc_counts_sub.rds`
    + Extracted from Analysis documented [here](dsc-log-fold-change/docs/pbmc_berge_null.html)

```{r}
df <- readRDS("/project/mstephens/data/external_public_supp/pbmc3k_filtered_gene_bc_matrices/pbmc3k_final.rds?dl=1")
pbmc_counts <- df@data
```

To cite this data:

Zheng, G. X. Y. et al. Massively parallel digital transcriptional profiling of single cells. Nat. Commun. 8, 14049 doi: 10.1038/ncomms14049 (2017)

To retrieve the original data from the source:

See the Seurat tutrial: http://satijalab.org/seurat/pbmc3k_tutorial.html
