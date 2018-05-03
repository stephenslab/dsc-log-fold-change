# dsc-log-fold-change

This is (or will be) a Dynamic Statistical Comparison
to estimating (or testing) the "log-fold-change" in mean between
two groups from count data.

Our intention is to initially focus on data from single-cell experiments.

# The goal:

The goal is to compare methods for estimating the log-fold-change between two groups.

Methods will input:

  - $Y1 a numeric matrix of data from group 1 (p by n, p genes in columns, n samples/cells in rows)
  - $Y2 a numeric matrix of data from group 2 (p by n, p genes in columns, n samples/cells in rows)
  
and optionally:

  - $X1 an n vector of covariates from group 1 (eg "library size")
  - $X2 an n vector of covariates from group 2

Methods will output:

  - $log_fold_change_est a vector of estimates of log(mu1/mu2) for each gene where mu1 is the mean of group 1 and mu2 is the mean of group 2
  - $s_hat a vector of standard error for the estimated $log_fold_change
  - $p a p-vector of p values testing whether each log-fold change is 0

# Data

We will create synthetic data (from real data) that have
known log-fold-change values, and compare the estimates with the real values.
We will also assess calibration of p values (eg on null data, we should get uniform p values)
and power.

To create data we will take a file containing count data and select
samples at random to create two groups. These will be "null" data.

Input:
  - file of data
  - n1 sample size for group 1
  - n2 sample size for group 2
  - p number of genes
  - pi0 proportion of nulls
  - g a distribution on non-zero effects

Output:
  - $Y1
  - $Y2
  - $log-fold-change (true value of log_fold_change for each gene)

# Methods

List methods we might want to use...

# Run DSC

The main DSC file is `benchmarks.dsc`. To see what is available:

```
./benchmark.dsc -h
```

and to run the benchmark:

```
./benchmark.dsc
```

Or to run a minimal test benchmark, eg

```
./benchmark.dsc --target "get_data * wilcoxon_test" --truncate --replicate 1 # default is in fact --replicate 1
```
