# dsc-two-sample-comparisons

This is (or will be) a Dynamic Statistical Comparison
to compare different approaches to comparing two groups, from count data.

Our intention is to initially focus on data from single-cell experiments.

# The goal:

The goal is to compare methods for estimating the log-fold-change between two groups.

Methods will input:

  - $Y1 a numeric matrix of data from group 1 (n by p, n samples/cells in rows, p genes in columns) 
  - $Y2 a numeric matrix of data from group 2 (n by p, n samples/cells in columns, p genes in columns) 
  
and optionally:

  - $X1 an n vector of covariates from group 1 (eg "library size")
  - $X2 an n vector of covariates from group 2

Methods will output: 

  - $log_fold_change a vector of estimates of log(mu1/mu2) for each gene where mu1 is the mean of group 1 and mu2 is the mean of group 2
  - $s_hat a vector of standard error for the estimated $log_fold_change
  - $p a p-vector of p values testing whether each log-fold change is 0

We will create synthetic data (from real data) that have
known log-fold-change values, and compare the estimates with the real values.
We will also assess calibration of p values (eg on null data, we should get uniform p values)
and power.

# Data 

We will start by creating null data by sampling two group at random from real data.)


# Methods

List methods we might want to use...



