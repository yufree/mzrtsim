
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mzrtsim

<!-- badges: start -->

[![Actions
Status](https://github.com/yufree/mzrtsim/workflows/Render%20and%20Deploy%20RMarkdown%20Website/badge.svg)](https://github.com/yufree/mzrtsim/actions)
<!-- badges: end -->

The goal of mzrtsim is to make batch effects simulation for LC/GC-MS
based peaks list data

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("yufree/mzrtsim")
```

## Batch effect simulation

You could use `mzrtsim` to make simulation.

``` r
library(mzrtsim)
data("monams1")
simdata <- mzrtsim(ncomp = 100, ncond = 2, ncpeaks = 0.05,
  nbatch = 3, nbpeaks = 0.1, npercond = 10, nperbatch = c(8, 5, 7), seed = 42, batchtype = 'mb', db=monams1)
```

Here we make a simulation of 100 compounds from selected database with
two conditions and three batches. 5 percentage of the peaks were
influnced by conditions and 10 percentage of the peaks were influnced by
batch effects. Three different type could be simulated: monotonic,
random and block. You could also bind batchtype, for example, ‘mb’ means
the simulation would contain both monotonic and block batch effects.
‘db’ means the spectra database to be used for simulation. ‘monams1’
means LC-MS spectra from MoNA amd ‘hmdbcms’ means GC-EI-MS spectra from
HMDB.

You could save the simulated data into multiple csv files by `simdata`
function. `simraw.csv` could be used for metaboanalyst. `simraw2.csv`
show the raw peaks list. `simcon.csv` show peaks influnced by conditions
only.`simbatchmatrix.csv` show peaks influnced by batch effects only.
`simbat.csv` show peaks influnced by batch effects and conditions.
`simcomp.csv` show independent peaks influnced by conditions and batch
effects. `simcompchange.csv` show the conditions changes of each groups.
`simblockbatchange.csv` show the block batch changes of each groups.
`simmonobatchange.csv` show the monotonic batch changes of each groups.

``` r
simdata(sim,name = "sim")
```

## Basic models of batch effects correction

Batch effects are the variances caused by factor other than the
experimental design. We could simply make a linear model for the
intensity of one peak:

\[Intensity =  Average + Condition + Batch + Error\]

Research is focused on condition contribution part and overall average
or random error could be estimated. However, we know little about the
batch contribution. Sometimes we could use known variables such as
injection order or operators as the batch part. However, in most cases
we such variable is unknown. Almost all the batch correction methods are
trying to use some estimations to balance or remove the batch effect.

For analytical chemistry, internal standards or pool quality control
samples are actually standing for the batch contribution part in the
model. However, it’s impractical to get all the internal standards when
the data is collected untargeted. For methods using internal standards
or pool quality control samples, the variations among those samples are
usually removed as median, quantile, mean or the ratios. Other ways like
quantile regression, centering and scaling based on distribution within
samples could be treated as using the stable distribution of peaks
intensity to remove batch effects. However, here we would turn to the
experiments without pooled quality control samples and correction
methods which make estimation of batch effects. Some combination methods
using pooled quality control samples or control samples to estimate the
batch effect like Remove Unwanted Variation(RUV) would not be covered
here.

### Distribution of intensity

Intensity collects from LC/GC-MS always showed a right-skewed
distribution. Log transformation is often necessary for further
statistical analysis. In some case, a Log-transformed intensity could be
visualized easily.

## Evaluation of batch correction

Various methods have been used for batch correction and evaluation. For
our simulation, we need some methods to check the performances of
correction. Difference analysis would be a common method for evaluation.
Then we could check whether this peak is true positive or false positive
by settings of the simulation.