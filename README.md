
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mzrtsim

<!-- badges: start -->

[![R-CMD-check](https://github.com/yufree/mzrtsim/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yufree/mzrtsim/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

The goal of mzrtsim is to make raw data and features table simulation
for LC/GC-MS based data

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("yufree/mzrtsim")
```

## Raw Data simulation

MS1 full scan data has been proved more complex than theoretical
prediction. Recently study showed that soft ionization will also contain
[fragment ions](https://www.nature.com/articles/s41596-023-00803-0) for
structure identification and contain lots of [redundant
peaks](https://pubs.acs.org/doi/10.1021/acs.analchem.7b02380). As shown
in the following figure, MS1 full scan data could contain different
types of ions from the same compound (red). Some peaks we could figure
out the sources while some peaks might be hard to interpret. In this
case, mzrtsim package will use experimental data instead of predicted
peaks from the compound formula to show the complexity of MS1 spectra.

<figure>
<img src="https://yufree.github.io/presentation/figure/peakcom.png"
alt="demo" />
<figcaption aria-hidden="true">demo</figcaption>
</figure>

You could use `simmzml` to generate one mzML file.

``` r
library(mzrtsim)
data("monams1")
simmzml(db=monams1, name = 'test')
```

You will find `test.mzML` and corresponding `test.csv` with m/z,
retention time and compound name of the peaks. Here the `monams1` and
`monahrms1` is from the MS1 data of MassBank of North America (MoNA) and
could be downloaded from their
[website](https://mona.fiehnlab.ucdavis.edu/downloads). You could also
use `hmdbcms` to simulate EI source data extracted from
[HMDB](https://hmdb.ca/downloads). Here we only use the MS1 full scan
data for simulation.

### Data extraction

For HMDB, you need to download their “GC-MS Spectra Files (XML) -
Experimental” data [here](https://hmdb.ca/downloads). Then you need the
following R code to construct the database for this package.

``` r
xmlpath <- list.files('hmdb_experimental_cms_spectra/',full.names = T,recursive = T)
library(xml2)
getcms <- function(i){
  tree = xml2::read_xml(xmlpath[i])
  accession = xml_text(xml2::xml_find_all(tree,'/c-ms/database-id'))
  rti = as.numeric(xml_text(xml2::xml_find_all(tree,'/c-ms/retention-index')))
  prec = as.numeric(xml_text(xml2::xml_find_all(tree,'/c-ms/derivative-exact-mass')))
  formula = xml_text(xml2::xml_find_all(tree,'/c-ms/derivative-formula'))
  ins = xml_text(xml2::xml_find_all(tree,'/c-ms/notes'))
  ins2 = xml_text(xml2::xml_find_all(tree,'/c-ms/instrument-type'))
  mode = xml_text(xml2::xml_find_all(tree,'/c-ms/ionization-mode'))
  idms = xml_text(xml2::xml_find_all(tree,'/c-ms/id'))
  instr = xml_text(xml2::xml_find_all(tree,'instrument-type'))
  masscharge = as.numeric(xml_text(xml2::xml_find_all(tree,'/c-ms/c-ms-peaks/c-ms-peak/mass-charge')))
  intensity =  as.numeric(xml_text(xml2::xml_find_all(tree,'/c-ms/c-ms-peaks/c-ms-peak/intensity')))
  idx <- intensity!=0
  mz <- list(name=accession, idms=idms, ionmode=mode, prec = prec, formula = formula, np = length(masscharge[idx]), rti = rti, instr = ins2, msm = ins, spectra=cbind.data.frame(mz=masscharge[idx],ins=intensity[idx]))
  return(mz)
}

library(parallel)
hmdbcms <- mcMap(getcms,1:length(xmlpath), mc.cores = 4)
saveRDS(hmdbcms,'hmdbcms.RDS')
```

For MoNA, you can download the “LC-MS Spectra” data in format of “msp”
from this [website](https://mona.fiehnlab.ucdavis.edu/downloads). Then
we need to filter the MS1 data to generate database for LC-MS and subset
data base for LC-HRMS.

``` r
# you need to install enviGCMS package by 'install.packages("enviGCMS")'.
msp <- enviGCMS::getMSP('MoNA-export-LC-MS_Spectra.msp')
idx <- sapply(msp, function(x) grepl('MS1',x$msm))
idx1 <- sapply(idx, function(x) length(x)>0)
monams1 <- msp[idx1]
idx2 <- sapply(monams1, function(x) grepl('MS1',x$msm))
monams1 <- monams1[idx2]

idx3 <- sapply(monams1, function(x) x$instr )
idx4 <- sapply(idx3, function(x) length(x)>0)
monams12 <- monams1[idx4]
saveRDS(monams12,'monams1.RDS')

idx5 <- sapply(monams12, function(x) x$instr )
monams13 <- monams12[grepl('FT|TOF',idx5)]
idx5 <- sapply(monams13, function(x) x$instr )

nn <- sapply(monams13,function(x) x$np)
median(as.numeric(nn))
mean(as.numeric(nn))
saveRDS(monams13,'monahrms1.RDS')
```

If you have your own database in “msp” format, you can generate the
custom database by the following code:

``` r
msp <- enviGCMS::getMSP('MoNA-export-LC-MS_Spectra.msp')
saveRDS(msp,'custom.RDS')
```

To use the “\*.RDS” database, you need to load them and refer their name
for ‘simmzml’ function.

``` r
custom <- readRDS('custom.RDS')
simmzml(db=custom, name = 'test')
```

To select certain compound for review/simulation, you could use name or
inchikey.

``` r
data(monahrms1)
# if you know the compound name
compoundname <- 'Triacylglycerol 14:0-16:0-16:0'
databasename <- sapply(monahrms1,function(x) x$name)
compoundidx <- databasename%in%compoundname
compound <- monahrms1[compoundidx]
compound
#> $`50051`
#> $`50051`$name
#> [1] "Triacylglycerol 14:0-16:0-16:0"
#> 
#> $`50051`$charge
#> [1] ""
#> 
#> $`50051`$ionmode
#> [1] "P"
#> 
#> $`50051`$prec
#> [1] 796.7
#> 
#> $`50051`$formula
#> [1] "C49H94O6"
#> 
#> $`50051`$exactmass
#> [1] 778.705
#> 
#> $`50051`$inchikey
#> [1] "JWVXCFSNEOMSHH-UHFFFAOYSA-N"
#> 
#> $`50051`$np
#> [1] "77"
#> 
#> $`50051`$ce
#> [1] "30 V"
#> 
#> $`50051`$rt
#> [1] 3729.8
#> 
#> $`50051`$column
#> [1] "ACQUITY UPLC BEH C18 column, Waters"
#> 
#> $`50051`$instr
#> [1] "LC-ESI-QTOF"
#> 
#> $`50051`$msm
#> [1] "MS1"
#> 
#> $`50051`$spectra
#>          mz        ins
#> 1  211.1904   1.550388
#> 2  221.2071   0.775194
#> 3  222.2356   0.387597
#> 4  239.2106   2.325581
#> 5  240.3326   1.550388
#> 6  245.4224   0.387597
#> 7  256.2933   0.387597
#> 8  257.2285   0.387597
#> 9  263.2312   0.387597
#> 10 264.1962   0.387597
#> 11 265.2154   0.387597
#> 12 285.2063   1.162791
#> 13 313.2623   0.775194
#> 14 488.2370   0.387597
#> 15 495.4385   2.325581
#> 16 496.4580   1.550388
#> 17 497.3890   0.775194
#> 18 497.7592   0.387597
#> 19 501.4811   0.387597
#> 20 517.1404   0.775194
#> 21 519.2001   0.775194
#> 22 519.8795   1.162791
#> 23 520.6265   2.325581
#> 24 522.2178   0.775194
#> 25 523.4456 100.000000
#> 26 523.7699   1.550388
#> 27 524.4470  41.860465
#> 28 525.4677   3.875969
#> 29 526.4229   0.775194
#> 30 527.4040   4.263566
#> 31 528.4025   0.387597
#> 32 529.1721   0.775194
#> 33 530.2812   0.775194
#> 34 530.9286   0.387597
#> 35 531.8790   0.387597
#> 36 532.3360   0.775194
#> 37 536.6188   0.387597
#> 38 538.9778   0.387597
#> 39 539.6244   0.387597
#> 40 541.8032   0.387597
#> 41 543.8115   0.387597
#> 42 545.5234   0.387597
#> 43 547.5812   1.937984
#> 44 548.4410   0.387597
#> 45 548.8799   0.387597
#> 46 549.5575   0.775194
#> 47 551.4700  38.759690
#> 48 552.5022  15.503876
#> 49 553.4861   6.976744
#> 50 554.8922   0.775194
#> 51 555.3589   0.387597
#> 52 555.7817   0.775194
#> 53 556.7309   1.162791
#> 54 557.2531   0.387597
#> 55 558.0118   0.387597
#> 56 559.4038   0.387597
#> 57 561.1909   0.387597
#> 58 561.7617   0.387597
#> 59 562.1678   0.387597
#> 60 562.6122   0.387597
#> 61 564.3788   0.387597
#> 62 565.3078   0.387597
#> 63 569.1461   0.387597
#> 64 572.5359   0.387597
#> 65 573.4073   0.387597
#> 66 575.5505   0.775194
#> 67 576.4756   0.775194
#> 68 577.4467   0.775194
#> 69 578.4893   0.775194
#> 70 579.5178   1.162791
#> 71 580.4246   1.162791
#> 72 595.8934   0.387597
#> 73 599.5848   0.387597
#> 74 601.4676   0.775194
#> 75 608.1764   0.387597
#> 76 661.2155   0.387597
#> 77 704.8252   0.387597
#> 
#> 
#> $`50067`
#> $`50067`$name
#> [1] "Triacylglycerol 14:0-16:0-16:0"
#> 
#> $`50067`$charge
#> [1] ""
#> 
#> $`50067`$ionmode
#> [1] "P"
#> 
#> $`50067`$prec
#> [1] 796.7
#> 
#> $`50067`$formula
#> [1] "C49H94O6"
#> 
#> $`50067`$exactmass
#> [1] 778.705
#> 
#> $`50067`$inchikey
#> [1] "JWVXCFSNEOMSHH-UHFFFAOYSA-N"
#> 
#> $`50067`$np
#> [1] "11"
#> 
#> $`50067`$ce
#> [1] "30 V"
#> 
#> $`50067`$rt
#> [1] 3727.02
#> 
#> $`50067`$column
#> [1] "ACQUITY UPLC BEH C18 column, Waters"
#> 
#> $`50067`$instr
#> [1] "LC-ESI-QTOF"
#> 
#> $`50067`$msm
#> [1] "MS1"
#> 
#> $`50067`$spectra
#>          mz   ins
#> 1  245.2631  12.5
#> 2  263.2137  12.5
#> 3  523.4880 100.0
#> 4  524.4532  75.0
#> 5  525.4467  25.0
#> 6  526.6809  12.5
#> 7  537.4563  25.0
#> 8  537.8970  12.5
#> 9  551.4921  37.5
#> 10 552.5740  12.5
#> 11 575.5162  37.5

# if you know the inchikey
compoundinchikey <- c('KQSFNXMDCOFFGW-GNDIVNLPSA-N','WRYJYFCCMSVEPQ-ORAXXRKOSA-N')
databaseinchikey <- sapply(monahrms1,function(x) x$inchikey)
compoundidx <- databaseinchikey%in%compoundinchikey
compound <- monahrms1[compoundidx]
compound
#> $`51812`
#> $`51812`$name
#> [1] "Araloside A"
#> 
#> $`51812`$charge
#> [1] ""
#> 
#> $`51812`$ionmode
#> [1] "N"
#> 
#> $`51812`$prec
#> [1] NA
#> 
#> $`51812`$formula
#> [1] "C47H74O18"
#> 
#> $`51812`$exactmass
#> [1] 926.4875
#> 
#> $`51812`$inchikey
#> [1] "KQSFNXMDCOFFGW-GNDIVNLPSA-N"
#> 
#> $`51812`$np
#> [1] "34"
#> 
#> $`51812`$ce
#> [1] ""
#> 
#> $`51812`$rt
#> [1] 912.347
#> 
#> $`51812`$column
#> [1] "Waters Atlantis T3 (2.1 x 150 mm, 5 um)"
#> 
#> $`51812`$instr
#> [1] "LC-ESI-ITTOF"
#> 
#> $`51812`$msm
#> [1] "MS1"
#> 
#> $`51812`$spectra
#>           mz        ins
#> 1   215.2013   0.976541
#> 2   485.2377   4.109053
#> 3   485.7246   0.976541
#> 4   486.2331   0.774378
#> 5   492.2482  13.645963
#> 6   492.7386   7.021813
#> 7   493.2507   2.910515
#> 8   493.7417   1.174814
#> 9   508.2165   0.774378
#> 10  925.4726 100.000000
#> 11  925.9696  82.191145
#> 12  926.2328   1.491368
#> 13  926.4668  64.071599
#> 14  926.9934  21.107286
#> 15  927.4908  10.436486
#> 16  927.9884   2.291380
#> 17  928.4569   1.493936
#> 18  982.4831   2.493279
#> 19  982.9651   1.161569
#> 20  988.4561   1.812136
#> 21 1039.4615   2.945574
#> 22 1040.4841   1.914932
#> 23 1041.4763   0.774378
#> 24 1388.7144  16.350562
#> 25 1389.2157  28.941044
#> 26 1389.7172  21.982023
#> 27 1390.2188  15.618221
#> 28 1390.7563   8.083492
#> 29 1391.2222   3.260675
#> 30 1391.7957   0.740115
#> 31 1852.9233   1.678967
#> 32 1853.4610   0.774378
#> 33 1854.3713   0.637322
#> 34 1973.4138   0.740115
#> 
#> 
#> $`51940`
#> $`51940`$name
#> [1] "Saikosaponin B2"
#> 
#> $`51940`$charge
#> [1] ""
#> 
#> $`51940`$ionmode
#> [1] "N"
#> 
#> $`51940`$prec
#> [1] NA
#> 
#> $`51940`$formula
#> [1] "C42H68O13"
#> 
#> $`51940`$exactmass
#> [1] 780.466
#> 
#> $`51940`$inchikey
#> [1] "WRYJYFCCMSVEPQ-ORAXXRKOSA-N"
#> 
#> $`51940`$np
#> [1] "3"
#> 
#> $`51940`$ce
#> [1] ""
#> 
#> $`51940`$rt
#> [1] NA
#> 
#> $`51940`$column
#> [1] "Waters Atlantis T3 (2.1 x 150 mm, 5 um)"
#> 
#> $`51940`$instr
#> [1] "LC-ESI-ITTOF"
#> 
#> $`51940`$msm
#> [1] "MS1"
#> 
#> $`51940`$spectra
#>         mz      ins
#> 1 112.9847  23.9250
#> 2 779.4593 100.0000
#> 3 780.4455  35.2949
```

### Multiple files with experiment design

You could stimulate two groups of raw data with different peak heights
for the same compounds. Retention time will follow a uniform
distribution. 100 compounds could be selected randomly and base peaks’
signal to noise ratio could be sampled from 100 to 1000. When signal to
noise ratio is fixed for all samples, the peaks height will still be the
only changed parameters for both groups. Each group contain 10 samples
and 30% compounds are changed between case and control groups.

``` r
dir.create('case')
dir.create('control')
# set different peak height for 100 compounds
ph1 <- c(rep(5,30),rep(10,40),rep(15,30))
ph2 <- c(rep(5,20),rep(10,30),rep(15,50))
# set retention time for 100 compounds. When the peak width is set as a single value, the simulated base peak width will be from a poisson distribution with this value as lambda. 
rt <- seq(10,590,length.out=100)
set.seed(1)
# select compounds from database
compound <- sample(c(1:4000),100)
set.seed(2)
# select signal to noise ratio for a larger dynamic range
sn <- sample(c(100:10000),100)
for(i in c(1:10)){
  simmzml(name=paste0('case/case',i),db=monahrms1,pheight = ph1,compound=compound,rtime = rt, sn=sn)
}

for(i in c(1:10)){
  simmzml(name=paste0('control/control',i),db=monahrms1,pheight = ph2,compound=compound,rtime = rt, sn=sn)
}
```

Then you could find 10 mzML files in case sub folder and another 10 mzML
files in control sub folder, as well as corresponding csv files with
m/z, retention time and compound name of the peaks.

If you prefer to simulate retention time shifts for different samples,
you could add a random number to the retention time of each compound.

``` r
rt0 <- rt
set.seed(42)
for(i in c(1:10)){
        rt <- rt0 + rnorm(100) 
        simmzml(name=paste0('case/case',i),db=monahrms1,pwidth = pw1,compound=compound,rtime = rt, sn=sn)
}
```

### Chromatography peaks

You could also use `simmzml` to stimulate tailing/leading peaks by
defining the tailing factor of the peaks. When the tailing factor is
lower than 1, the peaks are leading peaks. When the tailing factor is
larger than 1, the peaks are tailing peaks.

``` r
# leading peaks
simmzml(name='test',db=monahrms1,pwidth = 10,compound=1,rtime = 100, sn=10,tailingfactor = 0.8)
# tailing peaks
simmzml(name='test',db=monahrms1,pwidth = 10,compound=1,rtime = 100, sn=10,tailingfactor = 1.5)
```

### matrix stimulation

You could also input a m/z vector as matrix masses. Those masses will
generate background baseline signals. By default, the mass vector is
from matrix samples previous
[published](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-022-00586-8).

``` r
data(mzm)
simmzml(name='test',db=monahrms1,pwidth = 10,compound=1,rtime = 100, sn=10,matrixmz = mzm,matrix = TRUE)
```

### Retention time simulation

We included column types and experimental retention time (in seconds)
for LC-MS data and retention index for GC-MS data. You could use the
following code to access such information.

``` r
data(monahrms1)
# select compound by their index
compound <- c(1,2,3)
column <- sapply(monahrms1[compound],function(x) x$column)
# print column
column
#>                                 50050                                 50051 
#> "ACQUITY UPLC BEH C18 column, Waters" "ACQUITY UPLC BEH C18 column, Waters" 
#>                                 50052 
#> "ACQUITY UPLC BEH C18 column, Waters"
retentiontime <- sapply(monahrms1[compound],function(x) x$rt)
# print retention time
retentiontime
#>   50050   50051   50052 
#> 3969.40 3729.80 3712.13
data("hmdbcms")
# select compound by their index
compound <- c(1,2,300)
retentionindex <- sapply(hmdbcms[compound],function(x) x$rti)
# print retention index
retentionindex
#> [1] 1807.71 2048.73 2638.99
```

## Peaks list simulation

You could also use `mzrtsim` to make simulation of peak list.

Here we make a simulation of 100 compounds from selected database with
two conditions and three batches. 5 percentage of the peaks were
influenced by conditions and 10 percentage of the peaks were influenced
by batch effects. Three different type could be simulated: monotonic,
random and block. You could also bind batch type, for example, ‘mb’
means the simulation would contain both monotonic and block batch
effects. ‘db’ means the spectra database to be used for simulation as
metioned in raw data simulation section.

``` r
library(mzrtsim)
data("monams1")
simdata <- mzrtsim(ncomp = 100, ncond = 2, ncpeaks = 0.05,
  nbatch = 3, nbpeaks = 0.1, npercond = 10, nperbatch = c(8, 5, 7), seed = 42, batchtype = 'mb', db=monams1)
```

You could save the simulated data into multiple csv files by `simdata`
function. `simraw.csv` could be used for metaboanalyst. `simraw2.csv`
show the raw peaks list. `simcon.csv` show peaks influenced by
conditions only.`simbatchmatrix.csv` show peaks influnced by batch
effects only. `simbat.csv` show peaks influenced by batch effects and
conditions. `simcomp.csv` show independent peaks influenced by
conditions and batch effects. `simcompchange.csv` show the conditions
changes of each groups. `simblockbatchange.csv` show the block batch
changes of each groups. `simmonobatchange.csv` show the monotonic batch
changes of each group.

``` r
simdata(sim,name = "sim")
```
