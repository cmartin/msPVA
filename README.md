# msPVA : An R implementation of count-based multi-site population viability analysis

This package implements a stochastic simulation for a count-based multi-site population viability analysis as described in chapter 11 of Quantitative Conservation Ecology (Morris & Doak, 2002).

Code is highly inspired from the MatLab implementation described in the book.

## Installation
You have two options here.

### Install package from github.com
This allows you to easily install updates, have access to the function help files etc.

```r
library(devtools)
devtools::install_github("cmartin/msPVA")
library(msPVA)
```

### Or you download the source file and run it directly
In this case, you need download [the raw R file](https://raw.githubusercontent.com/cmartin/msPVA/master/R/MultisitePVA.R)
and then execute : 

```r
source("MultisitePVA.R")
```
This is a quicker, but dirtier way!

## Try some examples

With the same Clapper rail data as in Morris & Doak

### Definining parameters manually : 

```r
library(msPVA)
res <- simulate_ms_pva(
 leaving_prob = 0.2,
 reaching_prob = 0.5,
 growth_rate_means = c(0.043, -0.002, 0),
 growth_rate_vars = c(0.051, 0.041, 0.051),
 initial_pops = c(70, 26, 33),
 growth_rate_corrs = {matrix(
   c(	1.000,	0.995,   0.896,
   0.995,	1.000,   0.938,
   0.896,	0.938,   1.000),
   nrow = 3,
   ncol = 3,
   byrow = TRUE
 )},
 K = c(286, 60, 58),
 quasi_extinction_thresholds = c(20, 20, 20),
 n_years = 100,
 n_runs = 1000,
 with_progress_bar = FALSE 
)

print(res)
```

```
This is the deterministic lambda value :  1.0319573564599
And this is the mean stochastic lambda :  0.998326502312961
Below is mean and standard deviation of log lambda :
-0.00170197441398923 0.00737528068961553
```

```r
hist(res)
```

![](README_files/figure-html/unnamed-chunk-3-1.png) 

```r
plot(res)
```

![](README_files/figure-html/unnamed-chunk-3-2.png) 

### Or ask the package to calculate most of them automatically from a time-series : 
With a two-populations polar bear time series from [Stirling et al. 2004](http://arctic.journalhosting.ucalgary.ca/arctic/index.php/arctic/article/view/479/509)

```r
library(msPVA)
params <- calculate_params_from_file(
 system.file("extdata", "PolarBear_Stirling2004.csv", package = "msPVA")
)
```

```
Warning in calculate_params_from_file(system.file("extdata", "PolarBear_Stirling2004.csv", : There were gaps in the time series.
  Mean and variance were computed by linear regression.
  Correlation matrix was calculated assuming regular intervals.
```

```r
res <- do.call(
 "simulate_ms_pva",
 c(
   params,
   list(
     K = c(300,200),
     leaving_prob = 0.1,
     reaching_prob = 0.7,
     quasi_extinction_thresholds = c(20, 20),
     n_years = 100,
     n_runs = 1000,
     with_progress_bar = FALSE 
   )
 )
)
print(res)
```

```
This is the deterministic lambda value :  1.05786392886293
And this is the mean stochastic lambda :  0.994891450195734
Below is mean and standard deviation of log lambda :
-0.00513613229713802 0.00539279087662938
```

```r
hist(res)
```

![](README_files/figure-html/unnamed-chunk-4-1.png) 

```r
plot(res)
```

![](README_files/figure-html/unnamed-chunk-4-2.png) 

## Run it for a single population

```r
library(msPVA)
res <- simulate_ss_pva(
  growth_rate_means = 0.043,
  growth_rate_vars = 0.051,
  initial_pops = 70,
  K = 286,
  quasi_extinction_thresholds = 20,
  n_years = 50,
  n_runs = 100,
  with_progress_bar = FALSE 
)

print(res)
```

```
Over a 50 years span, the extinction risk of this population is 0.1
And the risk of decline is 0.17
```

```r
hist(res)
```

![](README_files/figure-html/unnamed-chunk-5-1.png) 

## Read the book
If you are to use this package, I highly recommend that you first read chapter 11 from [Quantitative Conservation Ecology (Morris & Doak, 2002)](http://www.sinauer.com/quantitative-conservation-biology-theory-and-practice-of-population-viability-analysis.html), so you understand limitations and assumptions from the underlying model. Managing animal populations should not be taken lightly.

## Problems
Please report any bugs to the [GitHub issue tracker](https://github.com/cmartin/msPVA/issues) and write any questions to <charles.martin1@uqtr.ca>

## Citation
If this code is useful to you, please cite as : 

```
Charles A. Martin (2015). msPVA: An R implementation of count-based multi-site population viability analysis. R package version 0.0.0.9002. https://github.com/cmartin/msPVA
```
