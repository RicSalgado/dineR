
# dineR

<img src="https://raw.githubusercontent.com/RicSalgado/dineR/master/raw/sticker/dineR.png" width="150" height="150" align="right"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/RicSalgado/diner/workflows/R-CMD-check/badge.svg)](https://github.com/RicSalgado/diner/actions)
[![CodeFactor](https://www.codefactor.io/repository/github/ricsalgado/diner/badge)](https://www.codefactor.io/repository/github/ricsalgado/diner)
[![Lifecycle](https://github.com/RicSalgado/dineR-dev/tree/dev/man/figures/lifecycle-stable.svg)
[![](https://coveralls.io/repos/github/google/benchmark/badge.svg?branch=master)](https://coveralls.io/github/google/benchmark)
[![](http://cranlogs.r-pkg.org/badges/grand-total/dineR?color=blue)](https://cran.r-project.org/package=dineR)
<!-- badges: end -->

## Overview

*dineR* is a R package, that aims to enable users of all backgrounds to easily, 
and computationally efficiently perform differential network estimation. Data can 
either be provided directly, or simulated. The differential network is then efficiently
estimated through the use of a selected loss function and either sequential or 
parallel optimization. 

## Installation

You can install the latest release version of dineR from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("dineR")
```

However, it is also possible to install the development version from [GitHub](https://github.com/RicSalgado/dineR) with: 

``` r
devtools::install_github("RicSalgado/dineR")
```

## Usage

A basic workflow of estimating a differential network using *dineR* is as follows:

``` r
# Load dineR into the current R session
library(dineR)
```


``` r
# Data Generation
n_X <- 100
n_Y <- n_X
p_X <- 100
p_Y <- p_X
#case <- "sparse"
case <- "asymsparse"

data <- data_generator(n = n_X, p = p_X, seed = 123)

X <- data$X
Y <- data$Y
diff_Omega <- data$diff_Omega
paste("The number of non-zero entries in the differential network is: ", sum(diff_Omega!=0))

# Estimation Preliminaries (All of the parameters are now optional as the function has pre-specified defaults)

loss <- "lasso"
nlambda <- 50
tuning <- "AIC"
stop_tol <- 1e-4
perturb <- F
correlation <- F
max_iter <- 500
lambda_min_ratio <- 0.5
#gamma <- 1 #Only if we use EBIC

# Estimation

result <- estimation(X, Y, loss = loss, nlambda = nlambda, tuning = tuning, stop_tol = stop_tol,
                      perturb = perturb, correlation = correlation,
                      max_iter = max_iter, lambda_min_ratio = lambda_min_ratio) 

# Results

print(result$path[[1]][1:5,1:5])
result$elapse
```

