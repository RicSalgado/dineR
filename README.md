# dineR

<img src="https://raw.githubusercontent.com/RicSalgado/dineR/master/raw/sticker/dineR.png" alt="dineR logo" width="150" height="150" align="right"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/RicSalgado/diner/workflows/R-CMD-check/badge.svg)](https://github.com/RicSalgado/diner/actions)
[![CodeFactor](https://www.codefactor.io/repository/github/ricsalgado/diner/badge)](https://www.codefactor.io/repository/github/ricsalgado/diner)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Coveralls coverage](https://coveralls.io/repos/github/google/benchmark/badge.svg?branch=master)](https://coveralls.io/github/google/benchmark)
[![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/dineR?color=blue)](https://cran.r-project.org/package=dineR)

<!-- badges: end -->

## Overview

*dineR* is a R package, that aims to enable users of all backgrounds to
easily, and computationally efficiently perform differential network
estimation. Data can either be provided directly, or simulated. The
differential network is then efficiently estimated through the use of a
selected loss function and either sequential or parallel optimization.

## Installation

You can install the latest release version of dineR from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("dineR")
```

However, it is also possible to install the development version from
[GitHub](https://github.com/RicSalgado/dineR) with:

``` r
devtools::install_github("RicSalgado/dineR")
```

## Usage

A basic workflow of estimating a differential network using *dineR* is
as follows:

``` r
# Load dineR into the current R session
library(dineR)
```

### Data Generation

*dineR* allows for any qualifying data to be used to estimate a
differential network.

For the purposes of this demonstration, we will simulate data using
*dineR*'s in-house functions.

``` r
# The sample size of the first sample:
n_X <- 100

# The sample size of the second sample:
n_Y <- n_X

# The number of features within our first sample:
p_X <- 100

# The number of features within our second sample:
p_Y <- p_X

# The covariance matrix structure of interest:
case <- "sparse"

# Generating the data:
data <- data_generator(n = n_X, p = p_X, seed = 123)

# Extracting the relevent samples:
X <- data$X
Y <- data$Y
```

An advantage of using the above approach to simulate data, over other
options within in R is that the data is guaranteed to be of the form
required by *dineR*, but it also produces the inverse covariance
matrices and as such it is possible to evaluate the accuracy of the
estimation approach.

To evaluate the analytical differential network simply call the
following:

``` r
diff_Omega <- data$diff_Omega
```

### Estimation

To perform differential network estimation, there are a variety of
options provided to the user to fit every use case and scenario. Each of
these options comes with an appropriately selected default allowing
users to get up and running with their estimation with greater ease.

The most basic function call in which only the two samples are specified
would then be as follows:

``` r
results <- estimation(X, Y)
```

For further details regarding each of the available options, the reader
is encouraged to review the function documentation and accompanying
literature.

### Outputs

Having completed the differential network estimation, there are a number
of outputs available to view to confirm the estimation has behaved as
anticipated.

The output of greatest interest however is easily the estimated
differential network which can be accessed as follows:

``` r
First_Estimate <- result$path[[1]][1:5, 1:5]
```

The above code extracts the first 5 rows, and 5 columns of the
differential network estimate obtained for the first value of the tuning
parameter $\lambda$.

### Getting more involved

The above example is an extremely high-level introduction into the topic
of differential network estimation through the use of *dineR* and it is
recommended any users who wish to explore the topic further consult the
available documentation for the relevant details.

Happy estimating.
