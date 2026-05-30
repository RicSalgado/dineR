# Data Generating Function

This function generates two multivariate normal samples, by means of
simulation. The samples can be represented as follows: \$\$X \sim N_p(0,
\Sigma_1)\$\$ \$\$Y \sim N_p(0, \Sigma_2)\$\$ with the only restriction
being that each sample has the same number of features \\p\\.

## Usage

``` r
data_generator(
  n_X = NULL,
  n_Y = NULL,
  p = NULL,
  Delta = NULL,
  case = "sparse",
  seed = NULL
)
```

## Arguments

- n_X:

  The number of observations to be generated for the first sample.

- n_Y:

  The number of observations to be generated for the second sample.

  - Only one of the above sample sizes need to be specified. In such a
    scenario, the sample size for the other sample is set to the same as
    the specified sample size.

- p:

  The dimensions/features for the samples.

- Delta:

  Optional parameter - Provides the differential network from which the
  sample covariance matrices must be derived.

- case:

  Optional parameter - Allows for the specification of the precision
  matrix structure. Possible cases are: "sparse" - Sparse Case or
  "asymsparse"- Asymptotically Sparse Case. Defaults to "sparse".

  - Sparse Case: \\\Omega_1 = (0.5^{\|i-j\|})^{-1}\\. That is,
    \\\\\Omega_1\\\_{1,1} = \\\Omega_1\\\_{p,p} = \frac{4}{3},
    \\\Omega_1\\\_{i,i} = \frac{5}{3}\\ for all other \\i\\.
    \\\\\Omega_1\\\_{i, {i + 1}} = \\\Omega_1\\\_{{i-1},i} =
    \frac{2}{3}\\ and \\\\\Omega_1\\\_{i, j} = 0\\ for all other \\i,
    j\\.

  - Asymptotically Sparse Case

- seed:

  Optional parameter - Allows a seed to be set for reproducibility.

## Value

A list of the various outputs, namely:

- case - The case used.

- seed_option - The seed used for simulation.

- X - The first multivariate normal sample.

- Y - The second multivariate normal sample.

- n_X - The number of observations simulated for X.

- n_Y - The number of observations simulated for Y.

- Sigma_X - The covariance matrix of X: \\\Sigma_X\\.

- Sigma_Y - The covariance matrix of Y: \\\Sigma_Y\\.

- Omega_X - The precision matrix of X: \\\Sigma_X^{-1} = \Omega_X\\.

- Omega_Y - The precision matrix of Y: \\\Sigma_Y^{-1} = \Omega_Y\\.

- Diff_Omega - The difference of the precision matrices: \\\Omega_X -
  \Omega_Y\\.

- Delta - The target differential network: \\\Delta\\.

## References

Tang, Z., Yu, Z. and Wang, C., 2020. A fast iterative algorithm for
high-dimensional differential network. Computational Statistics, 35(1),
pp.95-109.

## Examples

``` r
data <- data_generator(n_X = 100, p = 50, seed = 123)
data <- data_generator(n_X = 10, p = 50, case = "asymsparse")
```
