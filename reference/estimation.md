# Estimation

This function performs alternating direction method of multipliers
optimization for a variety of loss functions to estimate the
differential network given two samples of multivariate normal data.

## Usage

``` r
estimation(
  X,
  Y,
  lambdas = NULL,
  lambda_min_ratio = 0.3,
  nlambda = 10,
  a = NULL,
  loss = "lasso",
  tuning = "none",
  perturb = FALSE,
  stop_tol = 1e-05,
  max_iter = 500,
  correlation = FALSE,
  Delta_init = NULL,
  rho = NULL,
  gamma = NULL,
  cores = 1,
  verbose = FALSE
)
```

## Arguments

- X:

  The first multivariate normal sample.

- Y:

  The second multivariate normal sample.

- lambdas:

  Optional parameter - A list of the regularization values to be used
  within the loss functions.

- lambda_min_ratio:

  Optional parameter - Defines the smallest regularization values as
  this proportion of the largest regularization value. Defaults to 0.3.

- nlambda:

  Optional parameter - The number of regularization values considered.
  Defaults to 10.

- a:

  Optional parameter - The thresholding parameter used in SCAD and MCP
  loss functions. Defaults to 3.7 with SCAD, and 3 with MCP
  respectively.

- loss:

  Optional parameter - The loss function of choice to implement. The
  function allows for four choices, namely "lasso", "scad", "mcp" and
  "d-trace". Defaults to "lasso".

- tuning:

  Optional parameter - The tuning method selected to determine the
  optimal value for the regularization parameter. Options are "none",
  "AIC", "BIC" and "EBIC". Defaults to "none".

- perturb:

  Optional parameter - When set to TRUE perturbation as done by the
  CLIME software to improve performance is implemented. Options are TRUE
  or FALSE, with the function defaulting to FALSE.

- stop_tol:

  Optional parameter - The stop tolerance to determine whether
  convergence has occurred. Defaults to 1e-5.

- max_iter:

  Optional parameter - The maximum number of iterations that can be
  perform for any one regularization value. Defaults to 100.

- correlation:

  Optional parameter - Determines whether the sample correlation
  matrices should be used in the place of the sample covariance
  matrices. Choices are TRUE and FALSE with the function defaulting to
  FALSE.

- Delta_init:

  Optional parameter - Allows for the algorithm to provided an initial
  estimate of the differential network to ease computation.

- rho:

  Optional parameter - Allows the user to adjust the ADMM step-size.
  Defaults to 1.

- gamma:

  Optional parameter - Allows the user to adjust the EBIC value when
  EBIC is the selected tuning method. Defaults to 0.5.

- cores:

  Optional parameter - Allows the user to specify the number of cores
  used by the optimization. Defaults to 1, i.e sequential solving
  however appropriate values range from 2 to the minimum of the number
  of lambdas and the total number of available cores.

- verbose:

  Optional parameter - Allows the user to obtain a summary of the
  estimation results. Options are TRUE or FALSE, where FALSE indicates
  the summary is not provided. Defaults to FALSE.

## Value

A list of various outputs, namely:

- n_X - The number of observations in X.

- n_Y - The number of observations in Y.

- Sigma_X - The covariance matrix of X.

- Sigma_Y - The covariance matrix of Y.

- loss - The loss function implemented.

- tuning - The tuning method utilized.

- lip - The value of the lipschitz constant.

- iter - The iterations until convergence for each of the regularization
  values.

- elapse - The total system time (in seconds) elapsed from
  initialization to completion of the optimization.

- lambdas - The regularization parameter values used.

- sparsity - The level of sparsity of the differential network for each
  regularization value.

- path - The set of all differential networks for all regularization
  values considered.

- ic - The output obtained from any possible tuning.

- ic_index - The index at which the tuning is optimized.

- ic_value - The tuning method optimal value.

- chosen_lambda_ic - The regularization value that occurs at
  **ic_index**.

- loss_index - The index at which the loss function is optimized.

- loss_value - The loss function optimal value.

- chosen_lambda_loss - The regularization value that occurs at
  **loss_index**.

## References

Boyd, S., Parikh, N., Chu, E., Peleato, B. and Eckstein, J., 2011.
Distributed optimization and statistical learning via the alternating
direction method of multipliers. Foundations and Trends® in Machine
learning, 3(1), pp.1-122.

Chen, J. and Chen, Z., 2008. Extended Bayesian information criteria for
model selection with large model spaces. Biometrika, 95(3), pp.759-771.

Friedman, J., Hastie, T. and Tibshirani, R., 2008. Sparse inverse
covariance estimation with the graphical lasso. Biostatistics, 9(3),
pp.432-441.

Tang, Z., Yu, Z. and Wang, C., 2020. A fast iterative algorithm for
high-dimensional differential network. Computational Statistics, 35(1),
pp.95-109.

Yuan, H., Xi, R., Chen, C. and Deng, M., 2017. Differential network
analysis via lasso penalized D-trace loss. Biometrika, 104(4),
pp.755-770.

Zhang, T. and Zou, H., 2014. Sparse precision matrix estimation via
lasso penalized D-trace loss. Biometrika, 101(1), pp.103-120.

## Examples

``` r
data <- data_generator(n_X = 100, p = 50, seed = 123)
X <- data$X
Y <- data$Y

# Sequential (default)
result_seq <- estimation(X, Y, nlambda = 5, cores = 1)
#> --------------------------------------------------------------------------------
#> 
#> --------------------------------------------------------------------------------
result_seq$elapse
#> elapsed 
#>   0.727 

# \donttest{
# Parallel - set cores to a value greater than 1
# (detectCores() - 1 is recommended on your own machine)
result_par <- estimation(X, Y, nlambda = 5, cores = 2)
#> --------------------------------------------------------------------------------
#> 
#> --------------------------------------------------------------------------------
result_par$elapse
#> elapsed 
#>   1.153 
# }
```
