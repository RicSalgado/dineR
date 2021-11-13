library(testthat)
test_that("estimation behaves as expected", {

  library(MASS)
  library(progress)
  library(Matrix)

  n_X <- 100
  n_Y <- n_X
  p_X <- 100
  p_Y <- p_X
  case <- "sparse"

  data <- data_generator(n = n_X, p = p_X, seed = 123)

  X <- data$X
  Y <- data$Y
  diff_Omega <- data$diff_Omega

  loss <- "lasso"
  nlambda <- 3
  tuning <- "none"
  stop_tol <- 1e-4
  perturb <- F
  correlation <- F
  max_iter <- 500
  lambda_min_ratio <- 0.5

  result <- estimation(X, Y, loss = loss, nlambda = nlambda, tuning = tuning, stop_tol = stop_tol,
                       perturb = perturb, correlation = correlation,
                       max_iter = max_iter, lambda_min_ratio = lambda_min_ratio)

  expect_equal(length(result$lambdas), nlambda)
  expect_equal(result$n_X, n_X)
  expect_equal(ncol(result$Sigma_X), p_X)

})
