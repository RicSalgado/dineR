library(testthat)
test_that("data_generator produces expected data", {

  library(MASS)

  n_X <- 100
  n_Y <- n_X
  p_X <- 100
  p_Y <- p_X
  case <- "sparse"

  data <- data_generator(n = n_X, p = p_X, seed = 123)

  X <- data$X
  Y <- data$Y

  expect_equal(length(X), 10000)

})
