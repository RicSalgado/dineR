#' Data Generating Function
#'
#' This function generates two multivariate normal samples, by means of simulation.
#' The samples can be represented as follows:
#' \deqn{X \sim N_p(0, \Sigma_1)}
#' \deqn{Y \sim N_p(0, \Sigma_2)}
#' with the only restriction being that each sample has the same
#' number of features \eqn{p}.
#'
#' @param n_X The number of observations to be generated for the first sample.
#' @param n_Y The number of observations to be generated for the second sample.
#' \itemize{
#' \item Only one of the above sample sizes need to be specified. In such a scenario, the sample size for the other sample is set to the same as the specified sample size.}
#' @param p The dimensions/features for the samples.
#' @param Delta Optional parameter - Provides the differential network from which the sample covariance matrices must be derived.
#' @param case Optional parameter - Allows for the specification of the precision matrix structure. Possible cases are: "sparse" - Sparse Case or "asymsparse"- Asymptotically Sparse Case. Defaults to "sparse".
#' \itemize{
#' \item Sparse Case: \eqn{\Omega_1 = (0.5^{|i-j|})^{-1}}. That is, \eqn{\{\Omega_1\}_{1,1} = \{\Omega_1\}_{p,p} = \frac{4}{3}, \{\Omega_1\}_{i,i} = \frac{5}{3}} for all other \eqn{i}. \eqn{\{\Omega_1\}_{i, {i + 1}} = \{\Omega_1\}_{{i-1},i} = \frac{2}{3}} and \eqn{\{\Omega_1\}_{i, j} = 0} for all other \eqn{i, j}.
#' \item Asymptotically Sparse Case}
#' @param seed Optional parameter - Allows a seed to be set for reproducibility.
#'
#' @return A list of the various outputs, namely:
#' \itemize{
#' \item case - The case used.
#' \item seed_option - The seed used for simulation.
#' \item X - The first multivariate normal sample.
#' \item Y - The second multivariate normal sample.
#' \item n_X - The number of observations simulated for X.
#' \item n_Y - The number of observations simulated for Y.
#' \item Sigma_X - The covariance matrix of X: \eqn{\Sigma_X}.
#' \item Sigma_Y - The covariance matrix of Y: \eqn{\Sigma_Y}.
#' \item Omega_X - The precision matrix of X: \eqn{\Sigma_X^{-1} = \Omega_X}.
#' \item Omega_Y - The precision matrix of Y: \eqn{\Sigma_Y^{-1} = \Omega_Y}.
#' \item Diff_Omega - The difference of the precision matrices: \eqn{\Omega_X - \Omega_Y}.
#' \item Delta - The target differential network: \eqn{\Delta}.
#' }
#' @export
#'
#' @examples data <- data_generator(n_X = 100, p = 50, seed = 123)
#' @examples data <- data_generator(n_X = 10, p = 50, case = "asymsparse")
#'
#' @import MASS
#' @importFrom "stats" "cor" "cov" "qnorm" "sd" "toeplitz"
#'
#' @references Tang, Z., Yu, Z. and Wang, C., 2020. A fast iterative algorithm for high-dimensional differential network. Computational Statistics, 35(1), pp.95-109.

data_generator <- function(n_X = NULL, n_Y = NULL, p = NULL, Delta = NULL, case = "sparse", seed = NULL){

  # The first components of the function, are checks to ensure the validity of the supplied function arguments.

  # Sample Size

  if(is.null(n_X) & is.null(n_Y)){
    warning("The number of observations needs to be specified for at least one sample.")
    return(NULL)
  }else if(is.null(n_X)){
    n_X <- n_Y
  }else if(is.null(n_Y)){
    n_Y <- n_X
  }

  if(n_X < 1 | n_Y < 1){
    warning("The number of observations specified is too few.")
    return(NULL)
  }

  # Dimensions/Features

  if(is.null(p)){
    warning("The number of dimensions/features needs to be specified.")
    return(NULL)
  }else if(p < 2){
    warning("The number of dimensions specified is too few.")
    return(NULL)
  }

  # Differential Network

  if(is.null(Delta)){
    Delta <- matrix(0,p,p)
    Delta[1:2,1:2] <- matrix(c(0, -1, -1, 2), 2) #T his is the matrix used by Tang et al.
  }else if(nrow(Delta) != ncol(Delta)){
    warning("The provided differential network is not square.")
    return(NULL)
  }else if(!isSymmetric(Delta)){
    warning("The provided differential network is not symmetric.")
    return(NULL)
  }

  # Case

  cases <- c("sparse", "asymsparse")

  if(!is.element(case, cases)){
    warning("Please specify an appropriate case.")
    return(NULL)
  }

  # Seed

  if(length(seed) > 1){
    warning("Please provide a single appropriate seed.")
    return(NULL)
  }

  # Function

  results <- list() # Initialize a list to save each of the elements we want to access after function execution.

  # The function has randomness built in, thus a seed must be included.

  if(is.null(seed)){
    seed_option <- "No seed was specified."
  }else{
    seed_option <- seed
  }

  set.seed(seed)

  # The individual precision matrices can now be generated.

  if(case == 'sparse'){

    Omega_X <- matrix(0,p,p)

    ind <- row(Omega_X) - col(Omega_X) == 1
    Omega_X[ind] <- Omega_X[ind] + 2/3

    ind <- col(Omega_X) - row(Omega_X) == 1
    Omega_X[ind] <- Omega_X[ind] + 2/3

    diag(Omega_X) <- diag(Omega_X) + 5/3

    Omega_X[1,1] <- 4/3
    Omega_X[p,p] <- 4/3

  }else if(case == 'asymsparse'){ # Asymptotically Sparse Case

    Omega_X <- toeplitz(0.5^(0:(p-1)))

  }

  Omega_Y <- Delta + Omega_X

  # Using the above precision matrices, the covariance matrices are
  # solved for and then used to generate the data.

  Diff_Omega <- Omega_Y - Omega_X

  Sigma_X <- solve(Omega_X)
  Sigma_Y <- solve(Omega_Y)

  # The above returns the covariance matrix, as it provides the inverse of the input
  # and the inverse of the precision matrix is just the covariance matrix.

  # Using the above, normal data with a zero mean is generated.

  X <- mvrnorm(n = n_X, mu = rep(0, p), Sigma = Sigma_X)
  Y <- mvrnorm(n = n_Y, mu = rep(0, p), Sigma = Sigma_Y)

  results$case <- case
  results$seed <- seed_option
  results$X <- X
  results$Y <- Y
  results$n_X <- n_X
  results$n_Y <- n_Y
  results$Sigma_X <- Sigma_X
  results$Sigma_Y <- Sigma_Y
  results$Omega_X <- Omega_X
  results$Omega_Y <- Omega_Y
  results$Diff_Omega <- Diff_Omega
  results$Delta <- Delta

  return(results)

}
