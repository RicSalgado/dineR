# Data Generating Function

# The first components of this function are if statement to ensure the
# necessary dimensions are maintained.

#' Data Generator
#'
#' This functions generates two \eqn{n} by \eqn{p} size samples of multivariate normal
#' data. In doing this it also determines and provides the relevant covariance
#' matrices.
#'
#' @param n The number of observations generated.
#' @param p The number of dimensions for the generated samples.
#' @param Delta Optional parameter - Provides the differential network that will be used to obtain the sample covariance matrices.
#' @param case Optional parameter - Selects under which case the covariance matrices are determined. Possible cases are: "sparse" - Sparse Case or "asymsparse"- Asymptotically Sparse Case. Defaults to "sparse".
#' @param seed Optional parameter - Allows a seed to be set for reproducibility.
#'
#'
#' @return A list of various outputs, namely:
#' \itemize{
#' \item case - The case used.
#' \item seed_option - The seed provided.
#' \item X - The first multivariate normal sample.
#' \item Y - The second multivariate normal sample.
#' \item Sigma_X - The covariance matrix of X.
#' \item Sigma_Y - The covariance matrix of Y.
#' \item Omega_X - The precision matrix of X.
#' \item Omega_Y - The precision matrix of Y.
#' \item diff_Omega - The difference of precision matrices.
#' \item Delta - The target differential network.
#' }
#' @export
#'
#' @examples data <- data_generator(n = 100, p = 50, seed = 123)
#' @examples data <- data_generator(n = 10, p = 50, case = "asymsparse")

#' @import MASS
#' @importFrom "stats" "cor" "cov" "qnorm" "sd" "toeplitz"

data_generator <- function(n, p, Delta = NULL, case = "sparse", seed = NULL){

  if(n < 1){
    warning("The number of observations is too few.")
    return(NULL)
  }

  if(p < 2){
    warning("The number of dimensions is too few.")
    return(NULL)
  }

  cases <- c("sparse", "asymsparse")

  if(!is.element(case, cases)){
    warning("Please specify an appropriate case.")
    return(NULL)
  }

  if(length(seed) > 1){
    warning("Please provide an appropriate seed.")
    return(NULL)
  }

  results <- list() # This is a list where each of the elements we want to access after running the function

  if(is.null(Delta)){
    Delta <- matrix(0,p,p)
    Delta[1:2,1:2] <- matrix(c(0, -1, -1, 2), 2) #This is the matrix used by the source paper
  }

  if(nrow(Delta) != ncol(Delta)){
    warning("The provided differential network is not square.")
    return(NULL)
  }

  if(!isSymmetric(Delta)){
    warning("The provided differential network is not symmetric.")
    return(NULL)
  }

  # The function has randomness built in, thus a seed must be included

  if(is.null(seed)){
    seed_option <- "No seed was specified."
  }else{
    seed_option <- seed
  }

  set.seed(seed)

  # The individual precision matrices can now be generated

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

  # Using the precision matrices, the covariance matrices are
  # solved for and then used to generate the data

  diff_Omega <- Omega_Y - Omega_X

  Sigma_X <- solve(Omega_X)
  Sigma_Y <- solve(Omega_Y)

  # The above returns the covariance matrix, as it provides the inverse of the input
  # and the inverse of the precision matrix is just the covariance matrix

  # Using the above, normal data with a zero mean is generated

  X <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma_X)
  Y <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma_Y)

  results$case <- case
  results$seed <- seed_option
  results$X <- X
  results$Y <- Y
  results$Sigma_X <- Sigma_X
  results$Sigma_Y <- Sigma_Y
  results$Omega_X <- Omega_X
  results$Omega_Y <- Omega_Y
  results$diff_Omega <- diff_Omega
  results$Delta <- Delta

  return(results)

}
