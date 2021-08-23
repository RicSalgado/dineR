# Data Generating Function

# The first components of this function are if statement to ensure the 
# necessary dimensions are maintained.

data_generator = function(n, p, Delta = NULL, case = "sparse", seed = NULL){
  
  if(n < 1){
    cat("The number of observations is too few.")
    return(NULL)
  }
  
  if(p < 2){
    cat("The number of dimensions is too few.")
    return(NULL)
  }
  
  cases <- c("sparse", "asymsparse")
  
  if(!is.element(case, cases)){
    cat("Please specify an appropriate case.")
    return(NULL)
  }
  
  if(length(seed) > 1){
    cat("Please provide an appropriate seed.")
    return(NULL)
  }
  
  results <- list() # This is a list where each of the elements we want to access after running the function
  
  if(is.null(Delta)){
    Delta <- matrix(0,p,p)
    Delta[1:2,1:2] <- matrix(c(0, -1, -1, 2), 2) #This is the matrix used by the source paper
  } 
  
  if(nrow(Delta) != ncol(Delta)){
    cat("The provided differential network is not square.")
    return(NULL)
  }
  
  if(!isSymmetric(Delta)){
    cat("The provided differential network is not symmetric.")
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