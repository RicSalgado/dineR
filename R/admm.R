# This script contains "helper" functions, as well as the ADMM optimization

# This is the soft thresholding function that occurs within ADMM

soft <- function(X, Lambda){

  Y <- X
  Y <- (Y > Lambda)*(Y - Lambda) + (Y < (-Lambda)) * (Y + Lambda)
  return(Y)

}

# The loss function when not in a high-dimensional setting

small_p_loss_func <- function(Sigma_X, Sigma_Y, diff_Sigma, Delta){

  lf <- 0.5*sum(sum(Delta %*% (Sigma_X%*%Delta%*%Sigma_Y))) - sum(sum(Delta %*% diff_Sigma))
  return(lf)

}

# The loss function accounting for high-dimensional situations

big_p_loss_func <- function(X, Y, diff_Sigma, Delta){

  if(!isSymmetric(Delta)){
    warning("The Delta matrix provided is not symmetric.")
    return(NULL)
  }

  lf <- 0.5*sum(sum(Delta * t(X)%*%(X%*%Delta%*%t(Y)%*%Y))) - sum(sum(Delta * diff_Sigma))
  return(lf)

}

# The penalty function of the loss function

penalty_func <- function(Delta, Lambda){

  pl <- sum(sum(abs(Lambda * Delta)))
  return(pl)

}

# This is the actual ADMM algorithm

diffnet_lasso <- function(pSigma_X, pSigma_Y, p, pLambda, pX, pY, n_X, n_Y, epsilon_X, epsilon_Y,
                         lip, stop_tol, max_iter, pDelta){

  Sigma_X <- as.matrix(pSigma_X, p, p)
  Sigma_Y <- as.matrix(pSigma_Y, p, p)
  Lambda <- as.matrix(pLambda, p, p)

  diff_Sigma <- Sigma_X - Sigma_Y

  X <- as.matrix(pX, n_X, p)
  Y <- as.matrix(pY, n_Y, p)

  Delta <- as.matrix(pDelta, p, p)
  Delta_old <- Delta
  Delta_extra <- Delta

  grad_L <- matrix(NA, p, p)

  t <- 1
  t_old <- t

  err <- 0
  f_old <- 0
  f <- 0

  iter <- 0

  # Everything above is simply defining each of the components within the optimization algorithm

  if(n_X >= p || n_Y >= p){ # When the observations outnumber the variables i.e low dimensional setting

      f <- small_p_loss_func(Sigma_X, Sigma_Y, diff_Sigma, Delta) + penalty_func(Delta, Lambda)

    }else{ # If the above does not occur, then we need to use the adjusted loss function

      f <- big_p_loss_func(X, Y, diff_Sigma, Delta) + penalty_func(Delta, Lambda)

    }

  # We then perform the iterative procedure until convergence

  while((err > stop_tol) & (iter < max_iter) || iter == 0 ){

    Delta_extra <- Delta + (t_old -  1) / t * (Delta - Delta_old)
    Delta_old <- Delta

    # Compute the gradient
    if(n_X >= p || n_Y >= p){

      grad_L <- Sigma_X%*%Delta_extra%*%Sigma_Y
      grad_L <- (grad_L + t(grad_L)) / 2 - diff_Sigma

    }else{

      if(epsilon_X == 0 & epsilon_Y == 0){

        grad_L <- t(X) %*% (X %*% Delta_extra %*% t(Y)) %*% Y
        grad_L <- ( grad_L + t(grad_L) ) / 2 - diff_Sigma

      }else{

        grad_L <- t(X) %*% (X%*%Delta_extra%*%t(Y)) %*% Y  + epsilon_Y %*% t(X) %*% (X%*%Delta_extra) + epsilon_X %*% (Delta_extra%*%t(Y)) %*% Y + epsilon_X %*% epsilon_Y %*% Delta_extra
        grad_L <- ( grad_L + t(grad_L) ) / 2 - diff_Sigma

      }

    }

    # Apply the soft thresholding

    Delta <- soft(Delta_extra - grad_L/lip, Lambda/lip)

    f_old <- f

    if(n_X >= p || n_Y >= p){

      f <- small_p_loss_func(Sigma_X, Sigma_Y, diff_Sigma, Delta) + penalty_func(Delta, Lambda)

    }else{

      f <- big_p_loss_func(X, Y, diff_Sigma, Delta) + penalty_func(Delta, Lambda)

    }

    t_old <- t
    t <- (1 + sqrt(1 + 4* t^2)) / 2

    err <- (f_old - f)/(f_old + 1)
    err <- abs(err)

    iter <- iter + 1

  }

  output <- list()

  output$iter <- iter
  output$Delta <- Delta

  return(output)

}

compute_mcp_lambda <- function(Delta, lambda, a){

  weighted_lambda <- (Delta <= a*lambda) * (lambda - Delta/a)
  return(weighted_lambda)
}

diffnet_mcp <- function(pSigma_X, pSigma_Y, p,
                        lambda, a, pX, pY, n_X, n_Y,
                        epsilon_X, epsilon_Y,
                        lip, stop_tol, max_iter, pDelta){

  Delta <- as.matrix(pDelta, p, p)

  all_iter <- 0
  iter <- all_iter

  weighted_lambda <- compute_mcp_lambda(Delta, lambda, a)

  output <- diffnet_lasso(pSigma_X, pSigma_Y, p, weighted_lambda, pX, pY, n_X, n_Y,
                epsilon_X, epsilon_Y,
                lip, stop_tol, max_iter, pDelta)

  return(output)

}

BIG <- 1e20 # This is the upper bound for the scad penalty

compute_scad_lambda <- function(Delta, lambda, a){

 weighted_lambda <- lambda * ((Delta <= lambda) + (Delta > lambda) * (pmax(pmin(a*lambda - Delta, 0), BIG) / (a-1) / lambda))
 return(weighted_lambda)
}

diffnet_scad <- function(pSigma_X, pSigma_Y, p,
                        lambda, a, pX, pY, n_X, n_Y,
                        epsilon_X, epsilon_Y,
                        lip, stop_tol, max_iter, pDelta){

  Delta <- as.matrix(pDelta, p, p)
  all_iter <- 0
  iter <- all_iter

  weighted_lambda <- compute_scad_lambda(abs(Delta), lambda, a)

  output <- diffnet_lasso(pSigma_X, pSigma_Y, p, weighted_lambda, pX, pY, n_X, n_Y,
                          epsilon_X, epsilon_Y,
                          lip, stop_tol, max_iter, pDelta)

  return(output)

}

##############################################################################

loss_func <- function(Sigma_X, Sigma_Y, diff_Sigma, Delta){

  lf <- 0.25*sum(sum(Delta * (Sigma_X%*%Delta%*%Sigma_Y))) + 0.25*sum(sum(Delta * (Sigma_Y%*%Delta%*%Sigma_X))) - sum(sum(Delta * diff_Sigma))
  return(lf)

}

symmetric_admm <- function(pDelta0, pDelta3, pLambda0, pLambda3,
                           pSigmaX, pSigmaY, rho,
                           pC1, pC2, pUx, pUy, tol, p, lambda, max_iter){

  SigmaX <- as.matrix(pSigmaX, p, p)

  SigmaY <- as.matrix(pSigmaY, p, p)

  diff_Sigma <- SigmaX - SigmaY

  C1 <- as.matrix(pC1, p, p)
  C2 <- as.matrix(pC2, p, p)

  Ux <- as.matrix(pUx, p, p)
  Uy <- as.matrix(pUy, p, p)

  A <- matrix(NA, p, p)
  B <- matrix(NA, p, p)
  C <- matrix(NA, p ,p)

  Delta1 <- as.matrix(pDelta0, p, p)
  Delta2 <- as.matrix(pDelta0, p, p)
  Delta3 <- as.matrix(pDelta3, p, p)

  temp <- matrix(NA, p, p)

  Lambda1 <- as.matrix(pLambda0, p, p)

  Lambda2 <- as.matrix(pLambda0, p, p)

  Lambda3 <- as.matrix(pLambda3, p, p)

  iter <- 0
  f_old <- 0
  err <- 0

  f <- loss_func(SigmaX, SigmaY, diff_Sigma, Delta3) + penalty_func(Delta3, lambda)

  while(iter < max_iter){

    A <- SigmaX - SigmaY + rho*Delta2 + rho*Delta3 + Lambda3 - Lambda1

    Delta1 <- Uy %*% (C1 * (t(Uy)%*%A%*%Ux)) %*% t(Ux)

    B <- SigmaX - SigmaY + rho*Delta1 + rho*Delta3 + Lambda1 - Lambda2

    Delta2 <- Ux %*% (C2 *(t(Ux)%*%B%*%Uy)) %*% t(Uy)

    C <- (Lambda2/rho - Lambda3/rho + Delta1 + Delta2) / 2

    Delta3 <- soft(C, (lambda/rho) / 2)

    temp <- Delta3

    Delta3 <- (temp + t(temp))/2

    f_old <- f
    f <- loss_func(SigmaX, SigmaY, diff_Sigma, Delta3) + penalty_func(Delta3, lambda)

    err <- (f_old - f) / (f_old + 1)
    err <- abs(err)

    Lambda1 <- Lambda1 + rho * (Delta1 - Delta2)
    Lambda2 <- Lambda2 + rho * (Delta2 - Delta3)
    Lambda3 <- Lambda3 + rho * (Delta3 - Delta1)

    iter <- iter + 1

    if ((err < tol & iter != 0) || iter == max_iter){

      output <- list()
      output$Delta3 <- Delta3
      output$Lambda3 <- Lambda3
      output$iter <- iter

      return(output)

    }

  }

}
