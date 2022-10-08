#' Estimation
#'
#' This function performs alternating direction method of multipliers optimization
#' for a variety of loss functions to estimate the differential network given
#' two samples of multivariate normal data.
#'
#' @param X The first multivariate normal sample.
#' @param Y The second multivariate normal sample.
#' @param lambdas Optional parameter - A list of the regularization values to be used within the loss functions.
#' @param lambda_min_ratio Optional parameter - Defines the smallest regularization values as this proportion of the largest regularization value. Defaults to 0.3.
#' @param nlambda Optional parameter - The number of regularization values considered. Defaults to 10.
#' @param a Optional parameter - The thresholding parameter used in SCAD and MCP loss functions. Defaults to 3.7 with SCAD, and 3 with MCP respectively.
#' @param loss Optional parameter - The loss function of choice to implement. The function allows for four choices, namely "lasso", "scad", "mcp" and "d-trace". Defaults to "lasso".
#' @param tuning Optional parameter - The tuning method selected to determine the optimal value for the regularization parameter. Options are "none", "AIC", "BIC" and "EBIC". Defaults to "none".
#' @param perturb Optional parameter - When set to TRUE perturbation as done by the CLIME software to improve performance is implemented. Options are TRUE or FALSE, with the function defaulting to FALSE.
#' @param stop_tol Optional parameter - The stop tolerance to determine whether convergence has occurred. Defaults to 1e-5.
#' @param max_iter Optional parameter - The maximum number of iterations that can be perform for any one regularization value. Defaults to 100.
#' @param correlation Optional parameter - Determines whether the sample correlation matrices should be used in the place of the sample covariance matrices. Choices are TRUE and FALSE with the function defaulting to FALSE.
#' @param Delta_init Optional parameter - Allows for the algorithm to provided an initial estimate of the differential network to ease computation.
#' @param rho Optional parameter - Allows the user to adjust the ADMM step-size. Defaults to 1.
#' @param gamma Optional parameter - Allows the user to adjust the EBIC value when EBIC is the selected tuning method. Defaults to 0.5.
#' @param cores Optional parameter - Allows the user to specify the number of cores used by the optimization. Defaults to 1, i.e sequential solving however appropriate values range from 2 to the minimum of the number of lambdas and the total number of available cores.
#' @param verbose Optional parameter - Allows the user to obtain a summary of the estimation results. Options are TRUE or FALSE, where FALSE indicates the summary is not provided. Defaults to FALSE.
#'
#' @return A list of various outputs, namely:
#' \itemize{
#' \item n_X - The number of observations in X.
#' \item n_Y - The number of observations in Y.
#' \item Sigma_X - The covariance matrix of X.
#' \item Sigma_Y - The covariance matrix of Y.
#' \item loss - The loss function implemented.
#' \item tuning - The tuning method utilized.
#' \item lip - The value of the lipschitz constant.
#' \item iter - The iterations until convergence for each of the regularization values.
#' \item elapse - The total system time (in seconds) elapsed from initialization to completion of the optimization.
#' \item lambdas - The regularization parameter values used.
#' \item sparsity - The level of sparsity of the differential network for each regularization value.
#' \item path - The set of all differential networks for all regularization values considered.
#' \item ic - The output obtained from any possible tuning.
#' \item ic_index - The index at which the tuning is optimized.
#' \item ic_value - The tuning method optimal value.
#' \item chosen_lambda_ic - The regularization value that occurs at **ic_index**.
#' \item loss_index - The index at which the loss function is optimized.
#' \item loss_value - The loss function optimal value.
#' \item chosen_lambda_loss - The regularization value that occurs at **loss_index**.
#' }
#'
#' @export
#'
#' @examples data <- data_generator(n = 100, p = 50, seed = 123)
#' @examples X <- data$X
#' @examples Y <- data$Y
#' @examples result <- estimation(X,Y)
#'
#' @import doParallel
#' @import foreach
#' @import parallel
#' @import progress

estimation <- function(X, Y, lambdas = NULL, lambda_min_ratio = 0.3, nlambda = 10, a = NULL,
                       loss = "lasso", tuning = "none", perturb = FALSE, stop_tol = 1e-5,
                       max_iter = 500, correlation = FALSE, Delta_init = NULL, rho=NULL, gamma=NULL,
                       cores = 1, verbose = FALSE){

  # WARNING MESSAGES

  if(ncol(X) == ncol(Y)){
    p <- ncol(X)
  }else{
    warning("The dimensions of the matrices are inconsistent.")
    return(NULL)
  }

  if(is.null(gamma)){
    gamma <- 0.5
  }

  if(gamma > 1 || gamma < 0){
    warning("Please provide an appropriate value for gamma.")
    return(NULL)
  }

  if(!is.null(lambdas)){
    if(lengths(lambdas) > 1){
      warning("Please provide a 1-dimensional vector of lambda values.")
      return(NULL)
    }
  }

  if((lambda_min_ratio <= 0) || (length(lambda_min_ratio) > 1)){
    warning("Please provide a valid, single value for the smallest value of lambda.")
    return(NULL)
  }

  if((nlambda < 1) || (nlambda %% 1 != 0) || length(nlambda) > 1){
    warning("Please provide a valid number of lambda values.")
    return(NULL)
  }

  losses <- c("d-trace", "lasso", "mcp", "scad")

  if(!is.element(loss, losses)){
    warning("Please select an appropriate loss function.")
    return(NULL)
  }

  tuning_options <- c("none","AIC","BIC", "EBIC")

  if(!is.element(tuning, tuning_options)){
    warning("Please select an appropriate tuning procedure.")
    return(NULL)
  }

  if(is.null(tuning)){

    tuning <- "none"

  }

  perturb_options <- c(F, FALSE, T, TRUE)

  if(!is.element(perturb, perturb_options)){
    warning("Please select either TRUE or FALSE for whether to perturb.")
    return(NULL)
  }

  if((stop_tol < 0) || (length(stop_tol) > 1)){
    warning("Please provide a valid stop tolerance.")
    return(NULL)
  }

  if((max_iter < 1) || (length(max_iter) > 1)){
    warning("Please provide a valid maximum number of iterations.")
    return(NULL)
  }

  correlation_options <- c(F, FALSE, T, TRUE)

  if(!is.element(correlation, correlation_options)){
    warning("Please select either TRUE or FALSE for whether to use the correlation matrices.")
    return(NULL)
  }

  if(!is.null(Delta_init)){
    if(nrow(Delta_init) != p){
      warning("The provided data and differential network have inconsistent dimensions.")
      return(NULL)
    }
    if(nrow(Delta_init) != ncol(Delta_init)){
      warning("The provided differential network is not square.")
      return(NULL)
    }
    if(!isSymmetric(Delta_init)){
      warning("The provided differential network is not symmetric.")
      return(NULL)
    }
  }

  if((cores < 1) || (cores %% 1 != 0) || (cores > (detectCores())-1) || length(cores) > 1){
    warning("Please provide a valid number of cores.")
    return(NULL)
  }

  if(cores > nlambda){
    warning(paste("The number of cores you have specified is larger than the number of lambdas. Thus, only", cores, "core(s) are being used."))
  }

  #################################################################

  # HELPER FUNCTIONS

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

  L1_dts <- function(SigmaX, SigmaY, rho, lambda, Delta0 = NULL, Lambda0 = NULL,
                     Ux = NULL, Dx = NULL, Uy = NULL, Dy = NULL, C1 = NULL, C2 = NULL,
                     stop.tol = stop.tol, max.iter = 1e3){

    if(is.null(Delta0)) Delta0 <- solve(SigmaY+diag(nrow(SigmaY)))-solve(SigmaX+diag(nrow(SigmaX)))
    if(is.null(Lambda0)) Lambda0 <- matrix(1,nrow(SigmaX),ncol(SigmaX))

    p <- dim(Delta0)[1]
    k <- 0
    if(is.null(Ux) || is.null(Ux) || is.null(Uy) || is.null(Dy) || is.null(C1) || is.null(C2)){
      M1 <- eigen(SigmaX)
      M2 <- eigen(SigmaY)
      Ux <- M1$vectors
      Dx <- M1$values
      Uy <- M2$vectors
      Dy <- M2$values
      C1 <- matrix(0,p,p)
      C2 <- matrix(0,p,p)

      for (i in 1:p){
        for (j in 1:p){
          C1[i,j] <- 1/(Dy[j]*Dx[i]+2*rho)
          C2[i,j] <- 1/(Dy[i]*Dx[j]+2*rho)
        }
      }
    }

    Delta1 <- Delta0
    Delta2 <- Delta0
    Delta3 <- Delta0
    Lambda1 <- Lambda0
    Lambda2 <- Lambda0
    Lambda3 <- Lambda0

    l <- p^2
    iter <- 0
    admm <- symmetric_admm(Delta0, Delta3, Lambda0, Lambda0, SigmaX,
                           SigmaY, rho, C1, C2,
                           Ux, Uy, stop.tol, p, lambda, max.iter)

    Lambda3 <- matrix(admm$Lambda3, p ,p)
    Delta3 <- matrix(admm$Delta3, nrow=p, ncol=p)
    iter <- admm$iter

    result <- list()
    result$Delta <- Delta3
    result$Lambda3 <- Lambda3
    result$iter <- iter

    return(result)
  }

  #################################################################

  # PREAMBLE TO OPTIMIZATION

  # Create a vector to save the output in
  fit <- list()

  # Save the loss function used
  fit$loss <- loss[1]

  #Some preliminaries
  n_X <- nrow(X)
  n_Y <- nrow(Y)
  fit$n_X <- n_X
  fit$n_Y <- n_Y

  # This uses the correlation matrices instead of the covariance matrices
  if(correlation){
    fit$Sigma_X <- cor(X)
    fit$Sigma_Y <- cor(Y)

    X <- scale(X, center = T, scale = T) / sqrt(n_X-1)
    Y <- scale(Y, center = T, scale = T) / sqrt(n_Y-1)
  }else{ # This standardizes the data if we are not using the correlation matrices so scale does not affect the network
    fit$Sigma_X <- cov(X)*(1 - 1/n_X)
    fit$Sigma_Y <- cov(Y)*(1 - 1/n_Y)

    X <- scale(X, center = T, scale = F) / sqrt(n_X)
    Y <- scale(Y, center = T, scale = F) / sqrt(n_Y)
  }

  # This is the same perturbation as CLIME which was shown to improve performance for single precision matrices

  if(perturb){
    eigvals_X <- eigen(fit$Sigma_X, only.values=TRUE)$values
    eigvals_Y <- eigen(fit$Sigma_Y, only.values=TRUE)$values

    epsilon_X <- max(max(eigvals_X) - p*min(eigvals_X), 0) / (p-1)
    epsilon_Y <- max(max(eigvals_Y) - p*min(eigvals_Y), 0) / (p-1)

    fit$Sigma_X <- fit$Sigma_X + epsilon_X*diag(p)
    fit$Sigma_Y <- fit$Sigma_Y + epsilon_Y*diag(p)
  }else{
    epsilon_X <- 0
    epsilon_Y <- 0
  }

  # Calculate the lipschitz constant
  lip <- eigen(fit$Sigma_X, symmetric = T, only.values = T)$value[1]*eigen(fit$Sigma_Y, symmetric = T, only.values = T)$value[1]

  width <- getOption("width")
  separator <- strrep("-", width)
  message(paste0(separator))
  message("\n")

  fit$lip <- lip

  # If no vector of lambdas is given, they must be calculated

  if(!is.null(lambdas)){

    nlambda <- length(lambdas)

  }else{
    lambda_max <- max(abs(fit$Sigma_X - fit$Sigma_Y))
    lambda_min <- lambda_min_ratio*lambda_max
    lambdas <- exp(seq(log(lambda_max), log(lambda_min), length = nlambda))
  }

  lambdas <- sort(lambdas, decreasing = T)
  fit$lambdas <- lambdas

  if(is.null(a)) {
    if(loss[1] == "scad") a <- 3.7
    if(loss[1] == "mcp") a <- 3
  }

  # For each lambda ADMM will be performed

  fit$path <- vector("list", nlambda)
  fit$sparsity <- rep(0, nlambda)
  fit$iter <- rep(0, nlambda)

  # The initial guess of Delta is given, is used otherwise a zero matrix is used

  if(is.null(Delta_init)){

    Delta <- matrix(0, p, p)

  }
  else if(is.matrix(Delta_init)){

    Delta <- Delta_init
  }

  n_iter <- length(lambdas)

  pb <- progress_bar$new(format = "Lambda: :lambda [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = n_iter,
                         clear = FALSE,
                         width = 100)

  progress_lambda <- seq(1:n_iter)

  #################################################################

  # SEQUENTIAL OPTIMIZATION

  if(cores == 1){

    start <- proc.time()[3]
    for(i in 1:length(lambdas)){

      lambda <- lambdas[i] # Extract our chosen lambda

      iter <- 0
      Lambda <- matrix(lambda, p, p)

      if(loss[1] == "lasso"){
        out <- diffnet_lasso(fit$Sigma_X, fit$Sigma_Y, p,
                             Lambda, X, Y, n_X, n_Y,
                             epsilon_X, epsilon_Y,
                             lip, stop_tol, max_iter, Delta) # All of the parameters, including the selected lambda are given to the ADMM function

        Delta <- matrix(out$Delta, ncol=p)
        fit$iter[i] <- out$iter

      }
      if(loss[1] == "scad"){
        out <- diffnet_scad(fit$Sigma_X, fit$Sigma_Y, p,
                            lambda, a, X, Y, n_X, n_Y,
                            epsilon_X, epsilon_Y,
                            lip, stop_tol, max_iter, Delta)

        Delta <- matrix(out$Delta, ncol=p)
        fit$iter[i] <- out$iter

      }
      if(loss[1] == "mcp"){
        out <- diffnet_mcp(fit$Sigma_X, fit$Sigma_Y, p,
                           lambda, a, X, Y, n_X, n_Y,
                           epsilon_X, epsilon_Y,
                           lip, stop_tol, max_iter, Delta)

        Delta <- matrix(out$Delta, ncol=p)
        fit$iter[i] <- out$iter

      }
      if(loss[1] == "d-trace"){

        if(is.null(rho)){rho <- 1}

        M1 <- eigen(fit$Sigma_X)
        M2 <- eigen(fit$Sigma_Y)
        Ux <- M1$vectors
        Dx <- M1$values
        Uy <- M2$vectors
        Dy <- M2$values
        C1 <- matrix(0,p,p)
        C2 <- matrix(0,p,p)

        for (k in 1:p){
          for (j in 1:p){
            C1[k,j] <- 1/(Dy[j]*Dx[k]+2*rho)
            C2[k,j] <- 1/(Dy[k]*Dx[j]+2*rho)
          }
        }

        Delta0 <- solve(fit$Sigma_Y + diag(nrow(fit$Sigma_Y))) - solve(fit$Sigma_X + diag(nrow(fit$Sigma_X)))

        Lambda0 <- matrix(0, p, p)

        out <- L1_dts(fit$Sigma_X, fit$Sigma_Y, rho, lambda, Delta0, Lambda0, Ux, Dx, Uy, Dy, C1, C2, stop.tol = stop_tol, max.iter = max_iter)
        Delta <- matrix(out$Delta, ncol=p)
        fit$iter[i] <- out$iter

      }

      pb$tick(tokens = list(lambda = progress_lambda[i]))

      fit$path[[i]] <- Matrix::Matrix(Delta, sparse = T)
      fit$sparsity[i] <- sum(Delta != 0) / p / (p-1)

    }

    fit$elapse <-  proc.time()[3] - start

    # Now we have the ADMM results for each of the supplied lambdas, it is necessary to see which lambda results
    # in the best fit and select that as the final solution

    message(paste0(separator))

    if(max_iter %in% fit$iter){

      message("The ADMM did not converge for one or more lambdas.")

    }

    if(!is.null(tuning)){
      tuning_results <- selection(fit, gamma, tuning = tuning[1])
    }

    rm(X, Y, Delta)
    class(fit) <- "diffnet"

    output <- list()

    output$n_X <- fit$n_X
    output$n_Y <- fit$n_Y
    output$Sigma_X <- fit$Sigma_X
    output$Sigma_Y <- fit$Sigma_Y
    output$loss <- fit$loss
    output$tuning <- tuning
    output$lip <- fit$lip
    output$iter <- fit$iter
    output$elapse <- fit$elapse
    output$lambdas <- fit$lambdas
    output$sparsity <- fit$sparsity
    output$path <- fit$path
    output$ic <- tuning_results$ic

    tuning_output <- c("AIC", "BIC", "EBIC")

    if(is.element(tuning, tuning_output)){

      ic_index <- tuning_results$opt[[1]][[1]]
      ic_value <- tuning_results$opt[[2]][[1]]
      chosen_lambda_ic <- lambdas[ic_index]

      output$ic_index <- ic_index
      output$ic_value <- ic_value
      output$chosen_lambda_ic <- chosen_lambda_ic

      loss_index <- tuning_results$opt[[1]][[2]]
      loss_value <- tuning_results$opt[[2]][[2]]
      chosen_lambda_loss <- lambdas[loss_index]

      output$loss_index <- loss_index
      output$loss_value <- loss_value
      output$chosen_lambda_loss <- chosen_lambda_loss
    }

    if(verbose == TRUE){

      summary.estimation(output)

    }

    return(output)

  }

  #######################################################

  # PARALLEL OPTIMIZATION

  if(cores != 1){

    cluster <- makeCluster(cores)
    registerDoSNOW(cluster)

    progress <- function(n){
      pb$tick(tokens = list(lambda = progress_lambda[n]))
    }

    opts <- list(progress = progress)

    start <- proc.time()[3]

    parallel_output <- foreach(i = 1:length(lambdas), .combine = c, .options.snow = opts) %dopar% {

      iterations <- c()
      path <- c()

      lambda <- lambdas[i] # Extract our chosen lambda

      iter <- 0
      Lambda <- matrix(lambda, p, p)

      if(loss[1] == "lasso"){
        out <- diffnet_lasso(fit$Sigma_X, fit$Sigma_Y, p,
                             Lambda, X, Y, n_X, n_Y,
                             epsilon_X, epsilon_Y,
                             lip, stop_tol, max_iter, Delta) # All of the parameters, including the selected lambda are given to the ADMM function

        Delta <- matrix(out$Delta, ncol=p)
        iterations <- c(iterations, out$iter)

      }
      if(loss[1] == "scad"){
        out <- diffnet_scad(fit$Sigma_X, fit$Sigma_Y, p,
                            lambda, a, X, Y, n_X, n_Y,
                            epsilon_X, epsilon_Y,
                            lip, stop_tol, max_iter, Delta)

        Delta <- matrix(out$Delta, ncol=p)
        iterations <- c(iterations, out$iter)

      }
      if(loss[1] == "mcp"){
        out <- diffnet_mcp(fit$Sigma_X, fit$Sigma_Y, p,
                           lambda, a, X, Y, n_X, n_Y,
                           epsilon_X, epsilon_Y,
                           lip, stop_tol, max_iter, Delta)

        Delta <- matrix(out$Delta, ncol=p)
        iterations <- c(iterations, out$iter)

      }
      if(loss[1] == "d-trace"){

        if(is.null(rho)){rho <- 1}

        M1 <- eigen(fit$Sigma_X)
        M2 <- eigen(fit$Sigma_Y)
        Ux <- M1$vectors
        Dx <- M1$values
        Uy <- M2$vectors
        Dy <- M2$values
        C1 <- matrix(0,p,p)
        C2 <- matrix(0,p,p)

        for (k in 1:p){
          for (j in 1:p){
            C1[k,j] <- 1/(Dy[j]*Dx[k]+2*rho)
            C2[k,j] <- 1/(Dy[k]*Dx[j]+2*rho)
          }
        }

        Delta0 <- solve(fit$Sigma_Y + diag(nrow(fit$Sigma_Y))) - solve(fit$Sigma_X + diag(nrow(fit$Sigma_X)))

        Lambda0 <- matrix(0, p, p)

        out <- L1_dts(fit$Sigma_X, fit$Sigma_Y, rho, lambda, Delta0, Lambda0, Ux, Dx, Uy, Dy, C1, C2, stop.tol = stop_tol, max.iter = max_iter)
        Delta <- matrix(out$Delta, ncol=p)
        iterations <- c(iterations, out$iter)

      }

      path[[i]] <- Matrix::Matrix(Delta, sparse = T)
      #fit$path[[i]] <- Matrix::Matrix(Delta, sparse = T)
      #fit$sparsity[i] <- sum(Delta != 0) / p / (p-1)

      #return(fit)
      return(path)


    }

    #fit$iter <- parallel_output
    fit$path <- parallel_output

    fit$elapse <-  proc.time()[3] - start

    # Now we have the ADMM results for each of the supplied lambdas, it is necessary to see which lambda results
    # in the best fit and select that as the final solution

    message(paste0(separator))

    if(max_iter %in% fit$iter){

      message("The ADMM did not converge for one or more lambdas.")

    }

    if(!is.null(tuning)){
      tuning_results <- selection(fit, gamma, tuning = tuning[1])
    }

    rm(X, Y, Delta)
    class(fit) <- "diffnet"

    output <- list()

    output$n_X <- fit$n_X
    output$n_Y <- fit$n_Y
    output$Sigma_X <- fit$Sigma_X
    output$Sigma_Y <- fit$Sigma_Y
    output$loss <- fit$loss
    output$tuning <- tuning
    output$lip <- fit$lip
    output$iter <- fit$iter
    output$elapse <- fit$elapse
    output$lambdas <- fit$lambdas
    output$sparsity <- fit$sparsity
    output$path <- fit$path
    output$ic <- tuning_results$ic

    tuning_output <- c("AIC", "BIC", "EBIC")

    if(is.element(tuning, tuning_output)){

      ic_index <- tuning_results$opt[[1]][[1]]
      ic_value <- tuning_results$opt[[2]][[1]]
      chosen_lambda_ic <- lambdas[ic_index]

      output$ic_index <- ic_index
      output$ic_value <- ic_value
      output$chosen_lambda_ic <- chosen_lambda_ic

      loss_index <- tuning_results$opt[[1]][[2]]
      loss_value <- tuning_results$opt[[2]][[2]]
      chosen_lambda_loss <- lambdas[loss_index]

      output$loss_index <- loss_index
      output$loss_value <- loss_value
      output$chosen_lambda_loss <- chosen_lambda_loss
    }

    if(verbose == TRUE){

      summary.estimation(output)

    }

    stopCluster(cluster)
    return(output)

  }

}

summary.estimation <- function(x){

  width <- getOption("width")
  separator <- strrep("-", width)

  cat(paste0(separator))
  cat("\n")
  if(x$loss == "lasso"){
    cat("Model: Estimation was perform via Graphical LASSO.\n")
  }
  else if(x$loss == "mcp"){
    cat("Model: Estimation was perform via minimax concave penalty.\n")
  }
  else if(x$loss == "scad"){
    cat("Model: Estimation was perform via smoothly clipped absolute deviation penalty.\n")
  }else if(x$loss == "d-trace"){
    cat("Model: Estimation was perform via D-trace.\n")
  }

  cat("Number of lambdas:",length(x$lambda),"\n")
  cat("Graph dimension:",ncol(x$Sigma_X),"\n")
  cat("Range of lambda values considered:", min(x$lambdas), "----->", max(x$lambdas),"\n")
  cat("Levels of sparsity of Delta encountered:", min(x$sparsity), "----->", max(x$sparsity),"\n")
  cat(paste(separator))
  cat("\n")

  tuning_output <- c("AIC", "BIC", "EBIC")

  if(is.element(x$tuning, tuning_output)){
    cat("The results of the parameter selection are as follows:\n")
    width <- getOption("width")
    separator <- strrep("-", width)
    cat(paste0(separator))

    tab <- matrix(c(x$tuning, x$ic_index, x$ic_value, x$chosen_lambda_ic), ncol=1, byrow=TRUE)
    colnames(tab) <- c('Optimal Information Criteria')
    rownames(tab) <- c('Selected Information Criterion','Lambda Index','Information Criteria Value', 'Lambda Value')
    tab <- as.table(tab)

    print(tab)
    cat("\n")
    cat(paste0(separator))
    cat("\n")

    tab <- matrix(c(x$loss, x$loss_index, x$loss_value, x$chosen_lambda_loss), ncol=1, byrow=TRUE)
    colnames(tab) <- c('Minimum Loss Value')
    rownames(tab) <- c('Selected Loss Function','Lambda Index','Loss Function Value', 'Lambda Value')
    tab <- as.table(tab)

    print(tab)
    cat("\n")
    cat(paste0(separator))

  }

}
