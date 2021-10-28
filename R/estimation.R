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
#' @param lipschitz Optional parameter - Whether the lipshitz constant is printed or not. Choices are TRUE or FALSE, with the default being FALSE.
#' @param Delta_init Optional parameter - Allows for the algorithm to provided an initial estimate of the differential network to ease computation.
#' @param rho Optional parameter - Allows the user to adjust the ADMM step-size. Defaults to 1.
#' @param gamma Optional parameter - Allows the user to adjust the EBIC value when EBIC is the selected tuning method. Defaults to 0.5.
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
#' @examples estimation(X,Y)

#' @import progress

estimation <- function(X, Y, lambdas = NULL, lambda_min_ratio = 0.3, nlambda = 10, a = NULL,
                      loss = "lasso", tuning = "none", perturb = FALSE, stop_tol = 1e-5,
                      max_iter = 500, correlation = FALSE, lipschitz = FALSE, Delta_init = NULL, rho=NULL, gamma=NULL){

  # Warning messages
  if(ncol(X) == ncol(Y)){
    p <- ncol(X)
  }else{
    cat("The dimensions of the matrices are inconsistent.")
    return(NULL)
  }

  if(is.null(gamma)){
    gamma <- 0.5
  }

  if(gamma > 1 || gamma < 0){
    cat("Please provide an appropriate value for gamma.")
    return(NULL)
  }

  if(!is.null(lambdas)){
    if(lengths(lambdas) > 1){
      cat("Please provide a 1-dimensional vector of lambda values.")
      return(NULL)
    }
  }

  if((lambda_min_ratio <= 0) || (length(lambda_min_ratio) > 1)){
    cat("Please provide a valid, single value for the smallest value of lambda.")
    return(NULL)
  }

  if((nlambda < 1) || (nlambda %% 1 != 0) || length(nlambda) > 1){
    cat("Please provide a valid number of lambda values.")
    return(NULL)
  }

  losses <- c("d-trace", "lasso", "mcp", "scad")

  if(!is.element(loss, losses)){
    paste(c("Please select an appropriate loss function from the following list: ", losses))
    return(NULL)
  }

  tuning_options <- c("none","AIC","BIC", "EBIC")

  if(!is.element(tuning, tuning_options)){
    cat("Please select an appropriate tuning procedure.")
    paste(c("Please select an appropriate tuning procedure from the following list: ", tuning_options))
    return(NULL)
  }

  if(is.null(tuning)){

    tuning <- "none"

  }

  perturb_options <- c(F, FALSE, T, TRUE)

  if(!is.element(perturb, perturb_options)){
    cat("Please select either TRUE or FALSE for whether to perturb.")
    return(NULL)
  }

  if((stop_tol < 0) || (length(stop_tol) > 1)){
    cat("Please provide a valid stop tolerance.")
    return(NULL)
  }

  if((max_iter < 1) || (length(max_iter) > 1)){
    cat("Please provide a valid maximum number of iterations.")
    return(NULL)
  }

  correlation_options <- c(F, FALSE, T, TRUE)

  if(!is.element(correlation, correlation_options)){
    cat("Please select either TRUE or FALSE for whether to use the correlation matrices.")
    return(NULL)
  }

  lip_options <- c(F, FALSE, T, TRUE)

  if(!is.element(lipschitz, lip_options)){
    cat("Please select either TRUE or FALSE for whether to print the lipschitz constant.")
    return(NULL)
  }

  if(!is.null(Delta_init)){
    if(nrow(Delta_init) != p){
      cat("The provided data and differential network have inconsistent dimensions.")
      return(NULL)
    }
    if(nrow(Delta_init) != ncol(Delta_init)){
      cat("The provided differential network is not square.")
      return(NULL)
    }
    if(!isSymmetric(Delta_init)){
      cat("The provided differential network is not symmetric.")
      return(NULL)
    }
  }

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
  cat(paste0(separator))
  cat("\n")

  if(lipschitz){

    cat("The Lipschitz constant is: ", round(lip,3), "\n")

    }

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

  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = n_iter,
                         complete = "=",
                         incomplete = "-",
                         current = ">",
                         clear = FALSE,
                         width = 100)

  # This is entire estimation algorithm begins

  start <- proc.time()[3]
  for(i in 1:length(lambdas)){

    pb$tick()

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
          C1[k,j] = 1/(Dy[j]*Dx[k]+2*rho)
          C2[k,j] = 1/(Dy[k]*Dx[j]+2*rho)
        }
      }

      Delta0 <- solve(fit$Sigma_Y + diag(nrow(fit$Sigma_Y))) - solve(fit$Sigma_X + diag(nrow(fit$Sigma_X)))

      Lambda0 <- matrix(0, p, p)

      out <- L1_dts(fit$Sigma_X, fit$Sigma_Y, rho, lambda, Delta0, Lambda0, Ux, Dx, Uy, Dy, C1, C2, stop.tol = stop_tol, max.iter = max_iter)
      Delta <- matrix(out$Delta, ncol=p)
      fit$iter[i] <- out$iter

    }

    fit$path[[i]] <- Matrix::Matrix(Delta, sparse = T)
    fit$sparsity[i] <- sum(Delta != 0) / p / (p-1)

  }

  fit$elapse <-  proc.time()[3] - start

  # Now we have the ADMM results for each of the supplied lambdas, it is necessary to see which lambda results
  # in the best fit and select that as the final solution

  cat(paste0(separator))

  if(max_iter %in% fit$iter){

    cat("\n")
    cat("The ADMM did not converge for one or more lambdas.")
    cat("\n")

  }

  if(!is.null(tuning)){
    tuning_results <- selection(fit, gamma, tuning = tuning[1])
  }

  rm(X, Y, Delta)
  class(fit) <- "diffnet"

  print.diffnet <- function(x, ...)
  {

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
  }

  print.diffnet(fit)

  output <- list()

  output$n_X <- fit$n_X
  output$n_Y <- fit$n_Y
  output$Sigma_X <- fit$Sigma_X
  output$Sigma_Y <- fit$Sigma_Y
  output$loss <- fit$loss
  output$tuning <- fit$tuning
  output$lip <- fit$lip
  output$iter <- fit$iter
  output$elapse <- fit$elapse
  output$lambdas <- fit$lambdas
  output$sparsity <- fit$sparsity
  output$path <- fit$path
  output$ic <- tuning_results$ic

  tuning_output <- c("AIC", "BIC", "EBIC")

  if(is.element(tuning, tuning_output)){
    cat("The results of the parameter selection are as follows:\n")
    width <- getOption("width")
    separator <- strrep("-", width)
    cat(paste0(separator))

    ic_index <- tuning_results$opt[[1]][[1]]
    ic_value <- tuning_results$opt[[2]][[1]]
    chosen_lambda_ic <- lambdas[ic_index]

    tab <- matrix(c(tuning, ic_index, ic_value, chosen_lambda_ic), ncol=1, byrow=TRUE)
    colnames(tab) <- c('Optimal Information Criteria')
    rownames(tab) <- c('Selected Information Criterion','Lambda Index','Information Criteria Value', 'Lambda Value')
    tab <- as.table(tab)

    print(tab)
    cat("\n")
    cat(paste0(separator))
    cat("\n")

    output$ic_index <- ic_index
    output$ic_value <- ic_value
    output$chosen_lambda_ic <- chosen_lambda_ic

    loss_index <- tuning_results$opt[[1]][[2]]
    loss_value <- tuning_results$opt[[2]][[2]]
    chosen_lambda_loss <- lambdas[loss_index]

    tab <- matrix(c(loss, loss_index, loss_value, chosen_lambda_loss), ncol=1, byrow=TRUE)
    colnames(tab) <- c('Minimum Loss Value')
    rownames(tab) <- c('Selected Loss Function','Lambda Index','Loss Function Value', 'Lambda Value')
    tab <- as.table(tab)

    print(tab)
    cat("\n")
    cat(paste0(separator))

    output$loss_index <- loss_index
    output$loss_value <- loss_value
    output$chosen_lambda_loss <- chosen_lambda_loss
  }

  return(output)
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
