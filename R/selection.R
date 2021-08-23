selection <- function(output, gamma , Sigma_X = NULL, Sigma_Y = NULL, tuning = c("none","AIC","BIC","EBIC")){

    output$tuning <- tuning[1]
    
    n_X <- output$n_X
    n_Y <- output$n_Y
    
    if(is.null(Sigma_X)){
        
        Sigma_X <- output$Sigma_X
        
    } 
    if(is.null(Sigma_Y)){
        
        Sigma_Y <- output$Sigma_Y
        
    }
    
    p <- ncol(Sigma_X)
    
    res <- switch(tuning[1], # Default tuning is none
                  none = NA,
                  AIC = diffnet_ic(Sigma_X, Sigma_Y, output$path, output$loss, n_X + n_Y, 2),
                  BIC = diffnet_ic(Sigma_X, Sigma_Y, output$path, output$loss, n_X + n_Y, log(n_X + n_Y)),
                  EBIC = diffnet_ic(Sigma_X, Sigma_Y, output$path, output$loss, n_X + n_Y, log(n_X + n_Y) + 2*gamma*log(p))   
                  )
    
    if(tuning[1] != "none"){
        output$opt <- res$opt
        output$ic <- res$ic
    }

    class(output) <- "diffnet"
    rm(n_X, n_Y)
    return(output)
}

diffnet_ic <- function(Sigma_X, Sigma_Y, path, type_of_loss, n, penalty){
    
    lowtri <- which(lower.tri(path[[1]], diag=TRUE))
    df <- sapply(path, function(x){ sum(x[lowtri]!=0); })
    
    if(type_of_loss != "d-trace"){
        
        ic <- scale(n*sapply(path, lossFunction, Sigma_X, Sigma_Y), center = -penalty*df, scale = FALSE)
    
    }else{
        
        ic <- scale(n*sapply(path, lossFunctionDtrace, Sigma_X, Sigma_Y), center = -penalty*df, scale = FALSE)
        
    }
    
    opt <- vector("list", 2)

    opt_ind <- apply(ic, 1, which.min)
    names(opt_ind) <- c("min","F")
    opt[[1]] <- opt_ind

    opt_loss <- apply(ic, 1, min)
    names(opt_loss) <- c("min","F")
    opt[[2]] <- opt_loss
    
    optimal <- list()
    optimal$opt <- opt
    optimal$ic <- ic
    return(optimal)
    
}

lossFunction <- function(Delta, Sigma_X, Sigma_Y){
    
    err <- Sigma_X%*%Delta%*%Sigma_Y - Sigma_X + Sigma_Y
    return(c(max(abs(err)), sqrt(sum(err^2)))) # Max & Frobenius norms
        
}

lossFunctionDtrace <- function(Delta, Sigma_X, Sigma_Y){
    
  err <- (Sigma_X%*%Delta%*%Sigma_Y + Sigma_Y%*%Delta%*%Sigma_X) / 2 - Sigma_X + Sigma_Y
  return(c(max(abs(err)), sqrt(sum(err^2)))) # Max & Frobenius norms 

}
