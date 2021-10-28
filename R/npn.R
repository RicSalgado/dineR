#' NPN - Non paranormal Transformation
#'
#' This functions allows us to transform non-normal multivariate data to that of
#' non paranormal data.
#'
#' @param x The multivariate non-normal data to be transformed.
#' @param npn_func Optional parameter - The method of transformation to be applied. Can either be "shrinkage" or "truncation" but defaults to "shrinkage".
#' @param npn_thresh Optional parameter - The truncation threshold that is used when making use of truncation.
#' @param verbose Optional parameter - Prints additional output of the selected approach. Can either be "TRUE" or "FALSE" and defaults to "TRUE".
#'
#' @return Returns the transformed data matrix.
#' @export
#'
#' @examples data <- data_generator(n = 100, p = 50, seed = 123)
#' @examples X <- data$X
#' @examples X_transformed <- npn(X, npn_func = "truncation")

npn <- function(x, npn_func = "shrinkage", npn_thresh = NULL, verbose = TRUE){

  n <- nrow(x)
  d <- ncol(x)
  x_col <- colnames(x)
  x_row <- rownames(x)

  # Shrinkage transformation
	if(npn_func == "shrinkage"){
		if(verbose){

		  cat("Conducting the nonparanormal (npn) transformation via shrunkun ECDF...\n")

		}

		x <- qnorm(apply(x, 2, rank) / (n + 1))
		x <- x / sd(x[, 1])

		if(verbose){

		  cat("done.\n")

		}

		rm(n, d, verbose)
   	colnames(x) <- x_col
		rownames(x) <- x_row
	}

	# Truncation transformation
	if(npn_func == "truncation"){
		if(verbose){

		  cat("Nonparanomral transformation conducted via truncated ECDF...\n")

		}
	  if(is.null(npn_thresh)){

		  npn_thresh <- 1 / (4*(n^0.25)*sqrt(pi*log(n)))
		}

	  x <- qnorm(pmin(pmax(apply(x, 2, rank) / n, npn_thresh), 1 - npn_thresh))
    x <- x / sd(x[, 1])

    if(verbose){

      cat("done.\n")
    }
    rm(n, d, npn_thresh, verbose)

   	colnames(x) <- x_col
		rownames(x) <- x_row
	}

	if(npn_func == "skeptic"){
		if(verbose){

		  cat("Nonparanormal transformation conducted via skeptic....")

		}
		x <- 2*sin(pi / 6*cor(x, method="spearman"))

		if(verbose){

		  cat("done.\n")
		}
		rm(n, d, verbose)

   	colnames(x) <- x_col
		rownames(x) <- x_col
	}

	return(x)
}
