#' Test for the existence of a finite moment (finite mean)
#'
#' This function tests the existence of the first finite moment (finite mean) in a sample. 
#' The null hypothesis \eqn{H_0} is that the first moment is finite, and the alternative hypothesis \eqn{H_1} 
#' is that the moment is not finite (it could be either infinite or undefined).
#' 
#' The test is based on the method proposed by \insertCite{fedotenkov2013bootstrap;textual}{TIE}.
#' 
#' @param x A numeric vector representing the sample. Only positive values are supported.
#' @param xi A numeric value specifying the threshold for the test. Default: 0.999, see \insertCite{fedotenkov2013bootstrap;textual}{TIE}.
#' @param nn An integer specifying the size of bootstrap subsamples. Default: \eqn{\max(0.4 \cdot \log(n), 2)}.
#' @param B An integer specifying the number of bootstrap subsamples. Default: 1E5.
#' @return A numeric value representing the p-value for the test.
#' @details 
#' The function uses a bootstrap approach to estimate the existence of a finite moment. 
#' By default, the number of bootstrap samples \eqn{B} is set to \code{1E5}, and the size of each 
#' bootstrap subsample is determined as the maximum of \eqn{0.4 \cdot \log(n)} or 2, where \eqn{n} 
#' is the size of the input sample.
#' 
#' If testing for a higher-order moment (e.g., second, third, or \eqn{k}-th), 
#' raise the input sample values to the power of 2, 3, or \eqn{k} respectively before calling the function.
#' 
#' @examples 
#' # Example: Testing the first finite moment with 
#' # a sample of absolute Cauchy random variables
#' n = 1000
#' x = abs(stats::rcauchy(n))
#' boottail(x)
#' @importFrom Rdpack reprompt
#' @references \insertRef{fedotenkov2013bootstrap}{TIE}
#' @export
boottail <- function(x, xi = 0.999, nn = NULL, B = 1E5) {
  x <- as.vector(x)  # Make a vector, if x is in another format
  n <- length(x)     # Number of observations
  
  # Determine nn if not provided
  if (is.null(nn)) {
    nn <- max(round(0.4 * log(n)), 2)  # Size of bootstrap subsamples as in the paper. Minimal size of the subsample is equal to 2.
  }
  
  y <- sample(x, B * nn, replace = TRUE)
  Y <- matrix(y, nrow = B, ncol = nn)
  
  mu_star <- apply(matrix(Y, nrow = B, ncol = nn), 1, mean)
  MU <- mean(x)
  MM <- as.numeric(mu_star > xi * MU)
  
  return(pvalue = mean(MM))
}
