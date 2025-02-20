% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estim_functions.R
\name{FH_MLE2}
\alias{FH_MLE2}
\title{Feuerverger and Hall (1999) OLS Tail Index Estimator}
\usage{
FH_MLE2(par, X, k)
}
\arguments{
\item{par}{A numeric vector of length 3 containing the parameters \code{mu}, \code{D}, and \code{beta} to be optimized.}

\item{X}{A numeric vector containing the data sample.}

\item{k}{An optional integer specifying the number of top-order statistics to use. Defaults to \(2 \cdot \log(n)\), where \( n \) is the sample size.}
}
\value{
A numeric value representing the sum of squared residuals for the given parameters \code{mu}, \code{D}, and \code{beta}.
}
\description{
Computes the tail index estimator using the Ordinary Least Squares (OLS) approach proposed by Feuerverger and Hall (1999).
}
\details{
This function implements the OLS-based approach to tail index estimation:
- The parameters \code{mu}, \code{D}, and \code{beta} are optimized by minimizing the sum of squared residuals.
- Intermediate quantities are computed as follows:
  - \( U \): The scaled differences between order statistics:
    \[
    U = i \cdot \left(X_{(n-i+1)} - X_{(n-i)}\right), \, i = 1, \ldots, n-1.
    \]
  - \( V \): The logarithm of \( U \):
    \[
    V = \log(U).
    \]
- The sum of squared residuals is calculated as:
  \[
  \text{Residuals} = \sum_{i=1}^{k} \left(V_i - \mu - D \cdot \left(\frac{i}{n}\right)^{-\beta}\right)^2.
  \]

This function is designed for use in an optimization routine to find the best-fitting values for \code{mu}, \code{D}, and \code{beta}.
}
\examples{
set.seed(123)
n = 1000
x = (abs(stats::rcauchy(n)))^(2)
params <- c(0, 1, 0.5)  # Initial parameter guesses
FH_MLE2(params, x)

}
