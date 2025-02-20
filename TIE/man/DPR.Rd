% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estim_functions.R
\name{DPR}
\alias{DPR}
\title{DPR Estimator with Permutations}
\usage{
DPR(X, M)
}
\arguments{
\item{X}{A numeric vector containing the data sample.}

\item{M}{An optional integer specifying the number of permutations to perform. Defaults to the square root of the sample size.}
}
\value{
A numeric value representing the estimated tail index.
}
\description{
Computes the DPR estimator for the tail index of a heavy-tailed distribution using multiple permutations of the data.
}
\details{
The DPR estimator applies the \code{DPR1} function to multiple permutations of the input data. For each permutation, the tail index is computed, and the final estimate is obtained as the mean of the estimates across all permutations. This approach enhances robustness to variations in the data.
}
\examples{
set.seed(123)
n = 1000
x = (abs(stats::rcauchy(n)))^(2)
DPR(x)

}
