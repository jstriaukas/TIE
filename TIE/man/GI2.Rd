% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estim_functions.R
\name{GI2}
\alias{GI2}
\title{Gabaix & Ibragimov Harmonic Estimator}
\usage{
GI2(X, k)
}
\arguments{
\item{X}{A numeric vector containing the data sample.}

\item{k}{An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.}
}
\value{
A numeric value representing the estimated tail index.
}
\description{
Computes the Gabaix & Ibragimov Harmonic (Method 1) estimator for the tail index of a heavy-tailed distribution.
}
\details{
The Gabaix & Ibragimov Harmonic estimator uses harmonic sums of ranks as x-values and the logarithms of order statistics as y-values in a linear regression model. The reciprocal of the slope of the regression line is the estimated tail index.
}
\examples{
set.seed(123)
n = 1000
x = (abs(stats::rcauchy(n)))^(2)
GI2(x)

}
