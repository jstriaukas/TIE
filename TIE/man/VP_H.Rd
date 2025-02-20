% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estim_functions.R
\name{VP_H}
\alias{VP_H}
\title{Vaiciulis and Paulauskas Hill Estimator Modification}
\usage{
VP_H(X, k, r)
}
\arguments{
\item{X}{A numeric vector containing the data sample.}

\item{k}{An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.}

\item{r}{An optional numeric parameter controlling the transformation. Defaults to \code{-0.1}. When \code{r = 0}, the logarithmic transformation is used.}
}
\value{
A numeric value representing the estimated tail index.
}
\description{
Computes the Vaiciulis and Paulauskas (2013) modification of the Hill estimator for the tail index of a heavy-tailed distribution.
}
\details{
The Vaiciulis and Paulauskas modification of the Hill estimator introduces a parameter \code{r} for a generalized transformation of the scaled top-order statistics:
- When \code{r = 0}, the transformation is logarithmic (\code{log(Y)}).
- When \code{r != 0}, the transformation is \code{(Y^r - 1) / r}, where \code{Y} is the ratio of the top-order statistics to the threshold value.
- The tail index is calculated using the mean of the transformed values.

This modification provides a more flexible framework for estimating the tail index by allowing different transformations based on the value of \code{r}.
}
\examples{
set.seed(123)
n = 1000
x = (abs(stats::rcauchy(n)))^(2)
VP_H(x)

}
