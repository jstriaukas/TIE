% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estim_functions.R
\name{GM4}
\alias{GM4}
\title{Gomes and Martin (2002) Estimator - Type 2}
\usage{
GM4(X, k, tau)
}
\arguments{
\item{X}{A numeric vector containing the data sample.}

\item{k}{An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.}

\item{tau}{An optional numeric parameter for auxiliary adjustment. Defaults to \code{0}.}
}
\value{
A numeric value representing the estimated tail index.
}
\description{
Computes the second estimator from Gomes and Martin (2002) for the tail index of a heavy-tailed distribution.
}
\details{
This estimator refines the tail index calculation by incorporating the auxiliary function \(\rho_1\) and two Hill estimators:
- The auxiliary function \(\rho_1\) is computed using the \code{rho} function with \(\tau = 0\).
- The inverse of the Hill estimator (\( g_2 \)) is computed using \code{k/2} order statistics.
- The inverse of the Hill estimator (\( g_3 \)) is computed using \( k \) order statistics.
- The tail index is calculated as:
  \[
  \text{Tail Index} = \frac{1 - 2^{-\rho_1}}{g_3 - 2^{-\rho_1} \cdot g_2}.
  \]

This method uses the auxiliary function \(\rho_1\) to adjust the Hill estimators, reducing bias and improving accuracy.
}
\examples{
set.seed(123)
n = 1000
x = (abs(stats::rcauchy(n)))^(2)
GM4(x)

}
