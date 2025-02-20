% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estim_functions.R
\name{GMV3}
\alias{GMV3}
\title{Gomes, Miranda, and Viseu (Statistica Neerlandica) Jackknife Estimator}
\usage{
GMV3(X, k)
}
\arguments{
\item{X}{A numeric vector containing the data sample.}

\item{k}{An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.}
}
\value{
A numeric value representing the estimated tail index.
}
\description{
Computes the Jackknife tail index estimator from Gomes, Miranda, and Viseu (Statistica Neerlandica, equation 26, page 11).
}
\details{
This Jackknife estimator refines the tail index by combining generalized estimators with an auxiliary parameter \(\rho\):
- The auxiliary parameter \(\rho\) is calculated using the \code{rho} function, with \( k \cdot 1.1 \) or \( 0.9 \cdot n \) order statistics, whichever is smaller.
- Two generalized estimators are computed:
  - \( g_1 \): The generalized estimator with \(\alpha = 1\).
  - \( g_2 \): The generalized estimator with \(\alpha = 1 - \rho\).
- The tail index is then calculated as:
  \[
  \text{Tail Index} = \frac{\rho^2}{(1 - \rho)^2 \cdot g_1 - (1 - 2 \cdot \rho) \cdot g_2}.
  \]

This method reduces bias by combining estimators with a weighting scheme influenced by \(\rho\).
}
\examples{
set.seed(123)
n = 1000
x = (abs(stats::rcauchy(n)))^(2)
GMV3(x)

}
