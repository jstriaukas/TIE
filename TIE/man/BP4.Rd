% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estim_functions.R
\name{BP4}
\alias{BP4}
\title{Baek and Pipiras (2010) Tail Index Estimator (Equation 2.19)}
\usage{
BP4(X, k)
}
\arguments{
\item{X}{A numeric vector containing the data sample.}

\item{k}{An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.}
}
\value{
A numeric value representing the estimated tail index.
}
\description{
Computes the tail index using the Baek and Pipiras (2010) estimator based on their equation 2.19, which combines the Hill estimator and second moments.
}
\details{
This estimator refines the tail index by incorporating both the Hill estimator and the second moment of logarithmic differences:
- The Hill estimator is calculated as:
  \[
  H = \text{Hill}(X, k).
  \]
- The second moment of logarithmic differences is computed as:
  \[
  M_2 = \text{MM}(X, k, 2),
  \]
  where \(\text{MM}(X, k, \alpha)\) computes the mean of the \(\alpha\)-th power of the logarithmic differences.
- The tail index is then calculated as:
  \[
  \text{Tail Index} = -H^2 + \left(\frac{1}{H} + 1\right) \cdot \frac{2}{M_2}.
  \]

This method combines the Hill estimator with higher-order moments to improve tail index estimation accuracy.
}
\examples{
set.seed(123)
n = 1000
x = (abs(stats::rcauchy(n)))^(2)
BP4(x)

}
