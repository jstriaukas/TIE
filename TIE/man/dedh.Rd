% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estim_functions.R
\name{dedh}
\alias{dedh}
\title{\insertCite{dekkers1989;textual}{TIE} moment estimator}
\usage{
dedh(X, k)
}
\arguments{
\item{X}{A numeric vector containing the data sample.}

\item{k}{An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.}
}
\value{
A numeric value representing the estimated tail index.
}
\description{
Computes the \insertCite{dekkers1989;textual}{TIE}  moment estimator for the tail index of a heavy-tailed distribution.
}
\details{
The \insertCite{dekkers1989;textual}{TIE} estimator is based on the moments of the logarithmic differences of extreme order statistics. It is designed to reduce bias in tail index estimation by incorporating higher-order moments.
}
\examples{
set.seed(123)
n = 1000
x = (abs(stats::rcauchy(n)))^(2)
dedh(x)
}
\references{
\insertRef{dekkers1989}{TIE}
}
