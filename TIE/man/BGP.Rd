% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estim_functions.R
\name{BGP}
\alias{BGP}
\title{Brilhante et al. (2013) Tail Index Estimator}
\usage{
BGP(X, k, p)
}
\arguments{
\item{X}{A numeric vector containing the data sample.}

\item{k}{An optional integer specifying the number of top-order statistics to use. If not provided, it is calculated using the method from Gomes et al. (2016).}

\item{p}{An optional numeric parameter controlling the power transformation. If not provided, it is estimated using auxiliary parameters.}
}
\value{
A numeric value representing the estimated tail index.
}
\description{
Computes the tail index estimator from Brilhante et al. (2013), incorporating auxiliary bias adjustments and optimal parameter selection.
}
\details{
This estimator adjusts the tail index using auxiliary parameters and optimal threshold selection:
- The auxiliary parameters \(\rho\) and \(\beta\) are computed using the \code{rho} and \code{beta} functions.
- If \( k \) is not provided, it is calculated as:
  \[
  k = \left( \frac{(1 - \rho) \cdot n^{-\rho}}{\beta \cdot \sqrt{-2 \rho}} \right)^{\frac{2}{1 - 2 \rho}},
  \]
  with constraints \( 0.1 \cdot \sqrt{n} \leq k \leq 0.9 \cdot n \).
- If \( p \) is not provided, it is estimated as:
  \[
  p = \frac{\phi}{H},
  \]
  where \( \phi = \max\left(0, \min\left(\frac{\sqrt{2}}{2}, 1 - \frac{\rho}{2} - 0.5 \sqrt{\rho^2 - 4\rho + 2}\right)\right) \) and \( H \) is the inverse Hill estimator.
- The tail index is computed as:
  \[
  \alpha = \frac{1}{\gamma}, \quad \text{where} \quad \gamma = \frac{1 - \frac{1}{\text{mean}(U)}}{p},
  \]
  and \( U = \left( \frac{X_{(n-k+i)}}{X_{(n-k)}} \right)^p \) for \( i = 1, \ldots, k \).
}
\examples{
set.seed(123)
n = 1000
x = (abs(stats::rcauchy(n)))^(2)
BGP(x)

}
