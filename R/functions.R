#' Hill Estimator
#'
#' Computes the Hill estimator, a commonly used method for estimating the tail index of heavy-tailed distributions.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The Hill estimator is a semi-parametric estimator based on the order statistics of the data. It is used to estimate the tail index, which characterizes the heaviness of the tail of a distribution.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  #tail index = 0.5.
#' Hill(x)
#'
#' @export
Hill <- function(X, k) {
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  (sum(log(X[(n-k+1):n]))/k - log(X[n-k]))^(-1)
}

#' De Haan & Resnick Estimator
#'
#' Computes the De Haan & Resnick (1981) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The De Haan & Resnick estimator is based on the logarithmic differences of the top-order statistics. It provides an alternative to the Hill estimator for tail index estimation.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' dHR(x)
#'
#' @export
dHR <- function(X, k) {
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  log(k) / (log(X[n]) - log(X[n-k+1]))
}

#' Bacro Brito Estimator
#'
#' Computes the Bacro Brito (1995) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param alpha An optional numeric value between 0 and 1 specifying the adjustment factor. Defaults to 0.9.
#' @return A numeric value representing the estimated tail index.
#' @details The Bacro Brito estimator uses a refined \code{alpha} value to calculate the tail index, providing an alternative to other estimators like the Hill estimator.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' BB(x)
#'
#' @export
BB <- function(X, k, alpha) {
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(alpha)) alpha <- 0.9
  f <- floor(alpha * k)
  alpha <- f / k          # more exact alpha
  -log(alpha) / (log(X[n-f+1]) - log(X[n-k+1]))
}

#' Pickands Estimator
#'
#' Computes the Pickands (1975) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size. The value of \code{k} is adjusted to be divisible by 4.
#' @return A numeric value representing the estimated tail index.
#' @details The Pickands estimator is based on the ratios of spacings of extreme order statistics. It is a commonly used method for estimating the tail index of heavy-tailed distributions.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' Pickands(x)
#'
#' @export
Pickands <- function(X, k) {
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  k <- (k %/% 4 + 1) * 4  # Ensure k is divisible by 4
  log(2) / (log((X[n-k/4] - X[n-k/2]) / (X[n-k/2] - X[n-k])))
}

#' Dekkers, Einmahl, De Haan "Moment" Estimator
#'
#' Computes the Dekkers, Einmahl, De Haan (1989) "moment" estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The DEdH estimator is based on the moments of the logarithmic differences of extreme order statistics. It is designed to reduce bias in tail index estimation by incorporating higher-order moments.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' DEdH(x)
#'
#' @export
DEdH <- function(X, k) {
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  M1 <- sum(log(X[(n-k+1):n]) - log(X[n-k])) / k
  M2 <- sum((log(X[(n-k+1):n]) - log(X[n-k]))^2) / k
  (M1 + 1 - 0.5 * (1 - M1^2 / M2)^(-1))^(-1)
}

#' Aban & Meerschaert Shifted Hill's Estimator
#'
#' Computes the Aban & Meerschaert (2001) shifted Hill’s estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param eps An optional numeric value specifying the tolerance level for the iterative algorithm. Defaults to \code{1e-7}.
#' @return A numeric value representing the estimated tail index.
#' @details The Aban & Meerschaert shifted Hill’s estimator adjusts the Hill estimator by introducing a shift parameter to better handle small sample sizes and reduce bias. The iterative algorithm continues until the difference between the maximum and minimum shift parameter estimates is below \code{eps}.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' AB_Hill(x)
#'
#' @export
AB_Hill <- function(X, k, eps = 1e-7) {
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  s_min <- -10 * X[n-k]
  s_max <- X[n-k]
  repeat {
    s1 <- (s_min + s_max) / 2
    
    alpha <- (sum(log(X[(n-k+1):n] - s1)) / k - log(X[n-k] - s1))^(-1)
    Z <- sum((X[(n-k+1):n] - s1)^(-1))
    RHS <- (alpha + 1) * Z / k
    LHS <- alpha * (X[n-k] - s1)^(-1)
    
    if (RHS > LHS) s_min <- s1 else s_max <- s1
    if ((s_max - s_min) < eps) break
  }
  alpha
}


#' Peng Estimator (Hill Modification)
#'
#' Computes the Peng estimator, a modification of the Hill estimator, for estimating the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The Peng estimator refines the Hill estimator by incorporating additional moments and parameters to reduce bias. It adjusts the tail index estimation using intermediate order statistics.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' Peng_H(x)
#'
#' @export
Peng_H <- function(X, k) {
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  k1 <- round(n / (2 * log(n)))
  k2 <- round(n / log(n))
  g1 <- Hill(X, k1)^(-1)
  g2 <- Hill(X, k2)^(-1)
  g <- Hill(X, k)^(-1)
  M1 <- MM(X, k1)
  M2 <- MM(X, k2)
  M <- MM(X, k)
  M4 <- (M1 - 2 * g1^2) / (M2 - 2 * g2^2)
  if (M4 > 0) rho <- log(2)^(-1) * log(M4) else rho <- 1
  g7 <- g - (M - 2 * g^2) * (1 - rho) / (2 * g * rho)
  g7^(-1)
}

#' Peng Estimator (Pickands Modification)
#'
#' Computes the Peng estimator, a modification of the Pickands estimator, for estimating the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The Peng estimator modifies the Pickands estimator by adjusting the tail index computation using intermediate order statistics and introducing a parameter \code{rho} for bias correction.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' Peng_P(x)
#'
#' @export
Peng_P <- function(X, k) {
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  k1 <- round(n / log(n))
  k2 <- round(n / (2 * log(n)))
  k3 <- round(n / (4 * log(n)))
  g1 <- Pickands(X, k1)^(-1)
  g2 <- Pickands(X, k2)^(-1)
  g3 <- Pickands(X, k3)^(-1)
  g <- Pickands(X, k)^(-1)
  
  M <- (g2 - g3) / (g1 - g2)
  if (M > 0) {
    rho <- log(M) / log(2)
    g7 <- g - (g - g3) / (1 - 4^rho)
  } else {
    g7 <- g
  }
  g7^(-1)
}

#' Fialova, Jureckova, Picek (FJP) Estimator
#'
#' Computes the Fialova, Jureckova, Picek (2004) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param N An optional integer specifying the number of partitions. Defaults to 5.
#' @param alpha_0 An optional initial estimate of the tail index. Defaults to 1.2 times the Hill estimator.
#' @param delta An optional numeric parameter controlling the adjustment of \code{a_N}. Defaults to 0.9.
#' @return A numeric value representing the estimated tail index.
#' @details The FJP estimator divides the sample into partitions, computes means for each partition, and calculates the tail index using an adjusted value of \code{a_N} and the empirical proportion of values below \code{a_N}. If the proportion \code{F} is not in the range (0, 1), the initial estimate \code{alpha_0} is returned.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' FJP(x)
#'
#' @export
FJP <- function(X, N, alpha_0, delta) {
  n <- length(X)
  if (missing(N)) N <- 5
  if (missing(alpha_0)) alpha_0 <- Hill(X) * 1.2
  if (missing(delta)) delta <- 0.9
  
  r <- n %/% N
  r0 <- n %% N
  M <- matrix(X[1:(r * n)], nrow = r, ncol = N)
  XX <- apply(M, 1, mean)
  if (r0 > 0) XX <- c(XX, mean(X[n - r0:n]))
  
  a_N <- N^((1 - delta) / alpha_0)
  
  F <- mean(XX <= a_N)
  
  if ((F > 0) & (F < 1)) {
    alpha <- -log(1 - F) / log(a_N)
  } else {
    alpha <- alpha_0
  }
  alpha
}

#' Van Zyl Estimator
#'
#' Computes the Van Zyl (2015) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The Van Zyl estimator incorporates the Hill estimator along with location and scale parameters to compute the tail index. It uses an adjusted logarithmic mean of the top-order statistics for better accuracy.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' vanZyl(x)
#'
#' @export
vanZyl <- function(X, k) {
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  h_gamma <- Hill(X, k)
  h_mu <- X[n - k]
  h_sigma <- h_mu * h_gamma
  TT <- mean(log(X[(n - k + 1):n] * h_gamma / h_sigma + 1 - h_mu * h_gamma / h_sigma))
  (TT)^(-1)
}

#' Nuyts Estimator
#'
#' Computes the Nuyts (2015) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param eps An optional numeric value specifying the tolerance level for the iterative algorithm. Defaults to \code{1e-9}.
#' @return A numeric value representing the estimated tail index.
#' @details The Nuyts estimator uses a binary search algorithm to determine the tail index. The iterative process continues until the difference between the maximum and minimum estimates is below \code{eps}. If the difference exceeds a threshold, a second iteration is performed to refine the estimate.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' Nuyts(x)
#'
#' @export
Nuyts <- function(X, k, eps = 1e-9) {
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  max_a <- 20
  min_a <- 0
  repeat {
    a <- (max_a + min_a) / 2
    L <- mean(log(X[(n - k + 1):n]))
    R <- 1 / a + (log(X[n - k + 1]) * X[n - k]^(-a) - log(X[n]) * X[n]^(-a)) / 
      (X[n - k + 1]^(-a) - X[n]^(-a))
    if (L > R) max_a <- a else min_a <- a
    if (abs(max_a - min_a) < eps) break
  }
  if ((20 - max_a) > eps) {
    a1 <- a
  } else {
    max_a <- 20
    min_a <- 0
    repeat {
      a <- (max_a + min_a) / 2
      L <- mean(log(X[(n - k + 1):n]))
      R <- 1 / a + (log(X[n - k + 1]) * X[n - k]^(-a) - log(X[n]) * X[n]^(-a)) / 
        (X[n - k + 1]^(-a) - X[n]^(-a))
      if (L > R) min_a <- a else max_a <- a
      if (abs(max_a - min_a) < eps) break
    }
    a1 <- a
  }
  a1
}


VBCH1 <- function(theta, X, k) {
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  gamma <- theta[1]
  delta <- theta[2]
  rho <- theta[3]
  w <- theta[4]
  R1 <- (1 + delta)^2 / (gamma + 2) / gamma
  R2 <- 2 * delta * (1 - delta) * (1 - gamma * rho) / (gamma + 2 - gamma * rho) / gamma
  R3 <- delta^2 * (1 - gamma * rho)^2 / (gamma + 2 - 2 * gamma * rho) / gamma
  R <- R1 + R2 + R3
  Y <- X[(n - k + 1):n] / X[n - k]
  f1 <- (1 - delta) * Y^(-(1 + 1 / gamma)) / gamma
  f2 <- delta * (1 / gamma - rho) * Y^(-(1 + 1 / gamma - rho))
  ff <- sum(f1 + f2)
  w^2 * R - 2 * w * ff / k
}

#' Vandewalle, Beirlant, Christmann, and Hubert (VBCH) Estimator
#'
#' Computes the VBCH (2007) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The VBCH estimator minimizes the objective function defined by the \code{VBCH1} function using the \code{optim} function with the \code{L-BFGS-B} method. It uses starting values and constraints for the parameters (\code{gamma}, \code{delta}, \code{rho}, and \code{w}). The tail index is computed as the reciprocal of the optimized value of \code{gamma}.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' VBCH(x)
#'
#' @export
VBCH <- function(X, k) { 
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  pp <- optim(
    c(1, 0.1, 0.1, 0.5), VBCH1, X = X, k = k, method = "L-BFGS-B",
    lower = c(0.001, 0.001, 0.001, 0.001), upper = c(Inf, 1, Inf, 1)
  )
  1 / pp$par[1]
}

#' Tripathi, Kumar, Petropoulos (TKP1) Estimator
#'
#' Computes the TKP1 (2014) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The TKP1 estimator is based on the ratio of logarithmic differences of top-order statistics. It is designed to provide an efficient estimate of the tail index using equation 3.2 from the original paper.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' TKP1(x)
#'
#' @export
TKP1 <- function(X, k) {
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  (k - 3) / (sum(log(X[(n - k + 1):n] / X[n - k])))
}
#' Tripathi, Kumar, Petropoulos (TKP2) Estimator
#'
#' Computes the TKP2 (2014) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The TKP2 estimator is an extension of the TKP1 estimator that includes an additional term \code{max(0, log(X[n-k]))} in the denominator to improve robustness. It is based on equation 3.2 from the original paper.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' TKP2(x)
#'
#' @export
TKP2 <- function(X, k) {
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  (k - 3) / (sum(log(X[(n - k + 1):n] / X[n - k])) + max(0, log(X[n - k])))
}

#' Muller Estimator
#'
#' Computes the Muller estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param m An optional integer specifying the number of observations to use. Defaults to the square root of the sample size.
#' @param a1 An optional numeric value representing the first threshold parameter. Defaults to one-fourth of the range of the top \code{m} observations.
#' @param a2 An optional numeric value representing the second threshold parameter. Defaults to half of the range of the top \code{m} observations.
#' @return A numeric value representing the estimated tail index.
#' @details The Muller estimator calculates the tail index using two threshold parameters (\code{a1} and \code{a2}) and indicator functions applied to the data. The thresholds are adaptively determined based on the range of the top \code{m} observations if not provided.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' Muller(x)
#'
#' @export
Muller <- function(X, m, a1, a2) { 
  X <- sort(X)
  n <- length(X)
  if (missing(m)) m <- round(sqrt(n))
  if (missing(a1)) a1 <- (X[n - m + 1] - X[1]) / 4
  if (missing(a2)) a2 <- (X[n - m + 1] - X[1]) / 2
  
  X1 <- as.numeric(X < X[n - m + 1])
  X2 <- as.numeric(X >= (X[n - m + 1] - a1))
  X3 <- as.numeric(X >= (X[n - m + 1] - a2))
  M1 <- sum(X1 * X2)
  M2 <- sum(X2 * X3)
  
  (log(a2) - log(a1)) / (log(M1) - log(M2))
}


#' Muller and Rufibach Smoothed Pickands Estimator
#'
#' Computes the Muller and Rufibach (2009) smoothed Pickands estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The smoothed Pickands estimator employs logspline density estimation to smooth the quantiles of the data before calculating the tail index. It uses three quantiles (\code{X1}, \code{X2}, and \code{X4}) derived from the logspline density estimate to refine the Pickands ratio calculation.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' MRsP(x)
#'
#' @export
MRsP <- function(X, k) {              
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  r <- k / 4
  param <- logspline::logspline(X)
  X1 <- logspline::qlogspline(((n - r + 1) / n), param) 
  X2 <- logspline::qlogspline(((n - 2 * r + 1) / n), param) 
  X4 <- logspline::qlogspline(((n - 4 * r + 1) / n), param)
  log(2) / (log((X1 - X2) / (X2 - X4)))
}

#' Meerschaert and Scheffer Estimator
#'
#' Computes the Meerschaert and Scheffer (1998) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @return A numeric value representing the estimated tail index.
#' @details The Meerschaert and Scheffer estimator calculates the tail index using the variance of the sample around its mean. It is based on the relationship between the logarithm of the sample size and the logarithm of the sample variance.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' MS(x)
#'
#' @export
MS <- function(X) {
  n <- length(X)
  m <- mean(X)
  S <- sum((X - m)^2)
  (2 * log(n)) / max(0, log(S))
}

pl <- function(X) {          
  n <- length(X)
  I <- c(1:n)
  
  S <- cumsum(X^2) / I
  Y <- log(S)
  
  Yb <- mean(Y)
  Ib <- mean(log(I))
  
  D1 <- sum((Y - Yb) * (log(I) - Ib))
  D2 <- sum((log(I) - Ib)^2)
  2 / (D1 / D2 + 1)
}

#' Politis Estimator (2002)
#'
#' Computes the Politis (2002) tail index estimator for a heavy-tailed distribution using a bootstrap approach.
#'
#' @param X A numeric vector containing the data sample.
#' @return A numeric value representing the estimated tail index.
#' @details The Politis estimator applies a bootstrap procedure to compute the median of tail index estimates derived from resampled datasets. The tail index for each resampled dataset is calculated using the \code{pl} function, which estimates the power-law tail index.
#' 
#' The function uses \code{m = 1000} bootstrap iterations by default.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' Politis(x)
#'
#' @export
Politis <- function(X) { # Politis 2002
  n <- length(X)
  m <- 1000  # Number of bootstraps
  alpha <- vector(length = m)
  for (i in 1:m) {
    Y <- sample(X, n, replace = FALSE)
    alpha[i] <- pl(Y)
  }
  median(alpha)
}

#' McElroy and Politis BAS Estimator
#'
#' Computes the McElroy and Politis (2006) BAS estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param r An optional positive integer specifying the moment order. Defaults to \code{1}. \code{r} should be relatively large, such that the \code{2r}-th moment does not exist.
#' @return A numeric value representing the estimated tail index.
#' @details The McElroy and Politis BAS estimator calculates the tail index based on the \code{2r}-th moment of the sample. The method involves a logarithmic transformation of the sum of \code{2r}-th powers of the data and the sample size.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' MEP_BAS(x)
#'
#' @export
MEP_BAS <- function(X, r) {          
  n <- length(X)
  if (missing(r)) r <- 1
  S <- sum((X)^(2 * r))
  (2 * r * log(n)) / log(S)
}
#' McElroy and Politis CEN Estimator
#'
#' Computes the McElroy and Politis (2006) CEN estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @return A numeric value representing the estimated tail index.
#' @details The McElroy and Politis CEN estimator computes the tail index based on the second moments of the entire sample and the top \code{k} observations, where \code{k} is the square root of the sample size. This method relies on the logarithmic difference of these second moments.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' MEP_CEN(x)
#'
#' @export
MEP_CEN <- function(X) {          
  n <- length(X)
  k <- round(sqrt(n))
  
  S1 <- sum(X^(2))
  S2 <- sum(X[1:k]^(2))
  (2 * log(k)) / (log(S1) - log(S2))
}

#' McElroy and Politis CEN Estimator
#'
#' Computes the McElroy and Politis (2006) CEN estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @return A numeric value representing the estimated tail index.
#' @details The McElroy and Politis CEN estimator computes the tail index based on the second moments of the entire sample and the top \code{k} observations, where \code{k} is the square root of the sample size. This method relies on the logarithmic difference of these second moments.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' MEP_CEN(x)
#'
#' @export
MEP_CEN <- function(X) {          
  n <- length(X)
  k <- round(sqrt(n))
  
  S1 <- sum(X^(2))
  S2 <- sum(X[1:k]^(2))
  (2 * log(k)) / (log(S1) - log(S2))
}

#' McElroy and Politis RCEN Estimator
#'
#' Computes the McElroy and Politis (2006) RCEN (revised CEN) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @return A numeric value representing the estimated tail index.
#' @details The RCEN estimator divides the data into blocks of size \code{b^2}, where \code{b} is the cube root of the sample size (\code{n^(1/3)}). The last block is formed from the remaining observations. For each block, the CEN estimator is applied, and the final tail index is computed as the reciprocal of the mean of the reciprocal tail indices for all blocks.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' MEP_SCEN(x)
#'
#' @export
MEP_SCEN <- function(X) {          
  n <- length(X)
  b <- round(n^(1/3))
  M <- floor(n / (b^2)) - 1
  XX <- matrix(X[1:(M * b^2)], nrow = b^2, ncol = M)
  x <- X[(M * b^2 + 1):n]  # Last group formed out of b^2 last observations + remaining observations
  SS <- vector(length = M)
  for (i in 1:M) {
    SS[i] <- MEP_CEN(XX[, i])
  }
  ss <- MEP_CEN(x)
  mean((c(SS, ss))^(-1))^(-1)  # Averaging reciprocal tail indices
}

#' McElroy and Politis RCEN Estimator
#'
#' Computes the McElroy and Politis (2006) RCEN (revised CEN) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @return A numeric value representing the estimated tail index.
#' @details The RCEN estimator divides the sample into \code{b} blocks of size \code{sqrt(n)}, where \code{b} is the floor of the square root of the sample size. For each block, the estimator is calculated using the logarithmic difference of the second moments of the entire sample and the block. The final tail index is computed as the reciprocal of the mean of the block estimates.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' MEP_RCEN(x)
#'
#' @export
MEP_RCEN <- function(X) {          
  n <- length(X)
  b <- floor(sqrt(n))
  B <- vector(length = b)
  for (i in 1:b) {
    S1 <- sum(X^(2))
    S2 <- sum(X[((i - 1) * b + 1):(i * b)]^(2))
    B[i] <- (log(S1) - log(S2)) / (2 * log(b))
  }
  mean(B)^(-1)
}

#' McElroy and Politis SRCEN Estimator
#'
#' Computes the McElroy and Politis (2006) SRCEN (smoothed RCEN) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @return A numeric value representing the estimated tail index.
#' @details The SRCEN estimator divides the data into blocks of size \code{b^2}, where \code{b} is the cube root of the sample size (\code{n^(1/3)}). The last block is formed from the remaining observations. For each block, the RCEN estimator is applied, and the final tail index is computed as the reciprocal of the mean of the reciprocal tail indices for all blocks.
#' 
#' This approach smooths the RCEN estimates by averaging them, which provides robustness to the tail index estimation. The median of reciprocal tail indices may perform better in practice.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' MEP_SRCEN(x)
#'
#' @export
MEP_SRCEN <- function(X) {          
  n <- length(X)
  b <- round(n^(1/3))
  M <- floor(n / (b^2)) - 1
  XX <- matrix(X[1:(M * b^2)], nrow = b^2, ncol = M)
  x <- X[(M * b^2 + 1):n]  # Last group formed out of b^2 last observations + remaining observations
  SS <- vector(length = M)
  for (i in 1:M) {
    SS[i] <- MEP_RCEN(XX[, i])
  }
  ss <- MEP_RCEN(x)
  mean((c(SS, ss))^(-1))^(-1)  # Averaging reciprocal tail indices
}

#' Falk Estimator
#'
#' Computes the Falk (1994) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The Falk estimator is a weighted combination of two Pickands estimators: one using \code{k/2} order statistics and the other using \code{k} order statistics. The weights are set to \code{5/13} and \code{8/13}, respectively, to balance the contribution of the two estimators.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' Falk(x)
#'
#' @export
Falk <- function(X, k) {           
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  P1 <- Pickands(X, (k / 2))
  P2 <- Pickands(X, k)
  P1 * 5 / 13 + P2 * 8 / 13
}

#' Zipf Estimator (Quantile Plots Estimator)
#'
#' Computes the Zipf estimator, based on the Kratz and Resnick (1996) quantile plots approach, for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param w An optional vector of weights for the linear regression. Defaults to equal weights for all points.
#' @return A numeric value representing the estimated tail index.
#' @details The Zipf estimator uses a log-log plot of order statistics to estimate the tail index. A weighted linear regression is performed on the log-transformed data, where \code{w} specifies the weights for each point. The slope of the regression line is used to calculate the tail index.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' Zipf(x)
#'
#' @export
Zipf <- function(X, k, w) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(w)) w <- rep(1, k)
  
  x <- -log(c(n:1) / (n + 1))[(n - k + 1):n]
  y <- log(X[(n - k + 1):n])
  f <- lm(y ~ x, weights = w)
  1 / f$coefficients[2]
}

#' Aban & Meerschaert Hill Estimator
#'
#' Computes the Aban & Meerschaert (2001) adjusted Hill estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The Aban & Meerschaert Hill estimator adjusts the standard Hill estimator by multiplying it by \code{k / (k - 1)} to correct for bias. This estimator is especially useful for small sample sizes where the bias in the Hill estimator is significant.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' AM1(x)
#'
#' @export
AM1 <- function(X, k) {          
  if (missing(k)) k <- round(sqrt(length(X)))
  Hill(X, k) * k / (k - 1)
}

#' Aban & Meerschaert Simplified Estimator
#'
#' Computes the Aban & Meerschaert simplified estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The simplified Aban & Meerschaert estimator calculates the tail index using a weighted average of logarithmic values of the top-order statistics. The weights are derived from an adjusted cumulative sum of ranks. This estimator aims to provide a more straightforward alternative to the Hill estimator, incorporating a bias-correction mechanism.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' AM2(x)
#'
#' @export
AM2 <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  j <- c(n:1)
  
  a1 <- rev(cumsum(j^(-1)))
  a <- a1[(n - k + 1):n]
  x <- X[(n - k + 1):n]
  b_a <- mean(a)
  s <- b_a * (a - b_a) / (sum((a - b_a)^2))
  (sum(s * log(x)))^(-1)
}

#' Aban & Meerschaert Smoothing Regression Estimator
#'
#' Computes the Aban & Meerschaert smoothing regression estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The Aban & Meerschaert smoothing regression estimator uses a linear regression approach. The x-values for the regression are the negative logarithms of adjusted empirical ranks, and the y-values are the logarithms of the corresponding order statistics. The slope of the fitted regression line is used to calculate the tail index.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' AM3(x)
#'
#' @export
AM3 <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  
  x <- -log((c(n:1) - 0.5) / n)[(n - k + 1):n]
  y <- log(X[(n - k + 1):n])
  f <- lm(y ~ x)
  1 / f$coefficients[2]
}

#' Gabaix & Ibragimov Reversed Estimator
#'
#' Computes the Gabaix & Ibragimov reversed estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The Gabaix & Ibragimov reversed estimator applies a linear regression approach where the x-values are the logarithms of the order statistics and the y-values are the negative logarithms of adjusted empirical ranks. The slope of the fitted regression line is the tail index.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GI1(x)
#'
#' @export
GI1 <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  
  x <- -log((c(n:1) - 0.5) / n)[(n - k + 1):n]
  y <- log(X[(n - k + 1):n])
  f <- lm(x ~ y)
  f$coefficients[2]
}

#' Gabaix & Ibragimov Harmonic Estimator
#'
#' Computes the Gabaix & Ibragimov Harmonic (Method 1) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The Gabaix & Ibragimov Harmonic estimator uses harmonic sums of ranks as x-values and the logarithms of order statistics as y-values in a linear regression model. The reciprocal of the slope of the regression line is the estimated tail index.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GI2(x)
#'
#' @export
GI2 <- function(X, k) {          
  X <- sort(X, decreasing = TRUE)
  n <- length(X)
  
  if (missing(k)) k <- round(sqrt(n))
  I <- 1 / c(1:n)
  x1 <- -cumsum(I)
  
  H <- x1[1:(k - 1)]
  y <- log(X[2:k])
  
  f <- lm(y ~ H)
  1 / f$coefficients[2]
}

#' Gabaix & Ibragimov Harmonic Estimator (Method 2)
#'
#' Computes the Gabaix & Ibragimov Harmonic (Method 2) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The Gabaix & Ibragimov Harmonic estimator (Method 2) reverses the linear regression setup by using the harmonic sums of ranks as the dependent variable (\code{H}) and the logarithms of the order statistics as the independent variable (\code{y}). The slope of the regression line is the estimated tail index.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GI3(x)
#'
#' @export
GI3 <- function(X, k) {          
  X <- sort(X, decreasing = TRUE)
  n <- length(X)
  
  if (missing(k)) k <- round(sqrt(n))
  I <- 1 / c(1:n)
  x1 <- -cumsum(I)
  
  H <- x1[1:(k - 1)]
  y <- log(X[2:k])
  
  f <- lm(H ~ y)
  f$coefficients[2]
}

#' Beirlant, Vynckier, and Teugels (BVT) Estimator
#'
#' Computes the BVT estimator for the tail index of a heavy-tailed distribution using excess functions.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The BVT estimator calculates the tail index using the excess functions \code{UH}, which depend on the order statistics and the Hill estimator. A linear regression is performed on the logarithms of the excess functions and the negative logarithms of the empirical probabilities. The reciprocal of the slope of the regression line is the estimated tail index.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' BVT(x)
#'
#' @export
BVT <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  UH <- rep(0, k)
  
  for (j in 1:k) {
    UH[j] <- X[n - j] * (Hill(X, j)^(-1))
  }
  
  x <- -log(c(1:(k)) / n)
  y <- log(UH[1:(k)])
  f <- lm(y ~ x)
  1 / f$coefficients[2]
}

#' Beirlant, Dierckx, and Guillou (BDG) Estimator
#'
#' Computes the BDG (2005) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The BDG estimator calculates the tail index using the excess functions \code{UH}, derived from the order statistics and the Hill estimator. It involves a weighted mean of the logarithms of the excess functions, with weights based on the deviation of the empirical log-rank values from their mean.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' BDG(x)
#'
#' @export
BDG <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  UH <- rep(0, k)
  for (j in 1:k) {
    UH[j] <- X[n - j] * (Hill(X, j)^(-1))
  }
  
  I <- log((k + 1) / c(1:k))
  mI <- mean(I)
  (mean((I - mI) * log(UH)))^(-1)
}

#' Beirlant, Dierckx, and Guillou (BDG) Estimator - Hill Modification
#'
#' Computes the BDG (2005) estimator with the Hill modification for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The BDG Hill modification uses the excess functions \code{UH}, derived from the order statistics and the Hill estimator, as input to the Hill estimator. This method refines the original BDG estimator by further reducing bias.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' BDG_H(x)
#'
#' @export
BDG_H <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  UH <- rep(0, k)
  for (j in 1:k) {
    UH[j] <- X[n - j] * (Hill(X, j)^(-1))
  }
  Hill(UH, (k - 1))
}

#' Schultze and Steinebach Estimator
#'
#' Computes the Schultze and Steinebach (1996) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The Schultze and Steinebach estimator uses a linear regression approach without an intercept. The logarithms of scaled ranks are used as the independent variable, and the logarithms of the top-order statistics are used as the dependent variable. The ratio of the sums of squared logarithms and the cross-product of the logarithms provides the tail index estimate.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' SS(x)
#'
#' @export
SS <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  I <- n / c(1:k)
  Y1 <- sum(log(I)^2)
  Y2 <- sum(log(I) * log(X[n:(n - k + 1)]))
  Y1 / Y2
}

#' Schultze and Steinebach Estimator (Reversed Variables)
#'
#' Computes the Schultze and Steinebach (1996) estimator for the tail index of a heavy-tailed distribution, using reversed explanatory and dependent variables.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This version of the Schultze and Steinebach estimator reverses the roles of the explanatory and dependent variables compared to the original formulation. A linear regression is performed where the logarithms of the top-order statistics (\code{y}) are used as the explanatory variable, and the negative logarithms of scaled ranks (\code{x}) are the dependent variable. The slope of the regression line provides the tail index.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' SS2(x)
#'
#' @export
SS2 <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  
  x <- -log(c(n:1) / (n + 1))[(n - k + 1):n]
  y <- log(X[(n - k + 1):n])
  f <- lm(x ~ y)
  f$coefficients[2]
}

#' Beirlant, Vynckier, and Teugels (BVT) Estimator with Anchor Point
#'
#' Computes the BVT (1996) estimator for the tail index of a heavy-tailed distribution, using regression through the anchor point.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param w An optional vector of weights for the regression. Defaults to equal weights for all points.
#' @return A numeric value representing the estimated tail index.
#' @details The BVT estimator with anchor point performs a weighted linear regression through the origin. The x-values are adjusted log-ranks, and the y-values are adjusted logarithms of the top-order statistics. The slope of the regression line is used to compute the tail index as its reciprocal.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' BVT2(x)
#'
#' @export
BVT2 <- function(X, k, w) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(w)) w <- rep(1, k)
  
  x <- -log(c(n:1) / (n + 1))[(n - k + 1):n] + log((k + 1) / (n + 1))
  y <- log(X[(n - k + 1):n]) - log(X[n - k])
  f <- lm(y ~ 0 + x, weights = w)
  
  1 / f$coefficients[1]
}

#' Brito and Freitas Geometric Mean Estimator
#'
#' Computes the Brito and Freitas estimator for the tail index of a heavy-tailed distribution. This estimator is the geometric mean of two estimators: the Zipf estimator and the reversed Schultze and Steinebach estimator.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The Brito and Freitas estimator combines the Zipf and the reversed Schultze and Steinebach estimators by taking their geometric mean. This approach leverages the strengths of both methods for more robust tail index estimation.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' BF1(x)
#'
#' @export
BF1 <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  
  g1 <- Zipf(X, k)
  g2 <- SS2(X, k)
  
  sqrt(g1 * g2)
}

#' Danielsson, Jansen, and De Vries Estimator
#'
#' Computes the Danielsson, Jansen, and De Vries (1996) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The DJV estimator calculates the tail index using the ratio of the first and second moments of the logarithmic differences of the top-order statistics. This method is designed to provide a robust estimate of the tail index.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' DJV(x)
#'
#' @export
DJV <- function(X, k) {           
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  M1 <- sum(log(X[(n - k + 1):n]) - log(X[n - k])) / k
  M2 <- sum((log(X[(n - k + 1):n]) - log(X[n - k]))^2) / k
  2 * M1 / M2
}


#' DPR Estimator
#'
#' Computes the DPR estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param m An optional integer specifying the group size for dividing the sample. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The DPR estimator divides the sorted sample into groups of size \code{m}. For each group, the ratio of the second-to-last order statistic to the maximum order statistic is calculated. The last group includes the remaining observations. The final estimate is computed based on these ratios and their sum.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' DPR1(x)
#'
#' @export
DPR1 <- function(X, m) {          
  n <- length(X)
  if (missing(m)) m <- round(sqrt(n))
  r <- floor(n / m) - 1  # Last group is treated separately
  SS <- vector(length = (r + 1))  # Vector for M2/M1 ratios
  for (i in 1:r) {
    x1 <- sort(X[((i - 1) * m + 1):(i * m)])
    SS[i] <- x1[(m - 1)] / x1[m]
  }
  x <- sort(X[(m * r + 1):n])  # Last group formed from remaining observations
  n1 <- length(x)
  SS[r + 1] <- x[(n1 - 1)] / x[n1]
  ss <- sum(SS)
  ss / (r - ss)
}

#' DPR Estimator with Permutations
#'
#' Computes the DPR estimator for the tail index of a heavy-tailed distribution using multiple permutations of the data.
#'
#' @param X A numeric vector containing the data sample.
#' @param M An optional integer specifying the number of permutations to perform. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The DPR estimator applies the \code{DPR1} function to multiple permutations of the input data. For each permutation, the tail index is computed, and the final estimate is obtained as the mean of the estimates across all permutations. This approach enhances robustness to variations in the data.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' DPR(x)
#'
#' @export
DPR <- function(X, M) {   
  if (missing(M)) M <- round(sqrt(length(X)))  # Number of permutations
  V <- vector(length = M)
  n <- length(X)
  
  for (i in 1:M) {
    X1 <- sample(X, n, replace = FALSE)
    V[i] <- DPR1(X1)
  }
  mean(V)
}



#' Generalized DPR (GDPR) Estimator
#'
#' Computes the Generalized DPR (GDPR) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param m An optional integer specifying the group size for dividing the sample. Defaults to the square root of the sample size.
#' @param r An optional numeric parameter controlling the transformation. Defaults to \code{-0.1}. When \code{r = 0}, the logarithmic transformation is used.
#' @return A numeric value representing the estimated tail index.
#' @details The GDPR estimator generalizes the DPR1 estimator by introducing a parameter \code{r} to control the transformation applied to the group-based ratios (\code{SS}). 
#' - The input data is divided into groups of size \code{m}, with the last group formed from the remaining observations.
#' - For each group, the ratio of the second-to-last order statistic to the maximum order statistic is calculated.
#' - Depending on the value of \code{r}, either a logarithmic transformation (\code{-log(SS)}) or a generalized transformation (\code{(1 - SS^r) / r}) is applied.
#' - The final tail index is computed based on the mean of the transformed values.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GDPR1(x)
#'
#' @export
GDPR1 <- function(X, m, r) {          
  n <- length(X)
  if (missing(m)) m <- round(sqrt(n))
  if (missing(r)) r <- -0.1
  t <- floor(n / m) - 1  # Last group is treated separately
  SS <- vector(length = (t + 1))  # Vector for M2/M1 ratios
  for (i in 1:t) {
    x1 <- sort(X[((i - 1) * m + 1):(i * m)])
    SS[i] <- x1[(m - 1)] / x1[m]
  }
  x <- sort(X[(m * t + 1):n])  # Last group formed from remaining observations
  n1 <- length(x)
  SS[t + 1] <- x[(n1 - 1)] / x[n1]
  
  if (r == 0) {
    ZZ <- -log(SS)
  } else {
    ZZ <- (1 - SS^r) / r
  }
  
  pr <- mean(ZZ)
  (1 - r * pr) / pr
}

#' Generalized DPR (GDPR) Estimator with Permutations
#'
#' Computes the Generalized DPR (GDPR) estimator for the tail index of a heavy-tailed distribution using multiple permutations of the data.
#'
#' @param X A numeric vector containing the data sample.
#' @param M An optional integer specifying the number of permutations to perform. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The GDPR estimator applies the \code{GDPR1} function to multiple permutations of the input data. 
#' - For each permutation, the tail index is computed using the \code{GDPR1} method, which supports generalized transformations.
#' - The final estimate is obtained as the mean of the tail index estimates across all permutations.
#' This approach enhances robustness to variations in the data by averaging estimates derived from randomized orders of the sample.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GDPR(x)
#'
#' @export
GDPR <- function(X, M) {   
  if (missing(M)) M <- round(sqrt(length(X)))  # Number of permutations
  V <- vector(length = M)
  n <- length(X)
  
  for (i in 1:M) {
    X1 <- sample(X, n, replace = FALSE)
    V[i] <- GDPR1(X1)
  }
  mean(V)
}

#' Qi Estimator (2010)
#'
#' Computes the Qi (2010) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param m An optional integer specifying the group size for dividing the sample. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The Qi estimator divides the data into groups of size \code{m}. Within each group, the Hill estimator is applied to the top \code{k} order statistics, where \code{k = sqrt(m)}. The last group is formed from the remaining observations. The final tail index is computed as the mean of the group-specific Hill estimates.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' Qi1(x)
#'
#' @export
Qi1 <- function(X, m) {          
  n <- length(X)
  if (missing(m)) m <- round(sqrt(n))
  k <- round(sqrt(m))
  
  r <- floor(n / m) - 1  # Last group is treated separately
  HH <- vector(length = (r + 1))  # Vector for group-specific Hill estimates
  for (i in 1:r) {
    x1 <- X[((i - 1) * m + 1):(i * m)]
    HH[i] <- Hill(x1, k)
  }
  x1 <- X[(m * r + 1):n]  # Last group formed from remaining observations
  n1 <- length(x1)
  k <- round(sqrt(n1))
  HH[r + 1] <- Hill(x1, k)
  mean(HH)
}

#' Qi Estimator with Permutations
#'
#' Computes the Qi (2010) estimator for the tail index of a heavy-tailed distribution using multiple permutations of the data.
#'
#' @param X A numeric vector containing the data sample.
#' @param M An optional integer specifying the number of permutations to perform. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The Qi estimator applies the \code{Qi1} function to multiple permutations of the input data:
#' - For each permutation, the sample is divided into groups, and the Hill estimator is computed for the top-order statistics in each group.
#' - The final estimate is obtained as the mean of the tail index estimates across all permutations.
#' This approach improves robustness by incorporating randomness through permutations.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' Qi(x)
#'
#' @export
Qi <- function(X, M) {   
  if (missing(M)) M <- round(sqrt(length(X)))  # Number of permutations
  V <- vector(length = M)
  n <- length(X)
  
  for (i in 1:M) {
    X1 <- sample(X, n, replace = FALSE)
    V[i] <- Qi1(X1)
  }
  mean(V)
}



#' Vaiciulis Estimator (2012)
#'
#' Computes the Vaiciulis (2012) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param m An optional integer specifying the group size for dividing the sample. Defaults to the square root of the sample size.
#' @param l An optional integer specifying the maximum order of the power series. Defaults to 3.
#' @return A numeric value representing the estimated tail index.
#' @details The Vaiciulis estimator divides the input data into groups of size \code{m}. 
#' - Within each group, the logarithm of the ratio of the maximum order statistic to the second-to-last order statistic is computed and raised to successive powers up to \code{l+1}.
#' - The method calculates the averages of these transformed values across groups and applies a specific power series formula involving gamma functions to compute the tail index.
#' - The last group includes the remaining observations if the sample size is not divisible by \code{m}.
#'
#' This estimator involves higher-order corrections to improve the accuracy of tail index estimation.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' Vaiciulis(x)
#'
#' @export
Vaiciulis <- function(X, m, l) {          
  n <- length(X)
  if (missing(m)) m <- round(sqrt(n))
  if (missing(l)) l <- 3
  
  t <- floor(n / m) - 1  # Last group is treated separately
  FF <- matrix(nrow = (t + 1), ncol = (l + 1))  # Matrix for f(M2/M1)
  
  for (j in 1:(l + 1)) {
    for (i in 1:t) {
      x1 <- sort(X[((i - 1) * m + 1):(i * m)])
      FF[i, j] <- log(x1[m] / x1[(m - 1)])^j
    }
    
    x <- sort(X[(m * t + 1):n])  # Last group formed from remaining observations
    n1 <- length(x)
    FF[(t + 1), j] <- log(x[n1] / x[(n1 - 1)])^j
  }
  ss <- rowMeans(FF)
  
  D1 <- 0
  D2 <- 0
  for (i in 1:l) {
    D1 <- D1 + (-1)^(i + 1) * gamma(i + 1)^(-1) * ss[i]
  }
  for (i in 2:(l + 1)) {
    D2 <- D2 + (-1)^i * gamma(i + 1)^(-1) * ss[i]
  }
  D1 / D2
}


#' Vaiciulis Estimator (2009)
#'
#' Computes the Vaiciulis (2009) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param m An optional integer specifying the group size for dividing the sample. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The Vaiciulis (2009) estimator divides the data into overlapping groups of size \code{m}, computes the squared sums for consecutive groups, and calculates the ratio of the absolute difference to the sum of these squared sums. 
#' - The mean of these ratios (\code{ir}) is then used in a polynomial approximation to compute the tail index.
#' - The polynomial approximation incorporates trigonometric adjustment (\code{cos(pi * ir / 2)}) for accuracy.
#'
#' This estimator provides a robust approach for estimating the tail index by leveraging local squared-sum differences.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' Vaiciulis2009(x)
#'
#' @export
Vaiciulis2009 <- function(X, m) {          
  n <- length(X)
  if (missing(m)) m <- round(sqrt(n))
  
  IR <- vector(length = (n - 2 * m))
  
  for (i in 0:(n - 2 * m)) {
    R1 <- sum(X[(i + 1):(i + m)]^2)
    R2 <- sum(X[(i + m + 1):(i + 2 * m)]^2)
    IR[i] <- abs((R1 - R2)) / (R1 + R2)
  }
  ir <- mean(IR)
  (2.0024 - 0.4527 * ir + 0.4246 * ir^2 - 0.1386 * ir^3) * cos(pi * ir / 2)
}

#' Vaiciulis and Paulauskas Hill Estimator Modification
#'
#' Computes the Vaiciulis and Paulauskas (2013) modification of the Hill estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param r An optional numeric parameter controlling the transformation. Defaults to \code{-0.1}. When \code{r = 0}, the logarithmic transformation is used.
#' @return A numeric value representing the estimated tail index.
#' @details The Vaiciulis and Paulauskas modification of the Hill estimator introduces a parameter \code{r} for a generalized transformation of the scaled top-order statistics:
#' - When \code{r = 0}, the transformation is logarithmic (\code{log(Y)}).
#' - When \code{r != 0}, the transformation is \code{(Y^r - 1) / r}, where \code{Y} is the ratio of the top-order statistics to the threshold value.
#' - The tail index is calculated using the mean of the transformed values.
#'
#' This modification provides a more flexible framework for estimating the tail index by allowing different transformations based on the value of \code{r}.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' VP_H(x)
#'
#' @export
VP_H <- function(X, k, r) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(r)) r <- -0.1
  Y <- X[(n - k + 1):n] / X[n - k]
  
  if (r == 0) {
    ZZ <- log(Y)
  } else {
    ZZ <- (Y^r - 1) / r
  }
  
  lambda <- mean(ZZ)
  (1 + r * lambda) / lambda
}

#' Vaiciulis and Paulauskas Moment Estimator Modification
#'
#' Computes the Vaiciulis and Paulauskas (2013) modification of the moment estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param r An optional numeric parameter controlling the transformation. Defaults to \code{-0.1}. When \code{r = 0}, the logarithmic transformation is used.
#' @return A numeric value representing the estimated tail index.
#' @details The Vaiciulis and Paulauskas modification of the moment estimator introduces a parameter \code{r} for a generalized transformation of the scaled top-order statistics:
#' - When \code{r = 0}, the transformation is logarithmic (\code{log(Y)}).
#' - When \code{r != 0}, the transformation is \code{(Y^r - 1) / r}, where \code{Y} is the ratio of the top-order statistics to the threshold value.
#' - Higher-order moments (\code{M1, M12, M2}) are calculated and used to derive the tail index with a bias correction.
#'
#' This modification enhances the flexibility and robustness of the moment estimator by accommodating different transformations based on the value of \code{r}.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' VP_M(x)
#'
#' @export
VP_M <- function(X, k, r) {           
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(r)) r <- -0.1
  Y <- X[(n - k + 1):n] / X[n - k]
  
  if (r == 0) {
    M1 <- log(Y)
    M12 <- log(Y)
  } else {
    M1 <- (Y^r - 1) / r
    M12 <- (Y^(2 * r) - 1) / (2 * r)
  }
  M2 <- M1^2
  H1 <- mean(M1)
  H12 <- mean(M12)
  H2 <- mean(M2)
  
  (H1 + 1 - 0.5 * (1 - H1 * H12 / H2)^(-1))^(-1)
}

#' Vaiciulis and Paulauskas De Vries Estimator Modification
#'
#' Computes the Vaiciulis and Paulauskas (2013) modification of the De Vries estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param r An optional numeric parameter controlling the transformation. Defaults to \code{-0.1}. When \code{r = 0}, the logarithmic transformation is used.
#' @return A numeric value representing the estimated tail index.
#' @details The Vaiciulis and Paulauskas modification of the De Vries estimator introduces a parameter \code{r} for a generalized transformation of the scaled top-order statistics:
#' - When \code{r = 0}, the transformation is logarithmic (\code{log(Y)}).
#' - When \code{r != 0}, the transformation is \code{(Y^r - 1) / r}, where \code{Y} is the ratio of the top-order statistics to the threshold value.
#' - Higher-order moments (\code{M1, M12, M2}) are calculated and used to derive the tail index as \code{2 * H12 / H2}.
#'
#' This modification provides flexibility by allowing different transformations based on the value of \code{r}.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' VP_V(x)
#'
#' @export
VP_V <- function(X, k, r) {           
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(r)) r <- -0.1
  Y <- X[(n - k + 1):n] / X[n - k]
  
  if (r == 0) {
    M1 <- log(Y)
    M12 <- log(Y)
  } else {
    M1 <- (Y^r - 1) / r
    M12 <- (Y^(2 * r) - 1) / (2 * r)
  }
  M2 <- M1^2
  H1 <- mean(M1)
  H12 <- mean(M12)
  H2 <- mean(M2)
  
  2 * H12 / H2
}

#' Vaiciulis and Paulauskas Function \( G \)
#'
#' Computes the Vaiciulis and Paulauskas (2017) function \( G \), which is used in the context of tail index estimation for heavy-tailed distributions.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param r An optional numeric parameter controlling the power transformation. Defaults to \code{-0.1}.
#' @param u An optional numeric parameter controlling the logarithmic transformation. Defaults to \code{0}.
#' @return A numeric value representing the mean of the transformed top-order statistics based on the specified parameters \code{r} and \code{u}.
#' @details The function \( G \) transforms the scaled top-order statistics using a generalized formula:
#' - The scaled statistics are calculated as \( Y = X_{(n-k+1):n} / X_{(n-k)} \).
#' - The transformation is defined as \( g = Y^r \cdot (\log(Y))^u \), where \code{r} and \code{u} are parameters.
#' - The output is the mean of the transformed values (\( g \)).
#'
#' This function is used in various tail index estimation methodologies as a flexible tool for applying power and logarithmic transformations to order statistics.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' G(x)
#'
#' @export
G <- function(X, k, r, u) {           
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(r)) r <- -0.1
  if (missing(u)) u <- 0
  
  Y <- X[(n - k + 1):n] / X[n - k]
  g <- Y^r * log(Y)^u
  mean(g)
}

#' Vaiciulis and Paulauskas Gamma\(_2\) Estimator
#'
#' Computes the Vaiciulis and Paulauskas (2017) Gamma\(_2\) estimator (\(\gamma_2\)) for the tail index of a heavy-tailed distribution, as described on page 6 of their paper.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param r An optional numeric parameter controlling the transformation. Defaults to \code{-0.1}.
#' @return A numeric value representing the estimated tail index \(\gamma_2\).
#' @details This estimator relies on the transformation provided by the Vaiciulis and Paulauskas \( G \) function:
#' - The function \( G \) is called with \( r \) and logarithmic transformation parameter \( u = 1 \).
#' - The tail index \(\gamma_2\) is computed using the formula:
#'   \[
#'   \gamma_2 = \frac{D_2}{D_1}, \quad \text{where } D_1 = 2 \cdot G_1 \text{ and } D_2 = 2r \cdot G_1 + 1 + \sqrt{4r \cdot G_1 + 1}.
#'   \]
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' VP_2(x)
#'
#' @export
VP_2 <- function(X, k, r) {           
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(r)) r <- -0.1
  
  G1 <- G(X, k, r, 1)
  
  D1 <- 2 * G1
  D2 <- 2 * r * G1 + 1 + (4 * r * G1 + 1)^0.5
  D2 / D1
}

#' Vaiciulis and Paulauskas Gamma\(_3\) Estimator
#'
#' Computes the Vaiciulis and Paulauskas (2017) Gamma\(_3\) estimator (\(\gamma_3\)) for the tail index of a heavy-tailed distribution, as described on page 6 of their paper.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param r An optional numeric parameter controlling the transformation. Defaults to \code{-0.1}.
#' @return A numeric value representing the estimated tail index \(\gamma_3\).
#' @details This estimator uses the transformation provided by the Vaiciulis and Paulauskas \( G \) function:
#' - The function \( G \) is called with parameters \( r \) and \( u = 0 \) to compute \( G_0 \), and with \( u = 1 \) to compute \( G_1 \).
#' - The tail index \(\gamma_3\) is computed using the formula:
#'   \[
#'   \gamma_3 = \begin{cases} 
#'   \text{DJV}(X), & \text{if } r = 0, \\
#'   \frac{D_2}{D_1}, & \text{otherwise},
#'   \end{cases}
#'   \]
#'   where \( D_1 = r \cdot G_1 - G_0 + 1 \) and \( D_2 = r^2 \cdot G_1 \).
#'
#' If \( r = 0 \), the estimator defaults to the Danielsson, Jansen, and De Vries (DJV) estimator.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' VP_3(x)
#'
#' @export
VP_3 <- function(X, k, r) {           
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(r)) r <- -0.1
  
  G0 <- G(X, k, r, 0)
  G1 <- G(X, k, r, 1)
  
  D1 <- r * G1 - G0 + 1
  D2 <- (r^2) * G1
  
  if (r == 0) {
    a <- DJV(X)
  } else {
    a <- D2 / D1
  }
  a
}

#' Jurečková and Picek Estimator (2004)
#'
#' Computes the Jurečková and Picek (2004) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param m An optional integer specifying the group size for dividing the sample. Defaults to the square root of the sample size.
#' @param delta An optional numeric parameter controlling the scaling factor. Defaults to \code{0.05}.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator divides the input data into groups of size \code{m}, computes the maximum of each group, and derives the tail index using a quantile-based scaling formula:
#' - The input data is divided into \code{T} groups, and the maximum value of each group is stored in \code{SS}.
#' - The scaling factors \code{a1} and \code{a2} are calculated based on the group size (\code{m}) and parameter \code{delta}.
#' - The tail index is derived using the logarithm of \code{a1} and the quantile of the maxima vector \code{SS}.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' JP2004a(x)
#'
#' @export
JP2004a <- function(X, m, delta) {          
  n <- length(X)
  if (missing(m)) m <- round(sqrt(n))
  if (missing(delta)) delta <- 0.05
  
  T <- floor(n / m) - 1  # Last group is treated separately
  SS <- vector(length = (T + 1))  # Vector for maxima
  for (i in 1:T) {
    x1 <- sort(X[((i - 1) * m + 1):(i * m)])
    SS[i] <- x1[m]
  }
  x <- sort(X[(m * T + 1):n])  # Last group formed from remaining observations
  n1 <- length(x)
  SS[T + 1] <- x[n1]
  
  a1 <- m * T^(1 - delta)
  a2 <- (1 - T^(-(1 - delta)))
  
  log(a1) / log(quantile(SS, probs = a2))
}

#' Jurečková and Picek Estimator with Permutations (2004)
#'
#' Computes the Jurečková and Picek (2004) estimator for the tail index of a heavy-tailed distribution using multiple permutations of the data.
#'
#' @param X A numeric vector containing the data sample.
#' @param M An optional integer specifying the number of permutations to perform. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator applies the \code{JP2004a} function to multiple permutations of the input data:
#' - For each permutation, the sample is divided into groups, and the maxima of the groups are used to estimate the tail index.
#' - The final estimate is obtained as the mean of the tail index estimates across all permutations.
#'
#' Permutations enhance robustness by averaging the estimates derived from different random orders of the sample.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' JP2004(x)
#'
#' @export
JP2004 <- function(X, M) {   
  if (missing(M)) M <- round(sqrt(length(X)))  # Number of permutations
  V <- vector(length = M)
  n <- length(X)
  
  for (i in 1:M) {
    X1 <- sample(X, n, replace = FALSE)
    V[i] <- JP2004a(X1)
  }
  mean(V)
}


#' Brito, Cavalcante, and Freitas (BCF) Estimator
#'
#' Computes the Brito, Cavalcante, and Freitas (2016) estimator for the tail index of a heavy-tailed distribution using the method described in Equation 6 of their paper.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator computes the tail index by combining moment-based measures and log-rank statistics:
#' - Log-rank statistics (\(S_1, S_2\)) are calculated using the ranks of the data.
#' - Moment-based measures (\(M_1, M_2\)) are derived from the top-order statistics using the \code{MM} function.
#' - The tail index is calculated using the formula:
#'   \[
#'   \text{Tail Index} = \frac{1}{\sqrt{\frac{M_2 - M_1^2}{S_1 - S_2}}}.
#'   \]
#'
#' This method provides a robust estimate of the tail index by accounting for both moments and rank statistics.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' BCF(x)
#'
#' @export
BCF <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  
  ind <- 1:k
  S1 <- mean(log(n / ind)^2)
  S2 <- mean(log(n / ind))^2
  i <- S1 - S2
  
  M1 <- MM(X, k, 1)
  M2 <- MM(X, k, 2)
  1 / sqrt((M2 - M1^2) / i)
}



#' Brito, Cavalcante, and Freitas (BCF) Adjusted Estimator
#'
#' Computes the adjusted Brito, Cavalcante, and Freitas (2016) estimator for the tail index of a heavy-tailed distribution using the method described in Equation 12 of their paper.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use for the initial BCF estimation. Defaults to the square root of the sample size.
#' @param k2 An optional integer specifying the number of top-order statistics to use for the adjustment. Defaults to \(\max(0.9 \cdot n, k)\).
#' @return A numeric value representing the adjusted tail index estimate.
#' @details This estimator adjusts the BCF tail index by incorporating additional parameters:
#' - The initial tail index estimate (\(g\)) is computed using the \code{BCF} function.
#' - The adjustment parameters (\(\beta, \rho\)) are computed using the \code{beta} and \code{rho} functions, respectively.
#' - The adjusted tail index (\(g_2\)) is calculated using the formula:
#'   \[
#'   g_2 = g \cdot \left(1 - \frac{(n / k)^r \cdot \beta}{(1 - r)^2}\right),
#'   \]
#'   where \( g \) is the initial BCF estimate, \( \beta \) and \( \rho \) are adjustment parameters, and \( n \) and \( k \) are sample and group sizes.
#'
#' This method provides a bias-corrected estimate of the tail index.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' BCF2(x)
#'
#' @export
BCF2 <- function(X, k, k2) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(k2)) k2 <- round(max(0.9 * n, k))
  
  g <- BCF(X, k)^(-1)
  b <- beta(X, k2)
  r <- rho(X, k2)
  g2 <- g * (1 - (n / k)^r * b / (1 - r)^2)
  g2^(-1)
}

#' Brito, Cavalcante, and Freitas (BCF) Adjusted Estimator (Exponential Correction)
#'
#' Computes the adjusted Brito, Cavalcante, and Freitas (2016) estimator for the tail index of a heavy-tailed distribution using the method described in Equation 13 of their paper, with an exponential correction.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use for the initial BCF estimation. Defaults to the square root of the sample size.
#' @param k2 An optional integer specifying the number of top-order statistics to use for the adjustment. Defaults to \(\max(0.9 \cdot n, k)\).
#' @return A numeric value representing the adjusted tail index estimate.
#' @details This estimator adjusts the BCF tail index by incorporating an exponential correction:
#' - The initial tail index estimate (\(g\)) is computed using the \code{BCF} function.
#' - The adjustment parameters (\(\beta, \rho\)) are computed using the \code{beta} and \code{rho} functions, respectively.
#' - The adjusted tail index (\(g_2\)) is calculated using the formula:
#'   \[
#'   g_2 = g \cdot \exp\left(-\frac{(n / k)^r \cdot \beta}{(1 - r)^2}\right),
#'   \]
#'   where \( g \) is the initial BCF estimate, \( \beta \) and \( \rho \) are adjustment parameters, and \( n \) and \( k \) are sample and group sizes.
#'
#' This method provides a refined, bias-corrected estimate of the tail index using an exponential adjustment.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' BCF3(x)
#'
#' @export
BCF3 <- function(X, k, k2) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(k2)) k2 <- round(max(0.9 * n, k))
  
  g <- BCF(X, k)^(-1)
  b <- beta(X, k2)
  r <- rho(X, k2)
  g2 <- g * exp(-(n / k)^r * b / (1 - r)^2)
  g2^(-1)
}

#' Huisman, Koedijk, Kool, and Palm (HKKP) Estimator
#'
#' Computes the Huisman, Koedijk, Kool, and Palm (2001) estimator for the tail index of a heavy-tailed distribution using the method described in Equation 5 of their paper.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator refines the Hill estimator by regressing a sequence of Hill estimates over the number of top-order statistics:
#' - The Hill estimator is computed for each order statistic from 4 to \( k \).
#' - A linear regression is performed with the Hill estimates as the dependent variable and the corresponding indices as the independent variable.
#' - The tail index is extracted as the intercept of the regression line.
#'
#' This method smooths the variability of the Hill estimates, providing a more stable estimate of the tail index.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' HKKP(x)
#'
#' @export
HKKP <- function(X, k) {         
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  HH <- vector(length = (k - 3))
  ind <- c(4:k)
  for (i in 4:k) {
    HH[i - 3] <- Hill(X, i)
  }
  rg <- lm(HH ~ ind)
  rg$coefficients[1]
}


#' De Haan and Pereira Estimator
#'
#' Computes the De Haan and Pereira estimator for the tail index of a heavy-tailed distribution, as described in Equation (2) on page 40 of their paper.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator is based on the ratio of the sum of squares of the top \( k \) order statistics to the sum of squares of the remaining order statistics:
#' - The top-order statistics are used to compute \( D_1 \), which is the product of \( k \) and the square of the \( (n-k) \)-th order statistic.
#' - The lower-order statistics are used to compute \( D_2 \), which is the sum of the squares of the remaining observations.
#' - The intermediate parameter \( \beta \) is defined as \( \beta = D_1 / D_2 \).
#' - The tail index is derived as \( \frac{2}{\beta + 1} \).
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' dHP(x)
#'
#' @export
dHP <- function(X, k) {         
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  D1 <- k * X[n - k]^2
  ind <- 1:(n - k)
  D2 <- sum(X[ind]^2)
  beta <- D1 / D2
  2 / (beta + 1)
}

#' Fan Estimator
#'
#' Computes the Fan (2004) estimator for the tail index of a heavy-tailed distribution, as described in Equation (7) on page 19 of their paper.
#'
#' @param X A numeric vector containing the data sample.
#' @param m An optional integer specifying the subsample size for averaging. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator computes the tail index by averaging logarithmic transformations of subsamples:
#' - The sample is divided into \( N \) random subsamples of size \( m \), drawn with replacement.
#' - For each subsample, the logarithm of the sum of subsample values is divided by the logarithm of \( m \).
#' - The tail index is calculated as the reciprocal of the mean of these transformed values.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' Fan(x)
#'
#' @export
Fan <- function(X, m) {         
  X <- sort(X)
  n <- length(X)
  if (missing(m)) m <- round(sqrt(n))
  N <- 10000  # Number of subsamples
  Nm <- N * m
  XX <- matrix(sample(X, Nm, replace = TRUE), ncol = m, nrow = N)
  h <- log(rowSums(XX)) / log(m)
  mean(h)^(-1)
}

#' Danielsson, Jansen, and De Vries (DJV) Estimator
#'
#' Computes the Danielsson, Jansen, and De Vries (1996) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details The DJV estimator calculates the tail index using the ratio of the first and second moments of the logarithmic differences of the top-order statistics:
#' - \( M_1 \) is the first moment of the logarithmic differences of the top-order statistics.
#' - \( M_2 \) is the second moment of the logarithmic differences of the top-order statistics.
#' - The tail index is computed as:
#'   \[
#'   \text{Tail Index} = \frac{2 \cdot M_1}{M_2}.
#'   \]
#'
#' This estimator is particularly effective for heavy-tailed distributions, leveraging the log-transformed moments of the extreme order statistics.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' DJV(x)
#'
#' @export
DJV <- function(X, k) {           
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  M1 <- sum(log(X[(n - k + 1):n]) - log(X[n - k])) / k
  M2 <- sum((log(X[(n - k + 1):n]) - log(X[n - k]))^2) / k
  2 * M1 / M2
}





#' Segers Estimator (2001)
#'
#' Computes the Segers (2001) estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param p An optional numeric parameter controlling the power transformation. Defaults to \code{0.7}.
#' @return A numeric value representing the estimated tail index.
#' @details The Segers estimator computes the tail index using a power transformation of the ratios of order statistics:
#' - For the top \( k \) order statistics, the transformation is:
#'   \[
#'   M = \frac{1}{k} \sum_{i=1}^{k} \left( \frac{X_{(n-k)}}{X_{(n-k+i)}} \right)^p,
#'   \]
#'   where \( X_{(i)} \) denotes the \( i \)-th order statistic of the sorted sample.
#' - The tail index is then calculated as:
#'   \[
#'   \text{Tail Index} = \left( \frac{M^{-1} - 1}{p} \right)^{-1}.
#'   \]
#'
#' This method is particularly effective for analyzing heavy-tailed distributions with flexible power parameter \( p \).
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' Segers2001(x)
#'
#' @export
Segers2001 <- function(X, k, p) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(p)) p <- 0.7
  M <- sum((X[n - k] / X[(n - k + 1):n])^p) / k
  ((M^(-1) - 1) / p)^(-1)
}


#' Gomes and Martin (2001) Estimator - Type 1
#'
#' Computes the Gomes and Martin (2001) Type 1 estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param alpha An optional numeric parameter controlling the moment order. Defaults to \code{0.5}.
#' @return A numeric value representing the estimated tail index.
#' @details The Gomes and Martin Type 1 estimator computes the tail index using a combination of moments:
#' - The \( \alpha \)-th moment (\( M \)) and the first moment (\( H \)) of the logarithmic differences of the top-order statistics are computed using the \code{MM} function.
#' - The tail index is calculated as:
#'   \[
#'   \text{Tail Index} = \frac{\Gamma(\alpha + 1) \cdot H^{\alpha - 1}}{M},
#'   \]
#'   where \( \Gamma \) is the Gamma function, \( M \) is the \( \alpha \)-th moment, and \( H \) is the first moment.
#'
#' This estimator incorporates a flexible parameter \( \alpha \), allowing for tailored moment-based calculations.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GM1(x)
#'
#' @export
GM1 <- function(X, k, alpha) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(alpha)) alpha <- 0.5
  M <- MM(X, k, alpha)
  H <- MM(X, k, 1)
  (gamma(alpha + 1) * H^(alpha - 1)) / M
}

#' Gomes and Martin (2001) Estimator - Type 2
#'
#' Computes the Gomes and Martin (2001) Type 2 estimator for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param alpha An optional numeric parameter controlling the moment order. Defaults to \code{0.5}.
#' @return A numeric value representing the estimated tail index.
#' @details The Gomes and Martin Type 2 estimator computes the tail index using the \( \alpha \)-th moment of the logarithmic differences of the top-order statistics:
#' - The \( \alpha \)-th moment (\( M \)) is computed using the \code{MM} function.
#' - The tail index is calculated as:
#'   \[
#'   \text{Tail Index} = \left( \frac{\Gamma(\alpha + 1)}{M} \right)^{1 / \alpha},
#'   \]
#'   where \( \Gamma \) is the Gamma function and \( M \) is the \( \alpha \)-th moment.
#'
#' This estimator is a variation of the Type 1 estimator, providing a more direct approach based on the moment measure.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GM2(x)
#'
#' @export
GM2 <- function(X, k, alpha) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(alpha)) alpha <- 0.5
  M <- MM(X, k, alpha)
  (gamma(alpha + 1) / M)^(1 / alpha)
}

#' Caeiro et al. (2005) First Estimator
#'
#' Computes the first estimator from Caeiro et al. (2005) for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator refines the Hill estimator by incorporating auxiliary functions and adjustments:
#' - The auxiliary function \(\rho_1\) is computed using the \code{rho} function with \(\tau = 0\).
#' - The auxiliary function \(\beta_1\) is computed using the \code{beta} function.
#' - The first moment of the logarithmic differences of the top-order statistics (\( H \)) is computed using the \code{MM} function.
#' - The tail index is calculated using:
#'   \[
#'   \text{Tail Index} = \left[ H \cdot \left( 1 - \frac{\beta_1}{1 - \rho_1} \cdot \left( \frac{n}{k} \right)^{\rho_1} \right) \right]^{-1}.
#'   \]
#'
#' This estimator combines the Hill estimator with auxiliary functions to provide a more robust estimate of the tail index.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' CGP1(x)
#'
#' @export
CGP1 <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  rho1 <- rho(X, k, 0)
  beta1 <- beta(X, k)
  H <- MM(X, k, 1)
  (H * (1 - beta1 / (1 - rho1) * (n / k)^rho1))^(-1)
}




#' Caeiro et al. (2005) Second Estimator
#'
#' Computes the second estimator from Caeiro et al. (2005) for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator refines the Hill estimator by incorporating exponential adjustments based on auxiliary functions:
#' - The auxiliary function \(\rho_1\) is computed using the \code{rho} function with \(\tau = 0\).
#' - The auxiliary function \(\beta_1\) is computed using the \code{beta} function.
#' - The first moment of the logarithmic differences of the top-order statistics (\( H \)) is computed using the \code{MM} function.
#' - The tail index is calculated using:
#'   \[
#'   \text{Tail Index} = \left[ H \cdot \exp\left(-\frac{\beta_1}{1 - \rho_1} \cdot \left( \frac{n}{k} \right)^{\rho_1}\right) \right]^{-1}.
#'   \]
#'
#' This estimator combines the Hill estimator with exponential adjustments based on auxiliary functions for increased robustness.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' CGP2(x)
#'
#' @export
CGP2 <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  rho1 <- rho(X, k, 0)
  beta1 <- beta(X, k)
  H <- MM(X, k, 1)
  (H * exp(-beta1 / (1 - rho1) * (n / k)^rho1))^(-1)
}

#' Gomes et al. (2000) Jackknife Estimator - Type 1
#'
#' Computes the first Jackknife estimator from Gomes et al. (2000) for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This Jackknife estimator combines the Hill estimator and the Danielsson, Jansen, and De Vries (DJV) estimator:
#' - The inverse of the Hill estimator (\( g_1 \)) is computed using the \code{Hill} function.
#' - The inverse of the DJV estimator (\( g_2 \)) is computed using the \code{DJV} function.
#' - The tail index is calculated as:
#'   \[
#'   \text{Tail Index} = \left[ 2 \cdot g_2 - g_1 \right]^{-1}.
#'   \]
#'
#' This method uses a Jackknife approach to reduce bias by combining two different estimators.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GMN1(x)
#'
#' @export
GMN1 <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  g1 <- (Hill(X, k))^(-1)
  g2 <- (DJV(X, k))^(-1)
  (2 * g2 - g1)^(-1)
}

#' Gomes et al. (2000) Jackknife Estimator - Type 2
#'
#' Computes the second Jackknife estimator from Gomes et al. (2000) for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This Jackknife estimator combines the Hill estimator and the Gomes-Martin Type 2 estimator:
#' - The inverse of the Hill estimator (\( g_1 \)) is computed using the \code{Hill} function.
#' - The inverse of the Gomes-Martin Type 2 estimator (\( g_3 \)) is computed using the \code{GM2} function with \code{alpha = 2}.
#' - The tail index is calculated as:
#'   \[
#'   \text{Tail Index} = 4 \cdot g_3 - 3 \cdot g_1.
#'   \]
#'
#' This method uses a Jackknife approach to reduce bias by combining the Hill and Gomes-Martin estimators.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GMN2(x)
#'
#' @export
GMN2 <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  g1 <- (Hill(X, k))^(-1)
  g3 <- (GM2(X, k, 2))^(-1)
  4 * g3 - 3 * g1
}


#' Gomes et al. (2000) Jackknife Estimator - Type 3
#'
#' Computes the third Jackknife estimator from Gomes et al. (2000) for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This Jackknife estimator combines the Danielsson, Jansen, and De Vries (DJV) estimator and the Gomes-Martin Type 2 estimator:
#' - The inverse of the DJV estimator (\( g_2 \)) is computed using the \code{DJV} function.
#' - The inverse of the Gomes-Martin Type 2 estimator (\( g_3 \)) is computed using the \code{GM2} function with \code{alpha = 2}.
#' - The tail index is calculated as:
#'   \[
#'   \text{Tail Index} = \left[ 3 \cdot g_2 - 2 \cdot g_3 \right]^{-1}.
#'   \]
#'
#' This method uses a Jackknife approach to reduce bias by combining the DJV and Gomes-Martin estimators.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GMN3(x)
#'
#' @export
GMN3 <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  g2 <- (DJV(X, k))^(-1)
  g3 <- (GM2(X, k, 2))^(-1)
  (3 * g2 - 2 * g3)^(-1)
}

#' Gomes et al. (2000) Jackknife Estimator - Type 4
#'
#' Computes the fourth Jackknife estimator from Gomes et al. (2000) for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This Jackknife estimator combines two Hill estimators using different subsets of the top-order statistics:
#' - The inverse of the Hill estimator (\( g_1 \)) is computed using \code{k / 2} order statistics.
#' - The inverse of the Hill estimator (\( g_2 \)) is computed using \( k \) order statistics.
#' - The tail index is calculated as:
#'   \[
#'   \text{Tail Index} = \left[ 2 \cdot g_1 - g_2 \right]^{-1}.
#'   \]
#'
#' This method uses a Jackknife approach to reduce bias by combining Hill estimators computed from two different subsets of the data.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GMN4(x)
#'
#' @export
GMN4 <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  g1 <- (Hill(X, round(k / 2)))^(-1)
  g2 <- (Hill(X, k))^(-1)
  (2 * g1 - g2)^(-1)
}



#' Gomes and Martin (2002) Estimator - Type 1
#'
#' Computes the first estimator from Gomes and Martin (2002) for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param tau An optional numeric parameter for auxiliary adjustment. Defaults to \code{0}.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator refines the tail index calculation by incorporating the auxiliary function \(\rho_1\), the DJV estimator, and the Gomes-Martin Type 2 estimator:
#' - The auxiliary function \(\rho_1\) is computed using the \code{rho} function with \(\tau = 0\).
#' - The inverse of the DJV estimator (\( g_2 \)) is computed using the \code{DJV} function.
#' - The inverse of the Gomes-Martin Type 2 estimator (\( g_3 \)) is computed using the \code{GM2} function with \code{alpha = 2}.
#' - The tail index is calculated as:
#'   \[
#'   \text{Tail Index} = \frac{\rho_1}{-(2 - \rho_1) \cdot g_2 + 2 \cdot g_3}.
#'   \]
#'
#' This estimator uses the auxiliary function \(\rho_1\) to account for additional structure in the data, improving bias and variance properties.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GM3(x)
#'
#' @export
GM3 <- function(X, k, tau) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(tau)) tau <- 0
  rho1 <- rho(X, k, 0)
  g2 <- (DJV(X, k))^(-1)
  g3 <- (GM2(X, k, 2))^(-1)
  rho1 / (-(2 - rho1) * g2 + 2 * g3)
}

#' Gomes and Martin (2002) Estimator - Type 2
#'
#' Computes the second estimator from Gomes and Martin (2002) for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param tau An optional numeric parameter for auxiliary adjustment. Defaults to \code{0}.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator refines the tail index calculation by incorporating the auxiliary function \(\rho_1\) and two Hill estimators:
#' - The auxiliary function \(\rho_1\) is computed using the \code{rho} function with \(\tau = 0\).
#' - The inverse of the Hill estimator (\( g_2 \)) is computed using \code{k/2} order statistics.
#' - The inverse of the Hill estimator (\( g_3 \)) is computed using \( k \) order statistics.
#' - The tail index is calculated as:
#'   \[
#'   \text{Tail Index} = \frac{1 - 2^{-\rho_1}}{g_3 - 2^{-\rho_1} \cdot g_2}.
#'   \]
#'
#' This method uses the auxiliary function \(\rho_1\) to adjust the Hill estimators, reducing bias and improving accuracy.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GM4(x)
#'
#' @export
GM4 <- function(X, k, tau) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(tau)) tau <- 0
  rho1 <- rho(X, k, 0)
  g2 <- (Hill(X, round(k / 2)))^(-1)
  g3 <- (Hill(X, k))^(-1)
  (1 - 2^(-rho1)) / (g3 - 2^(-rho1) * g2)
}

#' Gomes et al. (2002) Jackknife Estimator - Type 1
#'
#' Computes the first Jackknife estimator from Gomes et al. (2002) for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator refines the tail index using a leave-one-out Jackknife approach:
#' - For each observation \( i \), the inverse Hill estimator (\( \gamma_{n1,i} \)) is computed using the sample with the \( i \)-th observation removed.
#' - The overall inverse Hill estimator (\( g \)) is computed for the full sample.
#' - The mean of the leave-one-out estimates (\( \bar{\gamma}_{n1} \)) is calculated.
#' - The tail index is then computed using:
#'   \[
#'   \text{Tail Index} = \left[ \frac{n \cdot g - (n - 1) \cdot \bar{\gamma}_{n1}}{1} \right]^{-1}.
#'   \]
#'
#' This method reduces bias by incorporating the Jackknife technique.
#'
#' @examples
#' set.seed(123)
#' n = 100
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GMNj1(x)
#'
#' @export
GMNj1 <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  
  gamma_n1 <- vector(length = n)
  for (i in 1:n) {
    ind <- rep(TRUE, n)
    ind[i] <- FALSE  # i-th element is removed
    X1 <- X[ind]
    gamma_n1[i] <- (Hill(X1, k))^(-1)
  }
  
  g <- (Hill(X, k))^(-1)
  bar_g <- mean(gamma_n1)
  (n * g - (n - 1) * bar_g)^(-1)
}

#' Gomes et al. (2002) Jackknife Estimator - Type 2
#'
#' Computes the second Jackknife estimator from Gomes et al. (2002) for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the largest even integer less than or equal to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator refines the tail index using a Jackknife approach and two Hill estimators:
#' - The inverse Hill estimator (\( g_1 \)) is computed using \( k \) order statistics.
#' - The inverse Hill estimator (\( g_2 \)) is computed using \( k/2 \) order statistics.
#' - The weight \( w \) is calculated as:
#'   \[
#'   w = \frac{\log(1 - k / n)}{\log(1 - k / (2 \cdot n))}.
#'   \]
#' - The tail index is then computed as:
#'   \[
#'   \text{Tail Index} = \frac{1 - w}{g_1 - g_2 \cdot w}.
#'   \]
#'
#' This method reduces bias by combining two Hill estimators with a weighting factor \( w \).
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GMNj2(x)
#'
#' @export
GMNj2 <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n) / 2) * 2
  g2 <- (Hill(X, round(k / 2)))^(-1)
  g1 <- (Hill(X, k))^(-1)
  w <- log(1 - k / n) / log(1 - k / (2 * n))
  (1 - w) / (g1 - g2 * w)
}

#' Gomes et al. (2002) Jackknife Estimator - Type 3
#'
#' Computes the third Jackknife estimator from Gomes et al. (2002) for the tail index of a heavy-tailed distribution.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the largest even integer less than or equal to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator refines the tail index using a Jackknife approach and two Hill estimators:
#' - The inverse Hill estimator (\( g_1 \)) is computed using \( k \) order statistics.
#' - The inverse Hill estimator (\( g_2 \)) is computed using \( k/2 \) order statistics.
#' - The tail index is then computed as:
#'   \[
#'   \text{Tail Index} = \frac{1 + k / n}{g_2 \cdot (2 + k / n) - g_1}.
#'   \]
#'
#' This method combines two Hill estimators with an adjustment term \( k / n \) to improve bias and variance properties.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GMNj3(x)
#'
#' @export
GMNj3 <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n) / 2) * 2
  g2 <- (Hill(X, round(k / 2)))^(-1)
  g1 <- (Hill(X, k))^(-1)
  (1 + k / n) / (g2 * (2 + k / n) - g1)
}

#' Gomes et al. (2005) Generalized Estimator
#'
#' Computes the generalized tail index estimator from Gomes et al. (2005), as described on page 3 of their paper.
#'
#' @param X A numeric vector containing the data sample.
#' @param s An optional numeric parameter controlling the order of the estimator. Defaults to \code{2} (with \code{s >= 1}).
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator generalizes the tail index calculation by introducing a parameter \( s \) that adjusts the influence of the order statistics:
#' - The tail index is calculated using:
#'   \[
#'   \text{Tail Index} = \left[ \frac{s^2}{k^s} \cdot \sum_{i=1}^{k} i^{s-1} \cdot \log\left( \frac{X_{(n-k+i)}}{X_{(n-k)}} \right) \right]^{-1},
#'   \]
#'   where \( X_{(i)} \) represents the \( i \)-th order statistic of the sorted sample, \( k \) is the number of top-order statistics, and \( s \) is the order parameter.
#'
#' This method allows for greater flexibility in estimating the tail index by adjusting the parameter \( s \), with \( s \geq 1 \).
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GPM(x)
#'
#' @export
GPM <- function(X, s, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(s)) s <- 2  # s >= 1
  ind <- c(1:k)
  (s^2 * sum(ind^(s - 1) * log(rev(X[(n - k + 1):n] / X[n - k]))) / k^s)^(-1)
}

#' Gomes et al. (2005) Jackknife Estimator
#'
#' Computes the Jackknife estimator from Gomes et al. (2005), as described on page 3 of their paper.
#'
#' @param X A numeric vector containing the data sample.
#' @param s An optional numeric parameter controlling the order of the estimator. Defaults to \code{2} (with \code{s >= 1}).
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator refines the generalized tail index estimator by using a Jackknife adjustment:
#' - The generalized estimators \( g_1 \) and \( g_2 \) are computed using \code{GPM} with order parameters \( 1 \) and \( s \), respectively.
#' - The auxiliary function \(\rho_1\) is computed using the \code{rho} function with \code{tau = 1}. It is recommended to use \(\rho_1 = -1\) for improved performance.
#' - The tail index is calculated as:
#'   \[
#'   \text{Tail Index} = \left[ -\frac{s \cdot (1 - \rho_1)}{\rho_1 \cdot (s - 1)} \cdot \left( g_1 - \frac{s - \rho_1}{s \cdot (1 - \rho_1)} \cdot g_2 \right) \right]^{-1}.
#'   \]
#'
#' This method combines the generalized estimators with a Jackknife adjustment for bias reduction.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GPMj(x)
#'
#' @export
GPMj <- function(X, s, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(s)) s <- 2  # s >= 1
  g1 <- GPM(X, 1, k)^(-1)
  g2 <- GPM(X, s, k)^(-1)
  rho1 <- rho(X, k, 1)  # it is better to use -1
  # rho1 <- -1
  (-s * (1 - rho1) / (rho1 * (s - 1)) * (g1 - (s - rho1) / (s * (1 - rho1)) * g2))^(-1)
}

#' Gomes et al. (2005) Best Linear Estimator
#'
#' Computes the best linear estimator for the tail index from Gomes et al. (2005), as described on page 7 of their paper with \(\rho = -1\).
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator uses a weighted combination of Hill estimators to calculate the tail index:
#' - For each \( i \) from 1 to \( k \), the inverse Hill estimator (\( H_i \)) is computed for the top \( i \) order statistics.
#' - The weighted sum of the Hill estimators is calculated as:
#'   \[
#'   \text{Tail Index} = \left[ \frac{6}{k^2 - 1} \sum_{i=1}^{k-1} H_i \cdot i - H_k \cdot \frac{2k - 1}{k + 1} \right]^{-1}.
#'   \]
#'
#' This method improves upon traditional Hill estimators by leveraging a linear weighting scheme for reduced bias, assuming \(\rho = -1\).
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GFM(x)
#'
#' @export
GFM <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  H <- vector(length = k)
  for (i in 1:k) {
    H[i] <- Hill(X, i)^(-1)
  }
  ind <- c(1:(k - 1))
  (sum(H[1:(k - 1)] * ind) * 6 / (k^2 - 1) - H[k] * (2 * k - 1) / (k + 1))^(-1)
}

#' Gomes et al. (2005) Best Linear Estimator (General \(\rho\))
#'
#' Computes the best linear estimator for the tail index from Gomes et al. (2005), as described on page 17 of their paper, accommodating a general \(\rho\).
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator uses a weighted combination of log-transformed order statistics, with weights adjusted for a general auxiliary parameter \(\rho\):
#' - The auxiliary function \(\rho\) is computed using the \code{rho} function.
#' - Weights are calculated as:
#'   \[
#'   a_i = \frac{(1 - \rho)^2}{\rho^2 \cdot k} \cdot \left( 1 - \frac{k \cdot (1 - 2 \cdot \rho)}{1 - \rho} \cdot \left( \left( \frac{i}{k} \right)^{1 - \rho} - \left( \frac{i - 1}{k} \right)^{1 - \rho} \right) \right),
#'   \]
#'   for \( i = 1, \ldots, k \).
#' - The \( (k+1) \)-th weight is:
#'   \[
#'   a_{k+1} = -\frac{1 - \rho}{\rho}.
#'   \]
#' - The tail index is calculated as:
#'   \[
#'   \text{Tail Index} = \left[ \sum_{i=1}^{k+1} a_i \cdot \log(X_{(n-k+i)}) \right]^{-1}.
#'   \]
#'
#' This estimator provides flexibility in modeling tail behavior by incorporating the auxiliary parameter \(\rho\).
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GFM2(x)
#'
#' @export
GFM2 <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  
  a <- vector(length = (k + 1))  # Vector for weights
  
  ind <- c(1:k)
  r <- rho(X)
  a1 <- ((1 - r) / r)^2 / k
  a2 <- k * (1 - 2 * r) / (1 - r)
  a[1:k] <- a1 * (1 - a2 * ((ind / k)^(1 - r) - ((ind - 1) / k)^(1 - r)))
  a[k + 1] <- -(1 - r) / r
  1 / (sum(a * log(X[n:(n - k)])))
}

#' Caeiro and Gomes (2002) Estimator
#'
#' Computes the tail index estimator from Caeiro and Gomes (2002).
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param alpha An optional numeric parameter for the moment order. If not provided, it is estimated from the auxiliary function \(\rho\). Defaults to \code{NULL}.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator combines moments with the auxiliary parameter \(\alpha\) to refine the tail index:
#' - If \(\alpha\) is not provided, it is estimated as:
#'   \[
#'   \alpha = \max\left(1, \frac{-\log(1 - r - \sqrt{(1 - r)^2 - 1})}{\log(1 - r)}\right),
#'   \]
#'   where \( r \) is computed using the \code{rho} function.
#' - The tail index is then calculated using:
#'   \[
#'   \text{Tail Index} = \frac{1}{\left( \frac{\Gamma(\alpha)}{\text{MM}(X, k, \alpha - 1)} \right) \cdot \sqrt{\frac{\text{MM}(X, k, 2\alpha)}{\Gamma(2\alpha + 1)}}}.
#'   \]
#'
#' This method combines higher-order moments with auxiliary parameters to improve the accuracy of the tail index estimation.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' CG(x)
#'
#' @export
CG <- function(X, k, alpha) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  
  if (missing(alpha)) {
    r <- rho(X, k)
    alpha <- -log(1 - r - sqrt((1 - r)^2 - 1)) / log(1 - r)
    alpha <- max(alpha, 1)
  }
  D1 <- gamma(alpha) / MM(X, k, (alpha - 1))
  D2 <- MM(X, k, (2 * alpha) / gamma(2 * alpha + 1))
  1 / (D1 * sqrt(D2))
}

#' Gomes, Miranda, and Viseu (Statistica Neerlandica) Estimator
#'
#' Computes the tail index estimator from Gomes, Miranda, and Viseu (Statistica Neerlandica, equation 5, page 3).
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param alpha An optional numeric parameter controlling the moment order. Defaults to \code{1.5}.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator refines the tail index by combining a weighted sum of differences between log-transformed order statistics:
#' - The weights depend on the moment order parameter \(\alpha\) and the rank of the order statistics.
#' - The tail index is calculated as:
#'   \[
#'   \text{Tail Index} = \left[ \frac{\alpha}{k} \cdot \sum_{i=1}^{k} \left( \frac{i}{k} \right)^{\alpha - 1} \cdot i \cdot \left( \log(X_{(n-i+1)}) - \log(X_{(n-i)}) \right) \right]^{-1}.
#'   \]
#'
#' This method uses the parameter \(\alpha\) to adjust the influence of order statistics, with \(\alpha = 1.5\) as a common default.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GMV1(x)
#'
#' @export
GMV1 <- function(X, k, alpha) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(alpha)) alpha <- 1.5
  
  ind <- c(1:k)
  U <- ind * (log(X[n - ind + 1]) - log(X[n - ind]))
  (alpha / k * sum((ind / k)^(alpha - 1) * U))^(-1)
}

#' Gomes, Miranda, and Viseu (Statistica Neerlandica) Estimator - Type 2
#'
#' Computes the second tail index estimator from Gomes, Miranda, and Viseu (Statistica Neerlandica, equation 6, page 3).
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param alpha An optional numeric parameter controlling the moment order. Defaults to \code{1.5}.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator refines the tail index by combining a weighted sum of differences between log-transformed order statistics, incorporating logarithmic adjustments:
#' - The weights depend on the moment order parameter \(\alpha\), the rank of the order statistics, and a logarithmic scaling factor.
#' - The tail index is calculated as:
#'   \[
#'   \text{Tail Index} = \left[ -\frac{\alpha^2}{k} \cdot \sum_{i=1}^{k} \left( \frac{i}{k} \right)^{\alpha - 1} \cdot \log\left(\frac{i}{k}\right) \cdot i \cdot \left( \log(X_{(n-i+1)}) - \log(X_{(n-i)}) \right) \right]^{-1}.
#'   \]
#'
#' This method uses the parameter \(\alpha\) to adjust the influence of order statistics, with \(\alpha = 1.5\) as a common default, and logarithmic scaling to further refine the weights.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GMV2(x)
#'
#' @export
GMV2 <- function(X, k, alpha) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(alpha)) alpha <- 1.5
  
  ind <- c(1:k)
  U <- ind * (log(X[n - ind + 1]) - log(X[n - ind]))
  (-alpha^2 / k * sum((ind / k)^(alpha - 1) * log(ind / k) * U))^(-1)
}


#' Gomes, Miranda, and Viseu (Statistica Neerlandica) Jackknife Estimator
#'
#' Computes the Jackknife tail index estimator from Gomes, Miranda, and Viseu (Statistica Neerlandica, equation 26, page 11).
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This Jackknife estimator refines the tail index by combining generalized estimators with an auxiliary parameter \(\rho\):
#' - The auxiliary parameter \(\rho\) is calculated using the \code{rho} function, with \( k \cdot 1.1 \) or \( 0.9 \cdot n \) order statistics, whichever is smaller.
#' - Two generalized estimators are computed:
#'   - \( g_1 \): The generalized estimator with \(\alpha = 1\).
#'   - \( g_2 \): The generalized estimator with \(\alpha = 1 - \rho\).
#' - The tail index is then calculated as:
#'   \[
#'   \text{Tail Index} = \frac{\rho^2}{(1 - \rho)^2 \cdot g_1 - (1 - 2 \cdot \rho) \cdot g_2}.
#'   \]
#'
#' This method reduces bias by combining estimators with a weighting scheme influenced by \(\rho\).
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GMV3(x)
#'
#' @export
GMV3 <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  
  r <- rho(X, min((k * 1.1), (0.9 * n)))  # k shall be higher for rho calculation
  g1 <- GMV1(X, k, 1)^(-1)
  g2 <- GMV1(X, k, (1 - r))^(-1)
  
  r^2 / ((1 - r)^2 * g1 - (1 - 2 * r) * g2)
}

#' Gomes, Miranda, and Viseu (Statistica Neerlandica) Jackknife Estimator - Type 2
#'
#' Computes the Jackknife tail index estimator from Gomes, Miranda, and Viseu (Statistica Neerlandica, equation 27, page 11).
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param eps An optional numeric parameter specifying the maximum error tolerance for calculating the optimal \(\alpha\). Defaults to \code{1e-8}.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator refines the tail index by iteratively solving for the optimal parameter \(\alpha\) and using it in a weighted Jackknife formula:
#' - The auxiliary parameter \(\rho\) is calculated using the \code{rho} function with \( k \cdot 1.1 \) or \( 0.9 \cdot n \) order statistics, whichever is smaller.
#' - The optimal \(\alpha\) is determined by solving the equation:
#'   \[
#'   T = 3\alpha^3 - 5\alpha^2 + \alpha \cdot (r^2 - r + 3) - (2r^2 - 2r + 1) = 0,
#'   \]
#'   using a binary search approach with a tolerance of \code{eps}.
#' - Two generalized estimators are computed:
#'   - \( g_1 \): The generalized estimator using \(\text{GMV1}\) with the optimal \(\alpha\).
#'   - \( g_2 \): The generalized estimator using \(\text{GMV2}\) with the optimal \(\alpha\).
#' - The tail index is then calculated as:
#'   \[
#'   \text{Tail Index} = \frac{r}{\alpha \cdot g_1 - (\alpha - r) \cdot g_2}.
#'   \]
#'
#' This method combines iterative optimization and weighting schemes to improve bias reduction.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GMV4(x)
#'
#' @export
GMV4 <- function(X, k, eps) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(eps)) eps <- 1e-8  # Default maximal error for alpha calculation
  
  r <- rho(X, min((k * 1.1), (0.9 * n)))  # k shall be higher for rho calculation
  # ------------------------------------------------------------ Calculation for optimal alpha
  a_min <- 1          
  a_max <- 100        # Maximal value of alpha
  
  repeat {
    alpha <- (a_min + a_max) / 2
    T <- 3 * alpha^3 - 5 * alpha^2 + alpha * (r^2 - r + 3) - (2 * r^2 - 2 * r + 1)
    if (T < 0) {
      a_min <- alpha
    } else {
      a_max <- alpha
    }
    if (abs(a_max - a_min) < eps) break
  }
  
  g1 <- GMV1(X, k, alpha)^(-1)
  g2 <- GMV2(X, k, alpha)^(-1)
  r / (alpha * g1 - (alpha - r) * g2)
}

#' Gomes, Miranda, and Viseu (Statistica Neerlandica) Estimator - Type 5 (\(\rho = -1\))
#'
#' Computes the tail index estimator from Gomes, Miranda, and Viseu (Statistica Neerlandica, page 15, \(\rho = -1\)).
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator refines the tail index by combining two generalized estimators:
#' - Two generalized estimators are computed:
#'   - \( g_1 \): The generalized estimator using \(\text{GMV1}\) with \(\alpha = 1\).
#'   - \( g_2 \): The generalized estimator using \(\text{GMV1}\) with \(\alpha = 2\).
#' - The tail index is calculated as:
#'   \[
#'   \text{Tail Index} = \left[ 4 \cdot g_1 - 3 \cdot g_2 \right]^{-1}.
#'   \]
#'
#' This method reduces bias by combining the two estimators in a weighted linear fashion, assuming \(\rho = -1\).
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GMV5(x)
#'
#' @export
GMV5 <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  
  g1 <- GMV1(X, k, 1)^(-1)
  g2 <- GMV1(X, k, 2)^(-1)
  
  (4 * g1 - 3 * g2)^(-1)
}

#' Beirlant, Figueiredo, Gomes, and Vandewalle (2008) Estimator
#'
#' Computes the improved reduced-bias tail index estimator from Beirlant, Figueiredo, Gomes, and Vandewalle (2008).
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator refines the tail index by combining reduced-bias methods with auxiliary parameter adjustments:
#' - The auxiliary parameter \(\rho\) is computed using the \code{rho} function with \(k_1 = n^{0.9}\) order statistics.
#' - Several quantities are calculated for bias reduction:
#'   - \( D_0, D_1, D_2 \): Weighted sums of log differences for \(k_1\) top-order statistics.
#'   - \(\beta\): Bias adjustment parameter calculated from the weighted sums.
#' - The optimal number of order statistics \(k_0\) is determined as:
#'   \[
#'   k_0 = \left( \frac{(1 - 2r) \cdot n^{-2r}}{-2r \cdot \beta^2} \right)^{1/(1-2r)}.
#'   \]
#'   \( k_0 \) is constrained to be within \([k/10, k_1]\).
#' - The tail index is calculated as:
#'   \[
#'   \text{Tail Index} = \left[ H - \beta \cdot \left(\frac{n}{k}\right)^r \cdot D_0 \right]^{-1},
#'   \]
#'   where \( H \) is the Hill estimator.
#'
#' This method reduces bias in the tail index estimation by leveraging auxiliary parameters and optimal order statistic selection.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' BFGV(x)
#'
#' @export
BFGV <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  
  k1 <- round(n^0.9)
  
  ind1 <- c(1:k1)
  U1 <- ind1 * (log(X[n - ind1 + 1]) - log(X[n - ind1]))
  
  r <- rho(X, k1, 0)
  d1 <- mean((ind1 / (k1 + 1))^(-r))  # for beta calculation
  
  D2 <- mean((ind1 / (k1 + 1))^(-2 * r) * U1)
  D1 <- mean((ind1 / (k1 + 1))^(-r) * U1)
  D0 <- mean(U1)
  beta <- ((k + 1) / (n + 1))^r * (d1 * D0 - D1) / (d1 * D1 - D2)
  
  k0 <- ((1 - 2 * r) * n^(-2 * r) / (-2 * r * beta^2))^(1 / (1 - 2 * r))
  k0 <- max(min(k0, k1), abs(k / 10))
  ind0 <- c(1:k0)  # for D calculation
  U0 <- ind0 * (log(X[n - ind0 + 1]) - log(X[n - ind0]))
  D0 <- mean((ind0 / (k0 + 1))^(-2 * r) * U0)
  
  H <- Hill(X, k)^(-1)
  (H - beta * (n / k)^r * D0)^(-1)
}


#' Fraga Alves, Gomes, de Haan, Neves (2009) Mixed Moment Estimator
#'
#' Computes the mixed moment tail index estimator from Fraga Alves, Gomes, de Haan, and Neves (2009), as described on page 3.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator refines the tail index using mixed moments and an auxiliary bias adjustment:
#' - The first moment \( M_1 \) is computed using the \code{MM} function with \( \alpha = 1 \).
#' - A bias correction term \( L_1 \) is calculated as:
#'   \[
#'   L_1 = \text{mean}\left( 1 - \frac{X_{(n-k)}}{X_{(n-i+1)}} \right), \, i = 1, \ldots, k.
#'   \]
#' - The auxiliary parameter \( \phi \) is then calculated as:
#'   \[
#'   \phi = \frac{M_1 - L_1}{L_1^2}.
#'   \]
#' - The tail index is computed as:
#'   \[
#'   \text{Tail Index} = \frac{1 + 2 \cdot \min(\phi - 1, 0)}{\phi - 1}.
#'   \]
#' - The adjustment ensures non-positive bias corrections to improve accuracy.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' AGHN(x)
#'
#' @export
AGHN <- function(X, k) {         
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  
  M1 <- MM(X, k, 1)
  ind <- c(1:k)
  L1 <- mean(1 - X[n - k] / X[n - ind + 1])
  phi <- (M1 - L1) / L1^2
  
  (1 + 2 * min((phi - 1), 0)) / (phi - 1)
}

#' Fraga Alves, Gomes, de Haan, Neves (2009) Mixed Moment Estimator with PORT Modification
#'
#' Computes the mixed moment tail index estimator with PORT modification from Fraga Alves, Gomes, de Haan, and Neves (2009), as described on page 3.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param p An optional numeric parameter controlling the PORT modification threshold. Defaults to \code{0.9}.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator refines the tail index using mixed moments, a bias adjustment, and the PORT modification:
#' - The PORT modification shifts the data by subtracting \( X_{(k_1+1)} \), where \( k_1 = \text{round}(p \cdot n) \), ensuring better tail behavior.
#' - The first moment \( M_1 \) is computed using the \code{MM} function with \( \alpha = 1 \).
#' - A bias correction term \( L_1 \) is calculated as:
#'   \[
#'   L_1 = \text{mean}\left( 1 - \frac{X_{(n-k)}}{X_{(n-i+1)}} \right), \, i = 1, \ldots, k.
#'   \]
#' - The auxiliary parameter \( \phi \) is then calculated as:
#'   \[
#'   \phi = \frac{M_1 - L_1}{L_1^2}.
#'   \]
#' - The tail index is computed as:
#'   \[
#'   \text{Tail Index} = \frac{1 + 2 \cdot \min(\phi - 1, 0)}{\phi - 1}.
#'   \]
#'
#' This method incorporates the PORT modification to improve robustness and reduce bias in the tail index estimation.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' AGHN_port(x)
#'
#' @export
AGHN_port <- function(X, k, p) {         
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(p)) p <- 0.9
  
  k1 <- round(p * n)
  X <- X - X[k1 + 1]
  
  M1 <- MM(X, k, 1)
  ind <- c(1:k)
  L1 <- mean(1 - X[n - k] / X[n - ind + 1])
  phi <- (M1 - L1) / L1^2
  
  (1 + 2 * min((phi - 1), 0)) / (phi - 1)
}

#' Brilhante et al. (2013) Tail Index Estimator
#'
#' Computes the tail index estimator from Brilhante et al. (2013), incorporating auxiliary bias adjustments and optimal parameter selection.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. If not provided, it is calculated using the method from Gomes et al. (2016).
#' @param p An optional numeric parameter controlling the power transformation. If not provided, it is estimated using auxiliary parameters.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator adjusts the tail index using auxiliary parameters and optimal threshold selection:
#' - The auxiliary parameters \(\rho\) and \(\beta\) are computed using the \code{rho} and \code{beta} functions.
#' - If \( k \) is not provided, it is calculated as:
#'   \[
#'   k = \left( \frac{(1 - \rho) \cdot n^{-\rho}}{\beta \cdot \sqrt{-2 \rho}} \right)^{\frac{2}{1 - 2 \rho}},
#'   \]
#'   with constraints \( 0.1 \cdot \sqrt{n} \leq k \leq 0.9 \cdot n \).
#' - If \( p \) is not provided, it is estimated as:
#'   \[
#'   p = \frac{\phi}{H},
#'   \]
#'   where \( \phi = \max\left(0, \min\left(\frac{\sqrt{2}}{2}, 1 - \frac{\rho}{2} - 0.5 \sqrt{\rho^2 - 4\rho + 2}\right)\right) \) and \( H \) is the inverse Hill estimator.
#' - The tail index is computed as:
#'   \[
#'   \alpha = \frac{1}{\gamma}, \quad \text{where} \quad \gamma = \frac{1 - \frac{1}{\text{mean}(U)}}{p},
#'   \]
#'   and \( U = \left( \frac{X_{(n-k+i)}}{X_{(n-k)}} \right)^p \) for \( i = 1, \ldots, k \).
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' BGP(x)
#'
#' @export
BGP <- function(X, k, p) {          
  X <- sort(X)
  n <- length(X)
  
  r <- rho(X)
  B <- beta(X)
  if (missing(k)) {
    k <- round(((1 - r) * n^(-r) / (B * sqrt(-2 * r)))^(2 / (1 - 2 * r)))  # k - Gomes et al. 2016
    if (k < (0.1 * sqrt(n))) k <- round(0.1 * sqrt(n) + 1)
    if (k > (0.9 * n)) k <- round(0.9 * n)
  } 
  if (missing(p)) {
    H <- 1 / Hill(X)
    r <- rho(X)
    phi <- 1 - r / 2 - 0.5 * sqrt(r^2 - 4 * r + 2)
    phi <- min(max(0, phi), sqrt(2) / 2)
    p <- phi / H
  }
  if (p == 0) {
    alpha <- Hill(X, k)
  } else {
    ind <- 1:k
    U <- (X[n - ind + 1] / X[n - k])^p
    gamma <- (1 - 1 / mean(U)) / p
    alpha <- 1 / gamma
  }
  alpha
}

#' Caeiro et al. (2016) Tail Index Estimator
#'
#' Computes the tail index estimator from Caeiro et al. (2016), incorporating auxiliary parameters and bias reduction techniques.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param p An optional numeric parameter controlling the power transformation. If not provided, it is calculated using auxiliary parameters.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator refines the tail index using auxiliary parameters (\(\rho\), \(\beta\), \(\phi\)) and a bias-adjusted formula:
#' - The auxiliary parameters are calculated as follows:
#'   - \(\rho\) is computed using the \code{rho} function.
#'   - \(\beta\) is computed using the \code{beta} function.
#'   - \(\phi\) is calculated based on the power transformation parameter \( p \) and the Hill estimator \( H \):
#'     - If \( p \) is not provided, \(\phi\) is computed as:
#'       \[
#'       \phi = 1 - \frac{\rho}{2} - 0.5 \sqrt{\rho^2 - 4\rho + 2}.
#'       \]
#'       \( \phi \) is constrained to \( [0, \sqrt{2}/2] \), and \( p \) is determined as:
#'       \[
#'       p = \phi / H.
#'       \]
#'     - If \( p \) is provided, \(\phi = p \cdot H\).
#' - The tail index is computed as:
#'   \[
#'   \text{Tail Index} = \left[ \gamma \cdot \left( 1 - \frac{\beta \cdot (1 - \phi)}{(1 - \rho - \phi)} \cdot \left(\frac{n}{k}\right)^{\rho} \right) \right]^{-1},
#'   \]
#'   where \( \gamma \) is calculated using power-transformed order statistics:
#'   \[
#'   \gamma = \frac{1 - \frac{1}{\text{mean}(U)}}{p},
#'   \]
#'   and \( U = \left( \frac{X_{(n-k+i)}}{X_{(n-k)}} \right)^p \) for \( i = 1, \ldots, k \).
#'
#' This method improves bias reduction by combining power transformations, auxiliary parameters, and the Hill estimator.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' CGBW(x)
#'
#' @export
CGBW <- function(X, k, p) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  
  r <- rho(X)
  B <- beta(X)
  H <- 1 / Hill(X)
  
  if (missing(p)) {
    phi <- 1 - r / 2 - 0.5 * sqrt(r^2 - 4 * r + 2)
    phi <- min(max(0, phi), sqrt(2) / 2)
    p <- phi / H
  } else {
    phi <- p * H
  }
  
  if (p == 0) {
    gamma <- 1 / Hill(X, k)
  } else {
    ind <- 1:k
    U <- (X[n - ind + 1] / X[n - k])^p
    gamma <- (1 - 1 / mean(U)) / p
  }
  (gamma * (1 - B * (1 - phi) / (1 - r - phi) * (n / k)^r))^(-1)
}



#' Gomes and Henriques-Rodrigues (2016) Tail Index Estimator with PORT Modification
#'
#' Computes the tail index estimator from Gomes and Henriques-Rodrigues (2016, REVSTAT), incorporating a PORT modification.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param q An optional numeric parameter controlling the PORT modification threshold. Defaults to \code{0.001}.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator refines the tail index by applying the PORT modification to the data before applying the \code{BGP} estimator:
#' - The PORT modification shifts the data by removing the first \( kk = \text{round}(n \cdot q + 1) \) smallest elements and subtracting \( X_{(kk)} \):
#'   \[
#'   Y = X_{(kk+1):n} - X_{(kk)}.
#'   \]
#' - The \code{BGP} estimator is then applied to the modified data \( Y \) using \( k \) as the number of top-order statistics.
#'
#' This method improves robustness and reduces bias in the tail index estimation by modifying the data to better fit the tail distribution.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GRM_port(x)
#'
#' @export
GRM_port <- function(X, k, q) {     
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(q)) q <- 0.001
  kk <- round(n * q + 1)
  Y <- X[(kk + 1):n] - X[kk]
  BGP(Y, k)
}

#' Gomes and Henriques-Rodrigues (2016) Tail Index Estimator with PORT Modification
#'
#' Computes the tail index estimator from Gomes and Henriques-Rodrigues (2016), using a PORT modification of the Caeiro et al. (2005) estimator.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param q An optional numeric parameter controlling the PORT modification threshold. Defaults to \code{0.001}.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator modifies the Caeiro et al. (2005) \code{CGP1} estimator by incorporating a PORT modification:
#' - The PORT modification shifts the data by removing the first \( kk = \text{round}(n \cdot q + 1) \) smallest elements and subtracting \( X_{(kk)} \):
#'   \[
#'   Y = X_{(kk+1):n} - X_{(kk)}.
#'   \]
#' - The \code{CGP1} estimator is then applied to the modified data \( Y \) using \( k \) as the number of top-order statistics.
#'
#' This approach improves robustness and reduces bias by adapting the \code{CGP1} estimator to better handle the tail of the distribution.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GHR_port(x)
#'
#' @export
GHR_port <- function(X, k, q) {            
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(q)) q <- 0.001
  kk <- round(n * q + 1)
  Y <- X[(kk + 1):n] - X[kk]
  CGP1(Y, k)
}

#' Alves (1995) Tail Index Estimator
#'
#' Computes the tail index estimator from Alves (1995), suitable for positive tail indexes.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param c An optional numeric constant to control the scaling parameter. Defaults to \code{5}.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator computes the tail index using the log ratio of order statistics:
#' - The constant \( c \) scales the order statistics used in the calculation. If \( c \cdot k > n \), it is adjusted to \( c = 2 \).
#' - The tail index is calculated as:
#'   \[
#'   \text{Tail Index} = \frac{\log(c)}{\log\left(\frac{X_{(n-k+1)}}{X_{(n-c \cdot k + 1)}}\right)}.
#'   \]
#' This method is designed for positive tail indexes and adjusts the scaling dynamically if the scaling factor exceeds the sample size.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' Alves1995(x)
#'
#' @export
Alves1995 <- function(X, k, c) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(c)) c <- 5
  if ((c * k) > n) c <- 2
  log(c) / log(X[n - k + 1] / X[n - c * k + 1])
}

#' Alves (1995) Tail Index Estimator
#'
#' Computes the tail index estimator from Alves (1995), suitable for both positive and negative tail indexes.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @param k0 An optional integer specifying the number of order statistics used for the auxiliary calculation. Defaults to \( k^{3/4} \).
#' @return A numeric value representing the estimated tail index.
#' @details This estimator allows for the calculation of both positive and negative tail indexes by analyzing ratios of differences between order statistics:
#' - The auxiliary parameter \( k_0 \) controls the subset of top-order statistics for the ratio calculation and defaults to \( k^{3/4} \).
#' - The tail index is calculated as:
#'   \[
#'   \text{Tail Index} = \left[ \frac{1}{k_0} \sum_{i=1}^{k_0} \log\left(\frac{X_{(n-k_0+i)} - X_{(n-k)}}{X_{(n-k_0)} - X_{(n-k)}}\right) \right]^{-1}.
#'   \]
#'
#' This method provides a flexible approach for estimating tail indexes, including those with negative values.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' Alves(x)
#'
#' @export
Alves <- function(X, k, k0) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(k0)) k0 <- round(k^(3/4))
  TT <- (X[(n - k0 + 1):n] - X[n - k]) / (X[n - k0] - X[n - k])
  (sum(log(TT)) / k0)^(-1)
}

#' Weis Tail Index Estimator
#'
#' Computes the tail index estimator using the Weis method, based on differences in log-transformed order statistics.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to \code{2 \cdot \log(n)}, where \( n \) is the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator calculates the tail index by analyzing the ratio of differences between selected order statistics:
#' - The default number of top-order statistics \( k \) is calculated as \( 2 \cdot \log(n) \).
#' - The tail index is computed as:
#'   \[
#'   \text{Tail Index} = \frac{\log(2)}{\log\left(\frac{X_{(k)} - X_{(1)}}{X_{(k/2)} - X_{(1)}}\right)}.
#'   \]
#' This method assumes that \( k \) is large enough to capture the tail behavior but small enough to avoid noise from the bulk of the distribution.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' Weis(x)
#'
#' @export
Weis <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(log(n)) * 2
  
  D1 <- log((X[k] - X[1]) / (X[(k / 2)] - X[1]))
  log(2) / D1
}


FH_MLE1 <- function(par, X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  D <- par[1]
  beta <- par[2]
  U <- (1:(n - 1)) * (log(X[n:2]) - log(X[(n - 1):1]))
  ind <- c(1:(n - 1)) / n
  W <- U * exp(-D * ind^(-beta))
  M1 <- mean(ind[1:k]^(-beta))
  M2 <- mean(W[1:k])
  D * M1 + log(M2)
}

#' Feuerverger and Hall (1999) Maximum Likelihood Tail Index Estimator
#'
#' Computes the tail index estimator using the Maximum Likelihood Estimator (MLE) approach proposed by Feuerverger and Hall (1999).
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This function optimizes the MLE parameters (\code{D} and \(\rho\)) using the \code{FH_MLE1} function and applies the estimated parameters to compute the tail index:
#' - The parameters \code{D} and \(\rho\) are estimated by minimizing the negative log-likelihood via \code{optim}.
#'   - The optimization starts with initial values of \code{D = 1} and \(\rho = 1\).
#'   - The optimization is performed using the "L-BFGS-B" method with constraints: \(\rho \geq 0\).
#' - Using the optimized parameters, the tail index is calculated as:
#'   \[
#'   \text{Tail Index} = \left[ \frac{1}{k} \sum_{i=1}^{k} U_i \cdot \exp\left(-D \cdot \text{ind}_i\right)^{-\rho} \right]^{-1},
#'   \]
#'   where:
#'   \[
#'   U = i \cdot \left( \log(X_{(n-i+1)}) - \log(X_{(n-i)}) \right), \quad \text{ind}_i = \frac{i}{n}.
#'   \]
#'
#' This method combines bias reduction through parameterized weighting and maximum likelihood estimation.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' FH_MLE(x)
#'
#' @export
FH_MLE <- function(X, k) { 
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  L <- optim(c(1, 1), FH_MLE1, X = X, k = k, method = "L-BFGS-B", lower = c(-Inf, 0), upper = c(Inf, Inf))
  D <- L$par[1]
  rho <- L$par[2]
  U <- (1:(n - 1)) * (log(X[n:2]) - log(X[(n - 1):1]))
  ind <- c(1:(n - 1)) / n
  mean(U[1:k] * exp(-D * ind[1:k])^(-rho))^(-1)
}

#' Feuerverger and Hall (1999) OLS Tail Index Estimator
#'
#' Computes the tail index estimator using the Ordinary Least Squares (OLS) approach proposed by Feuerverger and Hall (1999).
#'
#' @param par A numeric vector of length 3 containing the parameters \code{mu}, \code{D}, and \code{beta} to be optimized.
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to \(2 \cdot \log(n)\), where \( n \) is the sample size.
#' @return A numeric value representing the sum of squared residuals for the given parameters \code{mu}, \code{D}, and \code{beta}.
#' @details This function implements the OLS-based approach to tail index estimation:
#' - The parameters \code{mu}, \code{D}, and \code{beta} are optimized by minimizing the sum of squared residuals.
#' - Intermediate quantities are computed as follows:
#'   - \( U \): The scaled differences between order statistics:
#'     \[
#'     U = i \cdot \left(X_{(n-i+1)} - X_{(n-i)}\right), \, i = 1, \ldots, n-1.
#'     \]
#'   - \( V \): The logarithm of \( U \):
#'     \[
#'     V = \log(U).
#'     \]
#' - The sum of squared residuals is calculated as:
#'   \[
#'   \text{Residuals} = \sum_{i=1}^{k} \left(V_i - \mu - D \cdot \left(\frac{i}{n}\right)^{-\beta}\right)^2.
#'   \]
#'
#' This function is designed for use in an optimization routine to find the best-fitting values for \code{mu}, \code{D}, and \code{beta}.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' params <- c(0, 1, 0.5)  # Initial parameter guesses
#' FH_MLE2(params, x)
#'
#' @export
FH_MLE2 <- function(par, X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(log(n)) * 2
  mu <- par[1]
  D <- par[2]
  beta <- par[3]
  U <- (1:(n - 1)) * (X[n:2] - X[(n - 1):1])
  ind <- c(1:(n - 1)) / n
  V <- log(U)
  sum((V[1:k] - mu - D * ind[1:k]^(-beta))^2)
}

#' Gomes and Martin (2004) Variant of Feuerverger and Hall MLE Estimator
#'
#' Computes the tail index estimator using the Gomes and Martin (2004) variant of the Feuerverger and Hall Maximum Likelihood Estimator (MLE).
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator refines the Feuerverger and Hall MLE by incorporating additional moment calculations:
#' - Intermediate quantities are computed as follows:
#'   - \( U \): The scaled log differences:
#'     \[
#'     U = i \cdot \left(\log(X_{(n-i+1)}) - \log(X_{(n-i)})\right), \, i = 1, \ldots, n-1.
#'     \]
#'   - \( g_1 \): The Hill estimator for the tail index:
#'     \[
#'     g_1 = \frac{1}{\text{Hill}(X, k)}.
#'     \]
#'   - \( g_2 \): The mean of \( U \) over the top \( k \) order statistics:
#'     \[
#'     g_2 = \text{mean}(U_{1:k}).
#'     \]
#'   - \( g_3 \): The weighted sum of \( U \):
#'     \[
#'     g_3 = \sum_{i=1}^k \left(2i - k - 1\right) \cdot U_i.
#'     \]
#'   - \( g_4 \): The doubly weighted sum of \( U \):
#'     \[
#'     g_4 = \sum_{i=1}^k i \cdot \left(2i - k - 1\right) \cdot U_i.
#'     \]
#' - The tail index is calculated as:
#'   \[
#'   \text{Tail Index} = \left(g_1 - g_2 \cdot \frac{g_3}{g_4}\right)^{-1}.
#'   \]
#'
#' This method enhances the original MLE by incorporating additional weighting terms to reduce bias.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GM_FH(x)
#'
#' @export
GM_FH <- function(X, k) {      
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  U <- (1:(n - 1)) * (log(X[n:2]) - log(X[(n - 1):1]))
  g1 <- Hill(X, k)^(-1)
  g2 <- mean(U[1:k])
  g3 <- sum((2 * c(1:k) - k - 1) * U[1:k])
  g4 <- sum(c(1:k) * (2 * c(1:k) - k - 1) * U[1:k])
  (g1 - g2 * g3 / g4)^(-1)
}



#' Gomes and Martin (2004) Variant of Feuerverger and Hall OLS Estimator
#'
#' Computes the tail index estimator using the Ordinary Least Squares (OLS) approach proposed by Gomes and Martin (2004), based on the Feuerverger and Hall methodology.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator refines the Feuerverger and Hall OLS method by incorporating weighting and logarithmic transformations:
#' - Intermediate quantities are computed as follows:
#'   - \( U \): The scaled log differences:
#'     \[
#'     U = i \cdot \left(\log(X_{(n-i+1)}) - \log(X_{(n-i)})\right), \, i = 1, \ldots, n-1.
#'     \]
#'   - \( g_1 \): A weighted sum of log-transformed \( U \):
#'     \[
#'     g_1 = \frac{2(2k+1)}{k(k-1)} \cdot \sum_{i=1}^k \log(U_i).
#'     \]
#'   - \( g_2 \): Another weighted sum of log-transformed \( U \):
#'     \[
#'     g_2 = \frac{6}{k(k-1)} \cdot \sum_{i=1}^k i \cdot \log(U_i).
#'     \]
#' - The tail index is calculated as:
#'   \[
#'   \text{Tail Index} = \left[\exp\left(g_1 - g_2 + \gamma\right)\right]^{-1},
#'   \]
#'   where \(\gamma \approx 0.5772157\) is the Euler-Mascheroni constant.
#'
#' This method enhances the original OLS approach by introducing additional weights to reduce bias and improve estimation accuracy.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' GM_FH2(x)
#'
#' @export
GM_FH2 <- function(X, k) {      
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  U <- (1:(n - 1)) * (log(X[n:2]) - log(X[(n - 1):1]))
  g1 <- 2 * (2 * k + 1) / (k * (k - 1)) * sum(log(U[1:k]))
  g2 <- 6 / (k * (k - 1)) * sum(c(1:k) * log(U[1:k]))
  exp(g1 - g2 + digamma(1))^(-1)
}

#' SAG Tail Index Estimator
#'
#' Computes the tail index using the SAG estimator, which incorporates a threshold adjustment and applies the Hill estimator on the adjusted data.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size of the adjusted data.
#' @param q An optional numeric parameter controlling the threshold for the adjustment. Defaults to \code{0.1}.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator applies a threshold adjustment to the data before computing the tail index using the Hill estimator:
#' - The adjustment threshold is defined as \( m = \text{round}(n \cdot q + 1) \), where \( q \) controls the fraction of the sample used for adjustment.
#' - The adjusted data is calculated as:
#'   \[
#'   Y = X_{(m+1):n} - X_{(m)}.
#'   \]
#' - The Hill estimator is then applied to the adjusted data \( Y \), using \( k \) as the number of top-order statistics.
#'
#' This approach ensures that the tail index estimation focuses on the extreme values of the distribution, while removing potential bias introduced by smaller values.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' SAG(x)
#'
#' @export
SAG <- function(X, k, q) {          
  X <- sort(X)
  if (missing(q)) q <- 0.1
  m <- round(length(X) * q + 1)
  Y <- X[(m + 1):length(X)] - X[m]
  
  n <- length(Y)
  if (missing(k)) k <- round(sqrt(n))
  
  Hill(Y, k)
}

#' Smith Tail Index Estimator
#'
#' Computes the tail index using the Smith estimator, which applies a logarithmic transformation to exceedances above a threshold.
#'
#' @param X A numeric vector containing the data sample.
#' @param u An optional numeric value specifying the threshold. If not provided, the threshold is set as the \( (n - k) \)-th largest value, where \( k = \sqrt{n} \).
#' @return A numeric value representing the estimated tail index.
#' @details This estimator focuses on exceedances above a specified threshold to calculate the tail index:
#' - If the threshold \( u \) is not provided, it is set to:
#'   \[
#'   u = X_{(n-k)},
#'   \]
#'   where \( k = \sqrt{n} \).
#' - The exceedances \( Y \) above the threshold \( u \) are defined as:
#'   \[
#'   Y = X[X > u] - u.
#'   \]
#' - The transformed values \( Z \) are calculated as:
#'   \[
#'   Z = \log\left(\frac{Y}{u} + 1\right).
#'   \]
#' - The tail index is estimated as the inverse of the mean of \( Z \):
#'   \[
#'   \text{Tail Index} = \left(\text{mean}(Z)\right)^{-1}.
#'   \]
#'
#' This method effectively captures the behavior of the tail by focusing on values exceeding a dynamically chosen threshold.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' Smith(x)
#'
#' @export
Smith <- function(X, u) {          
  X <- sort(X)
  n <- length(X)
  if (missing(u)) {
    k <- round(sqrt(n))
    u <- X[n - k]
  }
  Y <- (X[X > u] - u)
  Z <- log(Y / u + 1)
  mean(Z)^(-1)
}

#' Baek and Pipiras (2010) Tail Index Estimator
#'
#' Computes the tail index using the Baek and Pipiras (2010) estimator, which is based on a regression model incorporating logarithmic and inverse transformations of order statistics.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator applies a linear regression model to transformed order statistics to estimate the tail index:
#' - The top \( k \) order statistics are selected for the calculation, where \( k \) defaults to \( \sqrt{n} \).
#' - The regression model is constructed using:
#'   - \( y \): The logarithm of scaled ranks:
#'     \[
#'     y = \log\left(\frac{i - 0.5}{k}\right), \, i = 1, \ldots, k.
#'     \]
#'   - \( X_1 \): The negative logarithm of scaled order statistics:
#'     \[
#'     X_1 = -\log\left(\frac{X_{(n-i+1)}}{X_{(n-k)}}\right).
#'     \]
#'   - \( X_2 \): The inverse of scaled order statistics:
#'     \[
#'     X_2 = \frac{X_{(n-k)}}{X_{(n-i+1)}}, \, i = 1, \ldots, k.
#'     \]
#' - A linear regression model is fitted:
#'   \[
#'   y = \beta_1 \cdot X_1 + \beta_2 \cdot X_2 + \varepsilon,
#'   \]
#'   and the tail index is estimated as \( \beta_1 \).
#'
#' This method incorporates both logarithmic and inverse transformations to improve the accuracy of the tail index estimation.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' BP(x)
#'
#' @export
BP <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  
  ind <- c(1:k)
  y <- log((ind - 0.5) / k)
  
  X1 <- -log(X[n - ind + 1] / X[n - k])
  X2 <- X[n - k] / X[n - ind + 1]
  f <- lm(y ~ X1 + X2)
  f$coefficients[2]
}

#' Baek and Pipiras (2010) Tail Index Estimator (Equation 2.16)
#'
#' Computes the tail index using the Baek and Pipiras (2010) estimator based on a regression model from their equation 2.16.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator applies a linear regression model to transformed order statistics to estimate the tail index:
#' - The top \( k \) order statistics are selected for the calculation, where \( k \) defaults to \( \sqrt{n} \).
#' - The regression model is constructed using:
#'   - \( X_1 \): The logarithm of scaled ranks:
#'     \[
#'     X_1 = \log\left(\frac{i - 0.5}{k}\right), \, i = 1, \ldots, k.
#'     \]
#'   - \( y \): The negative logarithm of scaled order statistics:
#'     \[
#'     y = -\log\left(\frac{X_{(n-i+1)}}{X_{(n-k)}}\right).
#'     \]
#'   - \( X_2 \): The inverse of scaled order statistics:
#'     \[
#'     X_2 = \frac{X_{(n-k)}}{X_{(n-i+1)}}, \, i = 1, \ldots, k.
#'     \]
#' - A linear regression model is fitted:
#'   \[
#'   y = \beta_1 \cdot X_1 + \beta_2 \cdot X_2 + \varepsilon,
#'   \]
#'   and the tail index is estimated as the inverse of \( \beta_2 \):
#'   \[
#'   \text{Tail Index} = (\beta_2)^{-1}.
#'   \]
#'
#' This method extends the original Baek and Pipiras estimator by incorporating additional transformations and regression coefficients for tail index estimation.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' BP1(x)
#'
#' @export
BP1 <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  
  ind <- c(1:k)
  X1 <- log((ind - 0.5) / k)
  
  y <- -log(X[n - ind + 1] / X[n - k])
  X2 <- X[n - k] / X[n - ind + 1]
  f <- lm(y ~ X1 + X2)
  (f$coefficients[2])^(-1)
}

#' Baek and Pipiras (2010) Tail Index Estimator (Equation 2.17)
#'
#' Computes the tail index using the Baek and Pipiras (2010) estimator based on a regression model from their equation 2.17, with an adjustment for bias correction.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator applies a linear regression model and incorporates bias correction to estimate the tail index:
#' - The top \( k \) order statistics are selected for the calculation, where \( k \) defaults to \( \sqrt{n} \).
#' - The Hill estimator is first computed:
#'   \[
#'   H = \frac{1}{\text{Hill}(X, k)}.
#'   \]
#' - The regression model is constructed using:
#'   - \( X_1 \): The logarithm of scaled ranks:
#'     \[
#'     X_1 = \log\left(\frac{i - 0.5}{k}\right), \, i = 1, \ldots, k.
#'     \]
#'   - \( y \): The negative logarithm of scaled order statistics:
#'     \[
#'     y = -\log\left(\frac{X_{(n-i+1)}}{X_{(n-k)}}\right).
#'     \]
#'   - \( X_2 \): The inverse of scaled order statistics:
#'     \[
#'     X_2 = \frac{X_{(n-k)}}{X_{(n-i+1)}}, \, i = 1, \ldots, k.
#'     \]
#' - A linear regression model is fitted:
#'   \[
#'   y = \beta_1 \cdot X_1 + \beta_2 \cdot X_2 + \varepsilon,
#'   \]
#'   where \( \beta_2 \) is used for bias correction.
#' - The bias correction term is calculated as:
#'   \[
#'   \text{Bias Correction} = \beta_2 \cdot \frac{1}{k} \sum_{i=1}^k \left(X_{(n-i+1)}^{-1} - X_{(n-k+1)}^{-1}\right).
#'   \]
#' - The tail index is estimated as:
#'   \[
#'   \text{Tail Index} = \left(H - \text{Bias Correction}\right)^{-1}.
#'   \]
#'
#' This method enhances the original Baek and Pipiras estimator by including a bias correction term to improve accuracy.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' BP2(x)
#'
#' @export
BP2 <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  
  H <- Hill(X, k)^(-1)
  
  ind <- c(1:k)
  X1 <- log((ind - 0.5) / k)
  
  y <- -log(X[n - ind + 1] / X[n - k])
  X2 <- X[n - k] / X[n - ind + 1]
  f <- lm(y ~ X1 + X2)
  beta <- f$coefficients[3]
  
  XX <- beta * sum(X[n - ind + 1]^(-1) - X[n - k + 1]^(-1)) / k
  (H - XX)^(-1)
}

#' Baek and Pipiras (2010) Tail Index Estimator (Equation 2.18)
#'
#' Computes the tail index using the Baek and Pipiras (2010) estimator based on their equation 2.18, which incorporates the Hill estimator and second moments.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator combines the Hill estimator and the second moment of logarithmic differences to compute the tail index:
#' - The Hill estimator is calculated as:
#'   \[
#'   H = \text{Hill}(X, k).
#'   \]
#' - The second moment of logarithmic differences is computed as:
#'   \[
#'   M_2 = \text{MM}(X, k, 2),
#'   \]
#'   where \(\text{MM}(X, k, \alpha)\) computes the mean of the \(\alpha\)-th power of the logarithmic differences.
#' - The tail index is then calculated as:
#'   \[
#'   \text{Tail Index} = \left[\frac{(2 + 1/H) \cdot 2}{M_2} - 2 \cdot H \cdot \sqrt{\frac{2}{M_2}}\right].
#'   \]
#'
#' This method leverages higher-order moments and the Hill estimator to enhance the accuracy of tail index estimation.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' BP3(x)
#'
#' @export
BP3 <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  
  H <- Hill(X, k)
  M2 <- MM(X, k, 2)
  
  (2 + 1 / H) * 2 / M2 - 2 * H * (2 / M2)^0.5
}

#' Baek and Pipiras (2010) Tail Index Estimator (Equation 2.19)
#'
#' Computes the tail index using the Baek and Pipiras (2010) estimator based on their equation 2.19, which combines the Hill estimator and second moments.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the square root of the sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator refines the tail index by incorporating both the Hill estimator and the second moment of logarithmic differences:
#' - The Hill estimator is calculated as:
#'   \[
#'   H = \text{Hill}(X, k).
#'   \]
#' - The second moment of logarithmic differences is computed as:
#'   \[
#'   M_2 = \text{MM}(X, k, 2),
#'   \]
#'   where \(\text{MM}(X, k, \alpha)\) computes the mean of the \(\alpha\)-th power of the logarithmic differences.
#' - The tail index is then calculated as:
#'   \[
#'   \text{Tail Index} = -H^2 + \left(\frac{1}{H} + 1\right) \cdot \frac{2}{M_2}.
#'   \]
#'
#' This method combines the Hill estimator with higher-order moments to improve tail index estimation accuracy.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' BP4(x)
#'
#' @export
BP4 <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  
  H <- Hill(X, k)
  M2 <- MM(X, k, 2)
  
  -H^2 + (1 / H + 1) * 2 / M2
}

#' Hosking-Wallis Tail Index Estimator (Method of Moments)
#'
#' Computes the tail index using the Hosking-Wallis method of moments estimator, which is based on the mean and variance of the top-order statistics.
#'
#' @param X A numeric vector containing the data sample.
#' @param k An optional integer specifying the number of top-order statistics to use. Defaults to the full sample size.
#' @return A numeric value representing the estimated tail index.
#' @details This estimator uses the mean and variance of the top-order statistics to compute the tail index:
#' - The top \( k \) order statistics are selected, where \( k \) defaults to the full sample size (\( n \)).
#' - The mean (\( \bar{x} \)) and variance (\( s^2 \)) of the selected order statistics are calculated:
#'   \[
#'   \bar{x} = \text{mean}(X_{(n-k+1):n}),
#'   \quad s^2 = \text{var}(X_{(n-k+1):n}).
#'   \]
#' - The tail index is then computed as:
#'   \[
#'   \text{Tail Index} = \left[0.5 \cdot \left(1 - \frac{\bar{x}^2}{s^2}\right)\right]^{-1}.
#'   \]
#'
#' This method leverages the relationship between the mean and variance of extreme values to estimate the tail index.
#'
#' @examples
#' set.seed(123)
#' n = 1000
#' x = (abs(stats::rcauchy(n)))^(2)  # tail index = 0.5.
#' HWm(x)
#'
#' @export
HWm <- function(X, k) {          
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- n #k <- round(sqrt(n))
  
  x <- mean(X[(n - k + 1):n])
  s2 <- var(X[(n - k + 1):n])
  
  (0.5 * (1 - x^2 / s2))^(-1)
}

