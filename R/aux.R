# Auxiliary functions:

MM <- function(X, k, alpha = 2) {
  n <- length(X)
  sum((log(X[(n - k + 1):n]) - log(X[n - k]))^alpha) / k
}

rho <- function(X, k, tau) {
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  if (missing(tau)) tau <- 0
  M1 <- MM(X, k, 1)
  M2 <- MM(X, k, 2)
  M3 <- MM(X, k, 3)
  if (tau == 0) {
    T <- (log(M1) - log(M2 / 2) / 2) / (log(M2 / 2) / 2 - log(M3 / 6) / 3)
  } else {
    T <- (M1^tau - (M2 / 2)^(tau / 2)) / ((M2 / 2)^(tau / 2) - (M3 / 6)^(tau / 3))
  }
  min(-abs(3 * (T - 1) / (T - 3)), 0)
}

beta <- function(X, k) {
  X <- sort(X)
  n <- length(X)
  if (missing(k)) k <- round(sqrt(n))
  I <- c(1:k) / k
  U <- c((n - 1):1) * (log(X[2:n]) - log(X[1:(n - 1)]))
  rho1 <- rho(X, k, 0)
  M0 <- (k / n)^rho1
  M1 <- mean(I^rho1)
  M2 <- mean(U[1:k])
  M3 <- mean(U[1:k] * I^rho1)
  M4 <- mean(U[1:k] * I^(2 * rho1))
  M0 * (M1 - M2) / (M1 * M2 - M3)
}