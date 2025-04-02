my_lasso <- function(y, X, lambda) {
  n <- length(y)
  p <- ncol(X)

  X_prime <- cbind(X, -X)

  Vmat <- t(X_prime) %*% X_prime
  dvec <- -t(X_prime) %*% y

  # We want sum of beta+ and beta- to be less than lambda
  Amat <- t(c(rep(1, p), rep(1, p)))
  bvec <- c(lambda)
  uvec <- c(rep(100, p), rep(100, p))


  fit <- LowRankQP(Vmat = Vmat,
                   dvec = dvec,
                   Amat = Amat,
                   bvec = bvec,
                   uvec = uvec)

  # Extract beta = beta+ - beta-
  beta_plus  <- fit$alpha[1:p]
  beta_minus <- fit$alpha[(p + 1):(2 * p)]

  beta <- beta_plus - beta_minus

  return(beta)

}
