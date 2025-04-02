my_rq <- function(y, X, tau) {
  n <- length(y)
  p <- ncol(X)

  # See slide 16 of class notes
  # u <- [yi - xi*beta]+
  # v <- [yi - xi*beta]-
  # We need 2*p instead of just p because we're splitting coefficients into
  # positive and negative

  #0*b+         #0*b-    # t * ui      #(1 - t)vi
  obj <- -1 * c(rep(0, p), rep(0, p), rep(tau, n), rep(1 - tau, n))

  # Constraints:
  # y_i = X_i b+ - X_i b- + u_i - v_i
  A <- cbind(X, -X, diag(n), -diag(n))
  b <- y

  # Using A3 for equality here
  result <- simplex(a = obj, A3 = A, b3 = b, maxi = TRUE)

  # Remember b+ is 1:p and b- is p+1:2p and y_i = X_i b+ - X_i b- + u_i - v_i
  solution <- result$soln

  b_plus <- solution[1:p]
  b_minus <- solution[(p + 1):(2 * p)]

  # Get final coefficients: beta = b+ - b- for xi
  beta_hat <- b_plus - b_minus

  beta_df <- as_tibble(as.list(beta_hat)) %>%
    mutate(tau = tau) %>%
    rename(intercept = x1) %>%
    select(tau, everything())


  return (beta_df)
}
