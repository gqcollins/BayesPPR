\dontrun{
  ################################
  ### univariate example
  ################################
  ## simulate data (Friedman function)
  f <- function(X){
    10*sin(pi*X[,1]*X[,2]) + 20*(X[,3] - .5)^2 + 10*X[,4] + 5*X[,5]
  }
  n <- 500 # number of observations
  p <- 10 #10 variables, only first 5 matter
  X <- matrix(runif(n*p), n, p)
  y <- f(X) + rnorm(n)

  ## fit BPPR
  fit <- bppr(X, y)

  ## prediction
  X_test <- matrix(runif(n*p), n, p)
  preds <- predict(fit, X_test, n_ignore = 9000) # posterior predictive samples
  true_f <- f(X_test)
  plot(true_f, colMeans(preds), xlab = 'true values', ylab='posterior predictive means')
  abline(a=0, b=1, col=2)
}

## minimal example for CRAN testing
fit <- bppr(matrix(1:10), 1:10, n_draws = 2)
