# Tests that the RJMCMC sampler is reproducible given a fixed random seed.

test_that("bppr() produces identical chains from the same seed", {
  make <- function() {
    set.seed(4242)
    n <- 60
    X <- matrix(runif(n * 3), n, 3)
    y <- X[, 1] + X[, 2]^2 + rnorm(n, sd = 0.1)
    set.seed(99)
    bppr(X, y, n_post = 25, n_burn = 10, n_adapt = 5, print_every = 0)
  }
  fit1 <- make()
  fit2 <- make()

  expect_equal(fit1$n_ridge, fit2$n_ridge)
  expect_equal(fit1$sd_resid, fit2$sd_resid)
  expect_equal(fit1$coefs, fit2$coefs)
  expect_equal(fit1$var_coefs, fit2$var_coefs)
})

test_that("bppr_pca() produces identical results from the same seed (serial)", {
  make <- function() {
    set.seed(2024)
    n <- 50; p <- 3; D <- 4
    X <- matrix(runif(n * p), n, p)
    Y <- sapply(1:D, function(k) X[, 1] * k + sin(pi * X[, 2]) + rnorm(n, 0.05))
    set.seed(7)
    bppr_pca(X, Y, n_pc = 2, n_cores = 1, n_post = 15, n_burn = 5, n_adapt = 0,
             print_every = 0)
  }
  fit1 <- make()
  fit2 <- make()
  expect_equal(fit1$fit_list[[1]]$coefs, fit2$fit_list[[1]]$coefs)
  expect_equal(fit1$fit_list[[2]]$sd_resid, fit2$fit_list[[2]]$sd_resid)
})
