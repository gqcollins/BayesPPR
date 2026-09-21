# Tests for the multivariate-response path: bppr_pca(), predict.bppr_pca(),
# and the internal pca_setup()/pca_reverse() helpers.

make_mv_data <- function(n = 60, p = 3, D = 5, seed = 13) {
  set.seed(seed)
  X <- matrix(runif(n * p), n, p)
  # Response dimensions are smooth functions of X so PCA can compress them.
  Y <- sapply(1:D, function(k) X[, 1] * k + sin(pi * X[, 2]) + rnorm(n, sd = 0.05))
  list(X = X, Y = Y)
}

test_that("bppr_pca() returns a well-formed bppr_pca object", {
  d <- make_mv_data()
  fit <- bppr_pca(d$X, d$Y, n_pc = 2, n_post = 15, n_burn = 5, n_adapt = 0,
                  print_every = 0)

  expect_s3_class(fit, "bppr_pca")
  expect_equal(fit$pca_Y$n_pc, 2)
  expect_length(fit$fit_list, 2)
  expect_s3_class(fit$fit_list[[1]], "bppr")
})

test_that("bppr_pca() selects n_pc automatically from prop_var", {
  d <- make_mv_data()
  fit <- bppr_pca(d$X, d$Y, n_pc = NULL, prop_var = 0.99, n_post = 10,
                  n_burn = 5, n_adapt = 0, print_every = 0)
  expect_true(fit$pca_Y$n_pc >= 1)
  expect_true(fit$pca_Y$n_pc <= min(nrow(d$Y), ncol(d$Y)))
})

test_that("predict.bppr_pca() returns a 3-D array of the right shape", {
  d <- make_mv_data()
  fit <- bppr_pca(d$X, d$Y, n_pc = 2, n_post = 15, n_burn = 5, n_adapt = 0,
                  print_every = 0)
  X_test <- matrix(runif(8 * 3), 8, 3)
  preds <- predict(fit, X_test)

  # Dimensions: draws x observations x response-dimension.
  expect_equal(length(dim(preds)), 3)
  expect_equal(dim(preds), c(fit$fit_list[[1]]$n_keep, nrow(X_test), ncol(d$Y)))
  expect_false(anyNA(preds))
})

test_that("predict.bppr_pca() accepts a single test row", {
  # Regression test: single-row newdata used to fail.
  d <- make_mv_data()
  fit <- bppr_pca(d$X, d$Y, n_pc = 2, n_post = 10, n_burn = 5, n_adapt = 0,
                  print_every = 0)
  preds <- predict(fit, matrix(runif(3), nrow = 1))
  expect_equal(dim(preds), c(fit$fit_list[[1]]$n_keep, 1, ncol(d$Y)))
})

test_that("predict.bppr_pca() recovers training response reasonably", {
  d <- make_mv_data()
  fit <- bppr_pca(d$X, d$Y, n_pc = 3, n_post = 20, n_burn = 10, n_adapt = 0,
                  print_every = 0)
  preds <- predict(fit, d$X)
  yhat <- apply(preds, c(2, 3), mean)      # obs x response-dim posterior mean
  expect_equal(dim(yhat), dim(d$Y))
  expect_true(cor(as.vector(yhat), as.vector(d$Y)) > 0.8)
})

test_that("bppr_pca() and predict validate par_type", {
  d <- make_mv_data()
  expect_error(
    bppr_pca(d$X, d$Y, n_pc = 2, par_type = "bogus", n_post = 5, n_burn = 0,
             n_adapt = 0, print_every = 0),
    "par_type"
  )
})

test_that("bppr_resume() works on a bppr_pca object", {
  d <- make_mv_data()
  fit <- bppr_pca(d$X, d$Y, n_pc = 2, n_post = 15, n_burn = 5, n_adapt = 0,
                  print_every = 0)
  fit2 <- bppr_resume(fit, n_post = 10, n_burn = 5, n_adapt = 0,
                      print_every = 0)
  expect_s3_class(fit2, "bppr_pca")
  expect_equal(fit2$fit_list[[1]]$n_keep, 10)

  fit3 <- bppr_resume(fit, append = TRUE, n_post = 10, n_burn = 5, n_adapt = 0,
                      print_every = 0)
  expect_equal(fit3$fit_list[[1]]$n_keep, fit$fit_list[[1]]$n_keep + 10)
})

# --- pca_setup() / pca_reverse() internals -----------------------------------

test_that("pca_setup() rejects invalid prop_var", {
  d <- make_mv_data()
  expect_error(BayesPPR:::pca_setup(d$X, d$Y, prop_var = -0.1), "prop_var")
  expect_error(BayesPPR:::pca_setup(d$X, d$Y, prop_var = 1.1), "prop_var")
})

test_that("pca_setup() rejects univariate Y", {
  set.seed(1)
  X <- matrix(runif(30), 15, 2)
  y <- matrix(rnorm(15), ncol = 1)
  expect_error(BayesPPR:::pca_setup(X, y), "univariate")
})

test_that("pca_setup() transposes Y when its rows match observations of X", {
  set.seed(2)
  n <- 10; D <- 4
  X <- matrix(runif(n * 2), n, 2)
  Yt <- matrix(rnorm(D * n), D, n)          # D x n -- needs transposing
  ps <- BayesPPR:::pca_setup(X, Yt)
  expect_equal(nrow(ps$Y), n)
  expect_equal(ncol(ps$Y), D)
})

test_that("pca_setup() errors on an X/Y dimension mismatch", {
  set.seed(3)
  X <- matrix(runif(20), 10, 2)
  Y <- matrix(rnorm(21), 7, 3)              # neither dim matches n = 10
  expect_error(BayesPPR:::pca_setup(X, Y), "mismatch")
})

test_that("pca_setup() caps n_pc at min(n, D) and warns", {
  set.seed(4)
  n <- 8; D <- 5
  X <- matrix(runif(n * 2), n, 2)
  Y <- matrix(rnorm(n * D), n, D)
  expect_warning(ps <- BayesPPR:::pca_setup(X, Y, n_pc = 99), "n_pc too large")
  expect_equal(ps$n_pc, min(n, D))
})

test_that("pca_setup() rejects n_pc < 1", {
  d <- make_mv_data()
  expect_error(BayesPPR:::pca_setup(d$X, d$Y, n_pc = 0), "at least 1")
})

test_that("pca_setup() handles prop_var = 1 by keeping all components", {
  set.seed(6)
  n <- 12; D <- 4
  X <- matrix(runif(n * 2), n, 2)
  Y <- matrix(rnorm(n * D), n, D)
  ps <- BayesPPR:::pca_setup(X, Y, n_pc = NULL, prop_var = 1)
  expect_equal(ps$n_pc, min(n, D))
})

test_that("pca_reverse() inverts the PCA standardisation", {
  d <- make_mv_data(D = 4)
  ps <- BayesPPR:::pca_setup(d$X, d$Y, n_pc = min(nrow(d$Y), ncol(d$Y)))
  # Reconstructing from all PCs of the scores recovers the original Y.
  recon <- BayesPPR:::pca_reverse(ps$Y_new, ps)     # D x n
  expect_equal(t(recon), d$Y, tolerance = 1e-8)
})
