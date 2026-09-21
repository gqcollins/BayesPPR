# Tests for predict.bppr(): output shape, idx_use handling and validation,
# and consistency between fitting and prediction.

fit_small <- function(n = 60, p = 4, seed = 3, ...) {
  set.seed(seed)
  X <- matrix(runif(n * p), n, p)
  y <- 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 + rnorm(n, sd = 0.1)
  fit <- bppr(X, y, n_post = 20, n_burn = 10, n_adapt = 0, print_every = 0, ...)
  list(fit = fit, X = X, y = y)
}

test_that("predict.bppr() returns one row per kept draw and one column per obs", {
  s <- fit_small()
  X_test <- matrix(runif(10 * 4), 10, 4)
  preds <- predict(s$fit, X_test)

  expect_true(is.matrix(preds))
  expect_equal(nrow(preds), s$fit$n_keep)
  expect_equal(ncol(preds), nrow(X_test))
  expect_false(anyNA(preds))
})

test_that("predict.bppr() predicting on training data recovers the response", {
  s <- fit_small()
  preds <- predict(s$fit, s$X)
  yhat <- colMeans(preds)
  # A reasonable fit should track the observed response.
  expect_true(cor(yhat, s$y) > 0.8)
})

test_that("predict.bppr() honours idx_use", {
  s <- fit_small()
  X_test <- matrix(runif(5 * 4), 5, 4)

  preds_all <- predict(s$fit, X_test)
  preds_some <- predict(s$fit, X_test, idx_use = c(1, 5, 10))
  expect_equal(nrow(preds_some), 3)
  expect_equal(preds_some[1, ], preds_all[1, ])
  expect_equal(preds_some[2, ], preds_all[5, ])
  expect_equal(preds_some[3, ], preds_all[10, ])
})

test_that("predict.bppr() accepts a single test row", {
  s <- fit_small()
  preds <- predict(s$fit, matrix(runif(4), nrow = 1))
  expect_equal(ncol(preds), 1)
  expect_equal(nrow(preds), s$fit$n_keep)
})

test_that("predict.bppr() rejects non-numeric newdata", {
  s <- fit_small()
  bad <- matrix(letters[1:8], 2, 4)
  expect_error(predict(s$fit, bad), "numeric")
})

test_that("predict.bppr() rejects wrong number of columns", {
  s <- fit_small(p = 4)
  expect_error(predict(s$fit, matrix(runif(6), 2, 3)), "columns")
})

test_that("predict.bppr() validates idx_use", {
  s <- fit_small()
  X_test <- matrix(runif(4), nrow = 1)

  expect_error(predict(s$fit, X_test, idx_use = 0), "between 1")
  expect_error(predict(s$fit, X_test, idx_use = s$fit$n_keep + 1), "between 1")
  expect_error(predict(s$fit, X_test, idx_use = 1.5), "whole numbers")
  expect_error(predict(s$fit, X_test, idx_use = numeric(0)), "non-empty")
  expect_error(predict(s$fit, X_test, idx_use = NA_integer_), "missing")
  expect_error(predict(s$fit, X_test, idx_use = "1"), "numeric")
})

test_that("predict.bppr() is consistent across categorical ridge functions", {
  # Regression test for stale cached basis functions: consecutive draws with
  # the same ridge count but different categorical features must not reuse a
  # stale basis. Compare batched prediction with per-draw prediction.
  set.seed(21)
  n <- 80
  X <- cbind(matrix(runif(n * 3), n, 3),
             sample(0:1, n, replace = TRUE),
             sample(0:1, n, replace = TRUE))
  y <- X[, 1] + X[, 4] - X[, 5] + rnorm(n, sd = 0.1)
  fit <- bppr(X, y, n_post = 30, n_burn = 10, n_adapt = 0, print_every = 0)

  X_test <- cbind(matrix(runif(6 * 3), 6, 3),
                  sample(0:1, 6, replace = TRUE),
                  sample(0:1, 6, replace = TRUE))

  preds_batch <- predict(fit, X_test)
  # Predicting one draw at a time forces a fresh basis every time.
  preds_each <- t(vapply(1:fit$n_keep,
                         function(k) predict(fit, X_test, idx_use = k)[1, ],
                         numeric(nrow(X_test))))
  expect_equal(preds_batch, preds_each, tolerance = 1e-10)
})
