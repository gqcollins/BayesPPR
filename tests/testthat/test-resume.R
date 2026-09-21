# Tests for bppr_resume(): continuing a chain, appending, and input validation.

fit_to_resume <- function(seed = 9, ...) {
  set.seed(seed)
  n <- 60
  X <- matrix(runif(n * 3), n, 3)
  y <- X[, 1] + X[, 2]^2 + rnorm(n, sd = 0.1)
  bppr(X, y, n_post = 20, n_burn = 10, n_adapt = 0, print_every = 0, ...)
}

test_that("bppr_resume() returns a bppr object with only the new draws by default", {
  fit <- fit_to_resume()
  fit2 <- bppr_resume(fit, n_post = 15, n_burn = 5, n_adapt = 0, print_every = 0)

  expect_s3_class(fit2, "bppr")
  expect_equal(fit2$n_keep, 15)
  expect_length(fit2$n_ridge, 15)
})

test_that("bppr_resume() appends new draws when append = TRUE", {
  fit <- fit_to_resume()
  fit2 <- bppr_resume(fit, append = TRUE, n_post = 15, n_burn = 5, n_adapt = 0,
                      print_every = 0)

  expect_equal(fit2$n_keep, fit$n_keep + 15)
  expect_length(fit2$n_ridge, fit$n_keep + 15)
  expect_length(fit2$coefs, fit$n_keep + 15)
  # The first draws should be the original ones carried through unchanged.
  expect_equal(fit2$n_ridge[seq_len(fit$n_keep)], fit$n_ridge)
})

test_that("bppr_resume() works with the flat prior", {
  fit <- fit_to_resume(prior_coefs = "flat")
  expect_no_error(
    fit2 <- bppr_resume(fit, n_post = 10, n_burn = 5, n_adapt = 0,
                        print_every = 0)
  )
  expect_s3_class(fit2, "bppr")
})

test_that("bppr_resume() handles a categorical ridge in the last draw", {
  # Regression test: resuming a fit whose last draw contained a categorical
  # ridge function used to corrupt basis bookkeeping and fail.
  set.seed(31)
  n <- 90
  X <- cbind(matrix(runif(n * 3), n, 3),
             sample(0:1, n, replace = TRUE),
             sample(0:1, n, replace = TRUE))
  y <- X[, 1] + X[, 4] - X[, 5] + rnorm(n, sd = 0.1)
  fit <- bppr(X, y, n_post = 40, n_burn = 10, n_adapt = 0, print_every = 0)

  expect_no_error(
    fit2 <- bppr_resume(fit, n_post = 10, n_burn = 5, n_adapt = 0,
                        print_every = 0)
  )
  expect_s3_class(fit2, "bppr")
})

test_that("bppr_resume() rejects objects of the wrong class", {
  expect_error(bppr_resume(list(a = 1), print_every = 0),
               "class 'bppr' or 'bppr_pca'")
})
