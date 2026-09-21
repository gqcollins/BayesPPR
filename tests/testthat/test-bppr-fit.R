# Tests for the main bppr() fitting function: structure of the returned
# object, input handling, and argument validation.

make_data <- function(n = 60, p = 4, seed = 1) {
  set.seed(seed)
  X <- matrix(runif(n * p), n, p)
  y <- 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 + rnorm(n, sd = 0.1)
  list(X = X, y = y)
}

test_that("bppr() returns a well-formed bppr object", {
  d <- make_data()
  fit <- bppr(d$X, d$y, n_post = 20, n_burn = 10, n_adapt = 5, print_every = 0)

  expect_s3_class(fit, "bppr")
  expect_equal(fit$n_keep, 20)
  expect_length(fit$n_ridge, 20)
  expect_length(fit$sd_resid, 20)
  expect_length(fit$coefs, 20)
  expect_length(fit$feat, 20)
  expect_length(fit$proj_dir, 20)
  expect_length(fit$knots, 20)

  # Recorded data should match the inputs.
  expect_equal(fit$X, d$X)
  expect_equal(fit$y, d$y)
  expect_equal(nrow(fit$X), length(fit$y))

  # Residual SDs are positive; number of ridge functions is a non-negative count.
  expect_true(all(fit$sd_resid > 0))
  expect_true(all(fit$n_ridge >= 0))
  expect_true(all(fit$n_ridge <= fit$n_ridge_max))
})

test_that("bppr() honours thinning", {
  d <- make_data()
  fit <- bppr(d$X, d$y, n_post = 20, n_burn = 5, n_adapt = 0, n_thin = 5,
              print_every = 0)
  expect_equal(fit$n_keep, 4)          # 20 / 5
  expect_length(fit$n_ridge, 4)
  expect_length(fit$coefs, 4)
})

test_that("bppr() rounds n_post down to a multiple of n_thin", {
  d <- make_data()
  fit <- bppr(d$X, d$y, n_post = 23, n_burn = 5, n_adapt = 0, n_thin = 5,
              print_every = 0)
  expect_equal(fit$n_post, 20)         # 23 - 23 %% 5
  expect_equal(fit$n_keep, 4)
})

test_that("bppr() accepts a data frame for X", {
  d <- make_data()
  df <- as.data.frame(d$X)
  expect_no_error(
    fit <- bppr(df, d$y, n_post = 10, n_burn = 5, n_adapt = 0, print_every = 0)
  )
  expect_s3_class(fit, "bppr")
})

test_that("bppr() supports the flat prior", {
  d <- make_data()
  fit <- bppr(d$X, d$y, prior_coefs = "flat", n_post = 15, n_burn = 5,
              n_adapt = 0, print_every = 0)
  expect_s3_class(fit, "bppr")
  # With the flat prior, burn-in is folded into the adapt phase.
  expect_equal(fit$n_burn, 0)
  expect_true(all(is.na(fit$var_coefs)))
})

test_that("bppr() rejects non-numeric X", {
  d <- make_data()
  Xchar <- d$X
  storage.mode(Xchar) <- "character"
  expect_error(bppr(Xchar, d$y, print_every = 0), "numeric")
})

test_that("bppr() rejects mismatched X/y dimensions", {
  d <- make_data(n = 60)
  expect_error(bppr(d$X, d$y[1:59], print_every = 0), "nrow\\(X\\)")
})

test_that("bppr() rejects missing values", {
  d <- make_data()
  Xna <- d$X
  Xna[1, 1] <- NA
  expect_error(bppr(Xna, d$y, print_every = 0), "missing")

  yna <- d$y
  yna[1] <- NA
  expect_error(bppr(d$X, yna, print_every = 0), "missing")
})

test_that("bppr() validates n_thin against n_post", {
  d <- make_data()
  expect_error(bppr(d$X, d$y, n_post = 10, n_thin = 20, print_every = 0),
               "n_thin")
})

test_that("bppr() validates prior_coefs", {
  d <- make_data()
  expect_error(
    bppr(d$X, d$y, prior_coefs = "bogus", n_post = 5, n_burn = 0, n_adapt = 0,
         print_every = 0),
    "prior_coefs"
  )
})

test_that("bppr() validates w_feat and w_n_act lengths", {
  d <- make_data(p = 4)
  expect_error(bppr(d$X, d$y, w_feat = rep(1, 3), print_every = 0), "w_feat")
  expect_error(
    bppr(d$X, d$y, n_act_max = 2, w_n_act = rep(1, 5), print_every = 0),
    "w_n_act"
  )
})

test_that("bppr() validates scale_proj_dir_prop range", {
  d <- make_data()
  expect_error(bppr(d$X, d$y, scale_proj_dir_prop = 0, print_every = 0),
               "scale_proj_dir_prop")
  expect_error(bppr(d$X, d$y, scale_proj_dir_prop = 1.5, print_every = 0),
               "scale_proj_dir_prop")
})

test_that("bppr() handles a constant column with few predictors", {
  # Regression test: fitting used to fail when X had a constant column and
  # ncol(X) <= 3.
  set.seed(7)
  n <- 60
  X <- cbind(runif(n), runif(n), rep(1, n))   # third column constant
  y <- X[, 1] + rnorm(n, sd = 0.1)
  expect_no_error(
    fit <- bppr(X, y, n_post = 10, n_burn = 5, n_adapt = 0, print_every = 0)
  )
  # A constant column carries zero feature weight.
  expect_equal(fit$w_feat[3], 0)
})

test_that("bppr() errors when all features are constant", {
  n <- 40
  X <- cbind(rep(1, n), rep(2, n))
  y <- rnorm(n)
  expect_error(bppr(X, y, print_every = 0), "usable features")
})

test_that("bppr() warns and caps n_act_max above the usable feature count", {
  d <- make_data(p = 3)
  expect_warning(
    fit <- bppr(d$X, d$y, n_act_max = 10, n_post = 5, n_burn = 0, n_adapt = 0,
                print_every = 0),
    "n_act_max"
  )
  expect_lte(fit$n_act_max, 3)
})

test_that("bppr() handles categorical features", {
  set.seed(11)
  n <- 80
  X <- cbind(matrix(runif(n * 3), n, 3),
             sample(0:1, n, replace = TRUE),
             sample(0:1, n, replace = TRUE))
  y <- X[, 1] + X[, 4] - X[, 5] + rnorm(n, sd = 0.1)
  expect_no_error(
    fit <- bppr(X, y, n_post = 20, n_burn = 10, n_adapt = 0, print_every = 0)
  )
  expect_s3_class(fit, "bppr")
})

test_that("bppr() is silent when print_every = 0", {
  d <- make_data()
  expect_silent(
    bppr(d$X, d$y, n_post = 10, n_burn = 5, n_adapt = 0, print_every = 0)
  )
})
