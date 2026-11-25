# Tests for bppr_sampler(): an updatable chain whose response can change.

sampler_data <- function(seed = 3, n = 80) {
  set.seed(seed)
  X <- matrix(runif(n * 3), n, 3)
  y <- X[, 1] + sin(4 * X[, 2]) + rnorm(n, sd = 0.1)
  list(X = X, y = y, n = n)
}

test_that("stepping a sampler by hand reproduces bppr()", {
  d <- sampler_data()
  set.seed(5)
  fit <- bppr(d$X, d$y, n_post = 40, n_burn = 0, print_every = 0)

  set.seed(5)
  s <- bppr_sampler(d$X, d$y)
  s$save()                              # iteration 1 is the initial state
  for (i in 2:40) {
    s$step(1)
    s$save()
  }
  fit2 <- s$as_bppr()
  expect_identical(fit2$n_ridge, fit$n_ridge)
  expect_identical(fit2$coefs, fit$coefs)
  expect_identical(fit2$sd_resid, fit$sd_resid)
  expect_identical(fit2$var_coefs, fit$var_coefs)
  expect_identical(predict(fit2, d$X), predict(fit, d$X))
})

test_that("set_y() matches a sampler rebuilt from the same state on the new y", {
  d <- sampler_data()
  set.seed(6)
  s <- bppr_sampler(d$X, d$y)
  s$step(100)
  expect_gt(s$n_ridge, 0)
  y2 <- d$y^2 - 3
  s$set_y(y2)

  s_new <- bppr_sampler(d$X, y2, bppr_init = s$current())
  K <- s$n_basis_total
  expect_equal(s$Bty[seq_len(K)], s_new$Bty[seq_len(K)], tolerance = 1e-12)
  expect_equal(c(s$qf_info$qf), c(s_new$qf_info$qf), tolerance = 1e-12)
  expect_equal(c(s$qf_info$ls_est), c(s_new$qf_info$ls_est), tolerance = 1e-10)
  expect_equal(s$ssy, s_new$ssy)
  expect_equal(c(s$sse), c(s_new$sse), tolerance = 1e-10)
  expect_equal(s$basis_mat, s_new$basis_mat)
})

test_that("the chain follows a changing response", {
  d <- sampler_data()
  set.seed(7)
  s <- bppr_sampler(d$X, d$y)
  s$step(200)
  s$set_y(d$y + 5)
  s$step(20)
  expect_equal(mean(s$fitted()), mean(d$y) + 5, tolerance = 0.05)
  expect_lt(sqrt(mean((s$fitted() - d$y - 5)^2)), 0.3)
})

test_that("s2_shape = s2_rate = 0 is the default prior, and other values change it", {
  d <- sampler_data()
  set.seed(8)
  a <- bppr(d$X, d$y, n_post = 30, n_burn = 30, print_every = 0)
  set.seed(8)
  b <- bppr(d$X, d$y, n_post = 30, n_burn = 30, print_every = 0, s2_shape = 0, s2_rate = 0)
  a$call <- b$call <- NULL
  expect_identical(a, b)

  # A strong prior centred far above the residual variance pulls s2 up
  set.seed(8)
  big <- bppr(d$X, d$y, n_post = 100, n_burn = 100, print_every = 0,
              s2_shape = 1000, s2_rate = 1000)
  expect_gt(mean(big$sd_resid^2), 10 * mean(a$sd_resid^2))
  expect_identical(big$s2_shape, 1000)
})

test_that("s2_fixed keeps the residual variance fixed", {
  d <- sampler_data()
  set.seed(9)
  fit <- bppr(d$X, d$y, n_post = 50, n_burn = 50, print_every = 0, s2_fixed = 0.02)
  expect_true(all(fit$sd_resid == sqrt(0.02)))
  expect_lt(sqrt(mean((colMeans(predict(fit, d$X)) - d$y)^2)), 0.3)
})

test_that("intercept = FALSE gives draws that average zero over the training x", {
  d <- sampler_data()
  yc <- d$y - mean(d$y)
  set.seed(10)
  s <- bppr_sampler(d$X, yc, intercept = FALSE)
  expect_identical(s$n_basis_total, 0L)
  expect_identical(s$fitted(), numeric(d$n))  # the empty model
  s$step(1)                                    # params drawn with no basis
  s$step(200)
  for (i in 1:20) {
    s$step(1)
    s$save()
    expect_lt(abs(mean(s$fitted())), 1e-10)
  }
  expect_gt(s$n_ridge, 0)
  fit <- s$as_bppr()
  expect_false(fit$intercept)
  expect_length(fit$basis_mn, 20)
  P <- predict(fit, d$X)
  expect_equal(P[20, ], s$fitted(), tolerance = 1e-10)
  expect_lt(max(abs(rowMeans(P))), 1e-10)
  expect_lt(sqrt(mean((colMeans(P) - yc)^2)), 0.3)

  # and through bppr(), including resume and the conditioning check
  set.seed(11)
  fit2 <- bppr(d$X, yc, n_post = 20, n_burn = 50, print_every = 0, intercept = FALSE)
  expect_lt(max(abs(rowMeans(predict(fit2, d$X)))), 1e-10)
  fit3 <- bppr_resume(fit2, n_post = 10, n_burn = 5, print_every = 0)
  expect_false(fit3$intercept)
  expect_lt(max(abs(rowMeans(predict(fit3, d$X)))), 1e-10)
})

test_that("death down to the empty model works without an intercept", {
  d <- sampler_data()
  set.seed(12)
  s <- bppr_sampler(d$X, d$X[, 1] - mean(d$X[, 1]), intercept = FALSE, n_ridge_mean = 0.5)
  s$step(50)
  expect_gt(s$n_ridge, 0)
  s$set_y(rnorm(d$n, sd = 0.01))          # nothing left to explain
  for (i in 1:300) s$step(1)
  expect_equal(s$n_ridge, 0)
  expect_length(s$coefs, 0)
  expect_identical(s$fitted(), numeric(d$n))
})

test_that("bppr_resume() takes a new y, through the live sampler or a rebuild", {
  d <- sampler_data()
  y2 <- d$y + 3
  set.seed(13)
  fit <- bppr(d$X, d$y, n_post = 20, n_burn = 100, print_every = 0, keep_sampler = TRUE)
  expect_s3_class(fit$sampler, "bppr_sampler")

  # live sampler: no rebuild, the adapted weights carry on
  r1 <- bppr_resume(fit, y = y2, n_post = 20, n_burn = 20, print_every = 0)
  expect_identical(r1$sampler, fit$sampler)
  expect_equal(r1$n_keep, 20)
  expect_identical(r1$y, y2)
  expect_equal(mean(colMeans(predict(r1, d$X))), mean(y2), tolerance = 0.05)
  r1a <- bppr_resume(r1, y = y2, append = TRUE, n_post = 10, n_burn = 0, print_every = 0)
  expect_equal(r1a$n_keep, 30)
  expect_length(r1a$coefs, 30)

  # rebuild from the last draw
  fit$sampler <- NULL
  r2 <- bppr_resume(fit, y = y2, n_post = 20, n_burn = 20, print_every = 0)
  expect_null(r2$sampler)
  expect_identical(r2$y, y2)
  expect_equal(mean(colMeans(predict(r2, d$X))), mean(y2), tolerance = 0.05)
})

test_that("a resumed chain starts with the right birth/death proposal ratio", {
  d <- sampler_data()
  set.seed(14)
  s <- bppr_sampler(d$X, d$y)
  s$step(100)
  expect_gt(s$n_ridge, 0)
  s2 <- bppr_sampler(d$X, d$y, bppr_init = s$current())
  expect_identical(s2$log_mh_bd, s$log_mh_bd)
})

test_that("bppr_sampler() validates its new arguments", {
  d <- sampler_data()
  expect_error(bppr_sampler(d$X, d$y, intercept = NA), "intercept")
  expect_error(bppr_sampler(d$X, d$y, s2_shape = -1), "non-negative")
  expect_error(bppr_sampler(d$X, d$y, s2_fixed = 0), "s2_fixed")
  s <- bppr_sampler(d$X, d$y)
  expect_error(s$set_y(d$y[-1]), "length")
  expect_output(print(s), "BayesPPR sampler")
})
