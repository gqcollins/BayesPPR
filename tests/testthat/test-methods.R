# Tests for the S3 methods print.bppr(), summary.bppr() and plot.bppr().

fit_for_methods <- function() {
  set.seed(5)
  n <- 60
  X <- matrix(runif(n * 3), n, 3)
  y <- X[, 1] + X[, 2]^2 + rnorm(n, sd = 0.1)
  bppr(X, y, n_post = 20, n_burn = 10, n_adapt = 0, print_every = 0)
}

test_that("print.bppr() reports the variable count and sample size", {
  fit <- fit_for_methods()
  out <- capture.output(print(fit))
  expect_true(any(grepl("Number of variables:\\s*3", out)))
  expect_true(any(grepl("Sample size:\\s*60", out)))
})

test_that("summary.bppr() reports ridge range and mean residual sd", {
  fit <- fit_for_methods()
  out <- capture.output(summary(fit))
  expect_true(any(grepl("Number of ridge functions", out)))
  expect_true(any(grepl("error sd", out)))
})

test_that("plot.bppr() runs without error and returns invisibly", {
  fit <- fit_for_methods()
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp)
  on.exit({ dev.off(); unlink(tmp) }, add = TRUE)

  expect_no_error(plot(fit))
  expect_no_error(plot(fit, pred = FALSE))
  expect_no_error(plot(fit, quants = NULL))
})

test_that("plot.bppr() rejects non-bppr input", {
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp)
  on.exit({ dev.off(); unlink(tmp) }, add = TRUE)
  expect_error(plot.bppr(list(a = 1)), "class bppr")
})

test_that("plot.bppr() restores graphical parameters on exit", {
  fit <- fit_for_methods()
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp)
  on.exit({ dev.off(); unlink(tmp) }, add = TRUE)

  before <- par("mfrow")
  plot(fit)
  expect_equal(par("mfrow"), before)
})
