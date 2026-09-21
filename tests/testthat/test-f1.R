test_that("f1 prediction matches previous validated values", {
  # Snapshot regression test: predictions from a stored fit must reproduce the
  # previously validated predictions bit-for-bit.
  eps <- 1e-12

  fit <- readRDS('../f1_fit.rda')       # previous model
  Xtest <- readRDS('../f1_Xtest.rda')   # x values
  oldpreds <- readRDS('../f1_preds.rda') # old predictions at x values
  newpreds <- predict(fit, Xtest)       # new predictions at x values

  diff <- max(abs(range(newpreds - oldpreds)))
  expect_lt(diff, eps)
})

test_that("predict.bppr() is deterministic (no randomness in prediction)", {
  fit <- readRDS('../f1_fit.rda')
  Xtest <- readRDS('../f1_Xtest.rda')
  # predict() draws no random numbers, so repeated calls must be identical.
  expect_identical(predict(fit, Xtest), predict(fit, Xtest))
})
