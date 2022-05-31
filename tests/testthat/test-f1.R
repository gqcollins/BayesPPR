test_that("f1 prediction matches previous validated values", {
  cat('Modified Friedman function test (prediction) \n')

  eps <- 1e-12

  fit <- readRDS('../f1_fit.rda') # previous model
  Xtest <- readRDS('../f1_Xtest.rda') # x values
  oldpreds <- readRDS('../f1_preds.rda') # old predictions at x values
  newpreds <- predict(fit, Xtest, n_ignore = 9000) # new predictions at x values
  diff <- max(abs(range(newpreds - oldpreds))) # difference between new and old predictions
  expect_that(diff, is_less_than(eps))
})
