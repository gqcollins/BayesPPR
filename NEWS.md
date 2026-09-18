# BayesPPR 0.2.0

# BayesPPR 0.1.0.9000

* Fixed corrupted basis bookkeeping when resuming a fit whose last draw
  contained a categorical ridge function (`bppr_resume()` failed outright).
* Fixed stale cached basis functions in `predict.bppr()`, which silently
  returned incorrect predictions when consecutive draws shared the same number
  of ridge functions but used different features for a categorical ridge.
* Fixed fitting failures when `X` contains a constant column and `ncol(X) <= 3`.
* `bppr()` now accepts a data frame for `X`, as documented.
* `bppr_pca()` and `predict.bppr_pca()` no longer fail when
  `parallel::detectCores()` returns `NA`, and now validate `par_type`.
* `pca_setup()` now bounds `n_pc` by `min(nrow(Y), ncol(Y))` instead of the
  maximum, and handles `prop_var = 1`.
* `predict.bppr_pca()` no longer fails when `newdata` has a single row.
* Added validation of `idx_use`, `newdata` dimensions, `w_feat`, `w_n_act`,
  and missing values in `X`/`y`.

# BayesPPR 0.1.0

* Initial CRAN submission.
