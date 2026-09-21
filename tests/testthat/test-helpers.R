# Unit tests for internal helper functions. These are pure/deterministic and
# so can be checked against exact expected values.

# --- relu() ------------------------------------------------------------------

test_that("relu() clamps negatives to zero and passes positives through", {
  relu <- BayesPPR:::relu
  expect_equal(relu(c(-3, -1, 0, 2, 5)), c(0, 0, 0, 2, 5))
  expect_equal(relu(0), 0)
})

# --- get_cat_basis() ---------------------------------------------------------

test_that("get_cat_basis() returns the single column unchanged for one feature", {
  get_cat_basis <- BayesPPR:::get_cat_basis
  x <- matrix(c(0, 1, 0, 1), ncol = 1)
  expect_equal(get_cat_basis(x), x)
})

test_that("get_cat_basis() encodes the 'any level active' interaction", {
  get_cat_basis <- BayesPPR:::get_cat_basis
  # basis = 1 - prod(1 - x_j): equals 0 only when every column is 0.
  X <- rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1))
  expect_equal(as.vector(get_cat_basis(X)), c(0, 1, 1, 1))
})

# --- get_mns_basis() ---------------------------------------------------------

test_that("get_mns_basis() has the expected shape and relu first column", {
  get_mns_basis <- BayesPPR:::get_mns_basis
  u <- seq(-1, 1, length.out = 20)
  # df_spline = 4 -> knot_quants gives df_spline + 1 = 5 knots, plus knot0 = 6.
  knots <- c(-0.9, quantile(u[u > -0.9], seq(0, 1, length.out = 5)))
  basis <- get_mns_basis(u, knots)

  expect_true(is.matrix(basis))
  expect_equal(nrow(basis), length(u))
  # First column is relu(u - knots[1]); non-negative and zero below the knot.
  expect_true(all(basis[, 1] >= 0))
  expect_true(all(basis[u < knots[1], 1] == 0))
})

test_that("get_mns_basis() with a single non-trivial knot is plain relu", {
  get_mns_basis <- BayesPPR:::get_mns_basis
  u <- c(-1, 0, 1, 2)
  # Only two knots -> df = 0, so basis is just relu(u - knots[1]).
  basis <- get_mns_basis(u, c(0, 5))
  expect_equal(as.vector(basis), c(0, 0, 1, 2))
})

# --- get_move_type() ---------------------------------------------------------

test_that("get_move_type() forces birth from an empty model", {
  get_move_type <- BayesPPR:::get_move_type
  expect_equal(get_move_type(0, 0, 10), "birth")
})

test_that("get_move_type() forces death from a full model with no quant ridges", {
  get_move_type <- BayesPPR:::get_move_type
  expect_equal(get_move_type(10, 0, 10), "death")
})

test_that("get_move_type() only proposes valid moves", {
  get_move_type <- BayesPPR:::get_move_type
  set.seed(1)
  # Full model, some quantitative ridges: death or change only (no birth).
  moves <- replicate(200, get_move_type(10, 3, 10))
  expect_true(all(moves %in% c("death", "change")))

  # Interior model, no quantitative ridges: birth or death only (no change).
  moves <- replicate(200, get_move_type(5, 0, 10))
  expect_true(all(moves %in% c("birth", "death")))

  # Interior model with quantitative ridges: all three possible.
  moves <- replicate(500, get_move_type(5, 2, 10))
  expect_true(all(moves %in% c("birth", "death", "change")))
})

# --- get_log_mh_bd() ---------------------------------------------------------

test_that("get_log_mh_bd() matches the number of available move types", {
  get_log_mh_bd <- BayesPPR:::get_log_mh_bd
  expect_equal(get_log_mh_bd(0, 0, 10), 0)          # empty
  expect_equal(get_log_mh_bd(10, 0, 10), 0)         # full, no quant
  expect_equal(get_log_mh_bd(10, 2, 10), log(2))    # full, with quant
  expect_equal(get_log_mh_bd(5, 0, 10), log(2))     # interior, no quant
  expect_equal(get_log_mh_bd(5, 2, 10), log(3))     # interior, with quant
})

# --- check_idx_use() ---------------------------------------------------------

test_that("check_idx_use() accepts valid indices and coerces to integer", {
  check_idx_use <- BayesPPR:::check_idx_use
  expect_equal(check_idx_use(c(1, 3, 5), 10), c(1L, 3L, 5L))
  expect_type(check_idx_use(c(1, 2), 10), "integer")
  expect_equal(check_idx_use(c(2, 4), 5), c(2L, 4L))
})

test_that("check_idx_use() rejects invalid indices", {
  check_idx_use <- BayesPPR:::check_idx_use
  expect_error(check_idx_use(numeric(0), 10), "non-empty")
  expect_error(check_idx_use("a", 10), "numeric")
  expect_error(check_idx_use(NA_real_, 10), "missing")
  expect_error(check_idx_use(1.5, 10), "whole numbers")
  expect_error(check_idx_use(0, 10), "between 1")
  expect_error(check_idx_use(11, 10), "between 1")
})

# --- get_qf_info() / append_qf_inv_chol() ------------------------------------

test_that("get_qf_info() computes the quadratic form for a well-conditioned system", {
  get_qf_info <- BayesPPR:::get_qf_info
  set.seed(2)
  B <- cbind(1, matrix(rnorm(30 * 2), 30, 2))
  y <- rnorm(30)
  BtB <- t(B) %*% B
  Bty <- t(B) %*% y

  info <- get_qf_info(BtB, Bty)
  expect_false(is.null(info))
  # Least-squares estimate matches lm without an extra intercept.
  ls <- solve(BtB, Bty)
  expect_equal(as.vector(info$ls_est), as.vector(ls), tolerance = 1e-8)
  # qf = Bty' (BtB)^{-1} Bty.
  expect_equal(as.vector(info$qf),
               as.vector(t(Bty) %*% ls), tolerance = 1e-8)
})

test_that("get_qf_info() returns NULL for a rank-deficient / ill-conditioned system", {
  get_qf_info <- BayesPPR:::get_qf_info
  # Duplicated column -> singular BtB -> chol fails -> NULL.
  B <- cbind(1, c(1, 2, 3), c(1, 2, 3))
  BtB <- t(B) %*% B
  Bty <- t(B) %*% c(1, 2, 3)
  expect_null(get_qf_info(BtB, Bty))
})

test_that("append_qf_inv_chol() gives the inverse of the Cholesky factor", {
  get_qf_info <- BayesPPR:::get_qf_info
  append_qf_inv_chol <- BayesPPR:::append_qf_inv_chol
  set.seed(3)
  B <- cbind(1, matrix(rnorm(20 * 2), 20, 2))
  BtB <- t(B) %*% B
  Bty <- t(B) %*% rnorm(20)
  info <- append_qf_inv_chol(get_qf_info(BtB, Bty), dim = 3)
  # chol %*% inv_chol should be the identity.
  expect_equal(info$chol %*% info$inv_chol, diag(3), tolerance = 1e-8)
})

# --- rps(): power-spherical draws --------------------------------------------

test_that("rps() returns a unit vector of the requested dimension", {
  # rps() is only ever invoked with d >= 2 (the d == 1 case is handled by a
  # direct sample of c(-1, 1) in bppr()).
  rps <- BayesPPR:::rps
  set.seed(4)
  for (d in c(2, 3, 5)) {
    x <- rps(rep(1 / sqrt(d), d), 0)
    expect_equal(length(x), d)
    expect_equal(sqrt(sum(x^2)), 1, tolerance = 1e-8)
  }
})

test_that("rps() with large kappa concentrates near the mean direction", {
  rps <- BayesPPR:::rps
  set.seed(5)
  # Use a non-axis mean direction: mu == e1 triggers a degenerate reflection.
  mu <- c(1, 1, 1) / sqrt(3)
  draws <- replicate(200, sum(mu * rps(mu, 1e4)))
  # Projection onto mu should be very close to 1 for large concentration.
  expect_true(mean(draws) > 0.99)
})

# --- dwallenius() ------------------------------------------------------------

test_that("dwallenius() returns a positive probability", {
  dwallenius <- BayesPPR:::dwallenius
  w <- rep(1, 5)
  expect_true(dwallenius(w, c(1, 2)) > 0)
  expect_true(dwallenius(w, 3) > 0)
})

test_that("dwallenius() is symmetric under equal weights", {
  dwallenius <- BayesPPR:::dwallenius
  w <- rep(1, 6)
  # Under equal weights, any two feature sets of the same size are exchangeable.
  expect_equal(dwallenius(w, c(1, 2)), dwallenius(w, c(4, 5)), tolerance = 1e-10)
  expect_equal(dwallenius(w, c(1, 2, 3)), dwallenius(w, c(4, 5, 6)),
               tolerance = 1e-10)
})

# --- get_n_cores_max() -------------------------------------------------------

test_that("get_n_cores_max() returns a positive integer", {
  n <- BayesPPR:::get_n_cores_max()
  expect_true(is.numeric(n))
  expect_true(n >= 1)
})
