########################################################################
## Persistent, updatable BayesPPR sampler
########################################################################

#' @title Updatable BayesPPR Sampler
#'
#' @description Sets up the state of a BayesPPR Markov chain so that it can be
#'   advanced a few iterations at a time, with a response that may change
#'   between calls. This is what \link{bppr} runs internally; use it directly
#'   when BayesPPR is one block of a larger Gibbs sampler whose response is
#'   itself a latent variable that is redrawn every outer iteration.
#'
#' @inheritParams bppr
#' @param intercept logical; if \code{TRUE} (default) the model has an
#'   intercept. If \code{FALSE} there is no intercept and the basis functions of
#'   every ridge function are centred at their means over the training inputs,
#'   so every draw of \eqn{f} averages exactly zero over the training inputs.
#' @param s2_shape,s2_rate shape and rate of the inverse-gamma prior on the
#'   residual variance. The default \code{0, 0} is the improper prior
#'   \eqn{1/\sigma^2} used by \link{bppr}.
#' @param s2_fixed if not \code{NULL}, the residual variance is fixed at this
#'   value and not sampled.
#' @param s2_init residual variance before the first draw.
#'
#' @details The chain state (the ridge functions, the basis matrix and its
#'   cross-products, and the current coefficients and variances) is kept in an
#'   environment. Only \eqn{B'y}, \eqn{y'y} and the quantities derived from them
#'   depend on the response, so \code{$set_y()} costs \eqn{O(nK)} for \eqn{K}
#'   basis functions and does not rebuild the basis.
#'
#'   The returned environment has these methods:
#'   \describe{
#'     \item{\code{set_y(y)}}{replace the response.}
#'     \item{\code{iterate(sample_params = TRUE)}}{one RJMCMC move (birth,
#'       death or change of a ridge function), followed, if
#'       \code{sample_params}, by Gibbs draws of the residual variance, the
#'       coefficients and their prior variance.}
#'     \item{\code{step(k = 1)}}{\code{k} calls of \code{iterate()}.}
#'     \item{\code{fitted()}}{the current mean function at the training inputs.}
#'     \item{\code{s2()}}{the current residual variance.}
#'     \item{\code{save()}}{append the current parameters to the stored draws.}
#'     \item{\code{n_saved()}}{number of stored draws.}
#'     \item{\code{accept_rate()}}{fraction of accepted RJMCMC proposals.}
#'     \item{\code{current()}}{the current state, as a list usable as
#'       \code{bppr_init}.}
#'     \item{\code{as_bppr(...)}}{the stored draws as a \code{"bppr"} object,
#'       usable with \link{predict.bppr}; \code{...} are extra fields.}
#'   }
#' @return An environment of class \code{"bppr_sampler"}.
#' @seealso \link{bppr}, \link{bppr_resume}.
#' @export
#' @examples
#' set.seed(1)
#' n <- 100
#' X <- matrix(runif(n * 2), n, 2)
#' s <- bppr_sampler(X, X[, 1]^2 + rnorm(n, sd = 0.1))
#' s$step(200)
#' for(it in 1:50){
#'   s$set_y(X[, 1]^2 + rnorm(n, sd = 0.1)) # a new response every iteration
#'   s$step(1)
#'   s$save()
#' }
#' fit <- s$as_bppr()
#' dim(predict(fit, X))
#'
bppr_sampler <- function(X, y, n_ridge_mean = 10, n_ridge_max = NULL, n_act_max = NULL, df_spline = 4,
                         prob_relu = 2/3, prior_coefs = "zs", shape_var_coefs = NULL, scale_var_coefs = NULL,
                         n_dat_min = NULL, scale_proj_dir_prop = NULL, adapt_act_feat = TRUE, w_n_act = NULL,
                         w_feat = NULL, intercept = TRUE, s2_shape = 0, s2_rate = 0, s2_fixed = NULL,
                         s2_init = 1, bppr_init = NULL){
  self <- environment()

  # Pre-processing
  X <- as.matrix(X)
  if(!is.numeric(X)){
    stop('X must be numeric. Categorical features should be coded as numeric.')
  }
  y <- as.numeric(y)
  n <- length(y)
  p <- ncol(X)
  if(nrow(X) != n){
    stop('nrow(X) must equal length(y)')
  }
  if(anyNA(X)  ||  anyNA(y)){
    stop('X and y must not contain missing values')
  }
  if(!is.logical(intercept)  ||  length(intercept) != 1  ||  is.na(intercept)){
    stop('intercept must be TRUE or FALSE')
  }
  if(s2_shape < 0  ||  s2_rate < 0){
    stop('s2_shape and s2_rate must be non-negative')
  }
  if(!is.null(s2_fixed)  &&  !(is.numeric(s2_fixed)  &&  length(s2_fixed) == 1  &&  s2_fixed > 0)){
    stop('s2_fixed must be NULL or a positive number')
  }
  jeffreys <- s2_shape == 0  &&  s2_rate == 0 # the improper 1/s2 prior of bppr()

  if(is.null(w_feat)){
    w_feat <- rep(1, p)
  }else if(length(w_feat) != p){
    stop('w_feat must have length ncol(X)')
  }

  mn_X <- sd_X <- numeric(p)
  X_st <- X
  feat_type <- character(p)
  for(j in 1:p){
    n_unique <- length(unique(X[, j]))
    if(n_unique <= 2){
      mn_X[j] <- 0
      sd_X[j] <- 1
      if(n_unique == 1){
        feat_type[j] <- ''
        w_feat[j] <- 0
      }else{
        feat_type[j] <- 'cat'
      }
    }else{
      mn_X[j] <- mean(X[, j])
      sd_X[j] <- sd(X[, j])
      X_st[, j] <- (X[, j] - mn_X[j]) / sd_X[j]
      if(n_unique <= df_spline){
        feat_type[j] <- 'disc'
      }else{
        feat_type[j] <- 'cont'
      }
    }
  }

  # Features with zero weight (e.g. constant columns) can never be proposed
  feat_usable <- which(w_feat > 0)
  p_use <- length(feat_usable)
  if(p_use == 0){
    stop('no usable features in X. All features are constant or have zero weight in w_feat.')
  }

  if(is.null(n_act_max)){
    n_cat <- sum(feat_type[feat_usable] == 'cat')
    n_act_max <- min(3, p_use - n_cat) + min(3, ceiling(n_cat/2))
  }else if(n_act_max > p_use){
    warning(paste0('n_act_max cannot exceed the number of usable features. Setting n_act_max = ', p_use))
  }
  n_act_max <- min(n_act_max, p_use)

  if(is.null(w_n_act)){
    w_n_act <- rep(1, n_act_max)
  }else if(length(w_n_act) != n_act_max){
    stop('w_n_act must have length n_act_max')
  }

  if(is.null(scale_proj_dir_prop)){
    proj_dir_prop_prec <- 1000 # scale_proj_dir_prop = 0.002
  }else if(scale_proj_dir_prop > 1  ||  scale_proj_dir_prop <= 0){
    stop("scale_proj_dir_prop must be in (0, 1]")
  }else{
    proj_dir_prop_prec <- 1/scale_proj_dir_prop
    proj_dir_prop_prec <- (proj_dir_prop_prec - 1) + sqrt(proj_dir_prop_prec * (proj_dir_prop_prec - 1))
  }

  if(is.null(n_dat_min)){
    n_dat_min <- min(20, 0.1 * n)
  }
  if(n_dat_min <= df_spline){
    warning('n_dat_min too small. If n_dat_min was set by default, df_spline is large compared to the sample size. Setting n_dat_min = df_spline + 1')
    n_dat_min <- df_spline + 1
  }
  p_dat_max <- 1 - n_dat_min / n # Maximum proportion of inactive datapoints in each ridge function

  if(is.null(n_ridge_max)){
    n_ridge_max <- min(150, floor(n/df_spline) - 2)
  }
  if(n_ridge_max <= 0){
    stop('n_ridge_max <= 0. If n_ridge_max was set by default, df_spline is too large compared to the sample size.')
  }

  knot_quants <- seq(0, 1, length.out = df_spline + 1) # Quantiles for knot locations (except initial knot)

  if(prior_coefs == 'zs'){
    if(is.null(shape_var_coefs)){
      shape_var_coefs <- 0.5
    }
    if(is.null(scale_var_coefs)){
      scale_var_coefs <- n/2
    }
  }else if(prior_coefs == 'flat'){
    var_coefs <- NA
    c_var_coefs <- 1
  }else{
    stop("prior_coefs must be either 'zs' or 'flat'")
  }
  proj_dir_mn <- lapply(1:n_act_max, function(a) rep(1/sqrt(a), a)) # prior mean for proj_dir (arbitrary, since prior precision is zero)

  # Basis bookkeeping. Segment 1 is the intercept (empty without one).
  n_int <- as.integer(intercept)
  n_basis_ridge <- n_int # Number of basis functions in each segment
  ridge_type <- character()
  j_quant <- integer() # which ridge functions are quantitative?
  n_quant <- 0 # how many ridge functions are quantitative?
  basis_mat <- matrix(1, n, n_int) # Current basis matrix
  basis_idx <- list(seq_len(n_int)) # Current indices of segments of basis functions
  BtB <- matrix(nrow = n_ridge_max * df_spline + 1, ncol = n_ridge_max * df_spline + 1)
  Bty <- numeric(n_ridge_max * df_spline + 1)

  # Without an intercept each ridge's basis is centred at its training means
  centre <- function(B){
    if(intercept) return(list(B = B, mn = NULL))
    mn <- colMeans(B)
    list(B = sweep(B, 2, mn), mn = mn)
  }

  # Initialization
  basis_mn <- list() # column means of each ridge's basis (intercept = FALSE only)
  if(is.null(bppr_init)){
    if(prior_coefs == 'zs'){
      var_coefs <- scale_var_coefs / shape_var_coefs
      c_var_coefs <- var_coefs / (var_coefs + 1)
    }
    sd_resid <- if(is.null(s2_fixed)) sqrt(s2_init) else sqrt(s2_fixed)
    coefs <- if(intercept) mean(y) else numeric(0)
    n_ridge <- 0
    n_act <- numeric() # Number of active features for jth ridge function
    feat <- list() # Features being used in jth ridge function
    knots <- list() # Location of knots for nsplines
    proj_dir <- list() # Ridge directions
  }else{
    if(prior_coefs == 'zs'){
      var_coefs <- bppr_init$var_coefs
      c_var_coefs <- var_coefs / (var_coefs + 1)
    }
    sd_resid <- if(is.null(s2_fixed)) bppr_init$sd_resid else sqrt(s2_fixed)
    coefs <- bppr_init$coefs
    n_ridge <- bppr_init$n_ridge
    n_act <- bppr_init$n_act
    feat <- bppr_init$feat
    knots <- bppr_init$knots
    proj_dir <- bppr_init$proj_dir
    if(n_ridge > 0){
      basis_idx_start <- n_int + 1
      for(j in 1:n_ridge){
        if(any(feat_type[feat[[j]]] == 'cont')){
          ridge_type[j] <- 'cont'
          n_basis_ridge[j + 1] <- df_spline
          j_quant <- c(j_quant, j)
          proj <- X_st[, feat[[j]], drop = FALSE] %*% proj_dir[[j]]
          rb <- centre(get_mns_basis(proj, knots[[j]])) # Get basis function
          basis_idx[[j + 1]] <- basis_idx_start:(basis_idx_start + df_spline - 1)
          basis_idx_start <- basis_idx_start + df_spline
        }else{
          n_basis_ridge[j + 1] <- 1
          if(any(feat_type[feat[[j]]] == 'disc')){
            ridge_type[j] <- 'disc'
            j_quant <- c(j_quant, j)
            rb <- centre(X_st[, feat[[j]], drop = FALSE] %*% proj_dir[[j]])
          }else{
            ridge_type[j] <- 'cat'
            rb <- centre(get_cat_basis(X_st[, feat[[j]], drop = FALSE]))
          }
          basis_idx[[j + 1]] <- basis_idx_start
          basis_idx_start <- basis_idx_start + 1
        }
        basis_mat <- cbind(basis_mat, rb$B)
        if(!intercept) basis_mn[[j]] <- rb$mn
      }
      n_quant <- length(j_quant)
    }
    if(!is.null(bppr_init$w_n_act)) w_n_act <- bppr_init$w_n_act
    if(!is.null(bppr_init$w_feat)) w_feat <- bppr_init$w_feat
  }

  n_basis_total <- sum(n_basis_ridge)
  idx_K <- seq_len(n_basis_total)
  BtB[idx_K, idx_K] <- t(basis_mat) %*% basis_mat # inner product of basis_mat with itself
  Bty[idx_K] <- t(basis_mat) %*% y # inner product of basis_mat with y
  qf_info <- get_qf_info(BtB[idx_K, idx_K, drop = FALSE], Bty[idx_K], intercept)
  if(is.null(qf_info)){
    stop('the initial basis matrix is singular or badly conditioned')
  }

  ssy <- c(t(y) %*% y) # Keep track of overall sse
  sse <- ssy - c_var_coefs * qf_info$qf
  log_mh_bd <- if(n_ridge > 0) get_log_mh_bd(n_ridge, n_quant, n_ridge_max) else 0
  log_mh_act_feat <- 0
  n_prop <- 0
  n_acc <- 0

  # Kept draws
  keep <- list(n_ridge = numeric(), n_act = list(), feat = list(), proj_dir = list(), knots = list(),
               coefs = list(), sd_resid = numeric(), var_coefs = numeric(), basis_mn = list())

  # Log ratio of marginal likelihoods (coefficients integrated out) between
  # residual sums of squares sse_prop and sse
  log_ml_ratio <- function(sse_prop, sse){
    if(!is.null(s2_fixed)){
      -(sse_prop - sse) / (2 * s2_fixed)
    }else if(jeffreys){
      -n/2 * (log(sse_prop) - log(sse))
    }else{
      -(n/2 + s2_shape) * (log(s2_rate + sse_prop/2) - log(s2_rate + sse/2))
    }
  }

  set_y <- function(y_new){
    y_new <- as.numeric(y_new)
    if(length(y_new) != n) stop('y must have length nrow(X)')
    if(anyNA(y_new)) stop('y must not contain missing values')
    y <<- y_new
    idx_K <- seq_len(n_basis_total)
    Bty[idx_K] <<- t(basis_mat) %*% y_new
    qf_info <<- update_qf_info(qf_info, Bty[idx_K])
    ssy <<- c(t(y_new) %*% y_new)
    sse <<- ssy - c_var_coefs * qf_info$qf
    invisible(self)
  }

  # Proposed basis for a quantitative ridge: projection, knots and basis
  propose_quant_basis <- function(feat_prop, proj_dir_prop, cont){
    proj_prop <- X_st[, feat_prop, drop = FALSE] %*% proj_dir_prop # Get proposed projection
    if(cont){
      max_knot0 <- quantile(proj_prop, p_dat_max)
      rg_knot0 <- (max_knot0 - min(proj_prop)) / prob_relu
      knot0_prop <- max_knot0 - rg_knot0 * runif(1)
      knots_prop <- c(knot0_prop, quantile(proj_prop[proj_prop > knot0_prop], knot_quants)) # Get proposed knots
      ridge_basis_prop <- get_mns_basis(proj_prop, knots_prop) # Get proposed basis functions
    }else{
      knots_prop <- NA
      ridge_basis_prop <- proj_prop
    }
    list(knots = knots_prop, basis = ridge_basis_prop)
  }

  birth <- function(){
    n_ridge_prop <- n_ridge + 1
    n_act_prop <- sample(n_act_max, 1, prob = w_n_act) # Propose number of active features
    if(adapt_act_feat){
      log_mh_act_feat <<- -(log(n_act_max) + log(w_n_act[n_act_prop]/sum(w_n_act))) # Nott, Kuk, Duc for n_act
      if(n_act_prop == 1){
        feat_prop <- feat_usable[sample(p_use, 1)] # Uniform over usable features, so no adjustment is needed
      }else{
        feat_prop <- sample(p, n_act_prop, prob = w_feat) # Propose features to include
        log_mh_act_feat <<- log_mh_act_feat - (lchoose(p_use, n_act_prop) + log(dwallenius(w_feat, feat_prop))) # Nott, Kuk, Duc for feat
      }
    }else{
      feat_prop <- sample(p, n_act_prop, prob = w_feat) # Propose features to include
    }

    if(all(feat_type[feat_prop] == 'cat')){ # Are all of the proposed features categorical?
      ridge_type_prop <- 'cat'
      n_quant_prop <- n_quant
      proj_dir_prop <- knots_prop <- NA
      ridge_basis_prop <- get_cat_basis(X_st[, feat_prop, drop = FALSE])
      n_basis_prop <- 1
    }else{
      n_quant_prop <- n_quant + 1
      if(n_act_prop == 1){
        proj_dir_prop <- matrix(sample(c(-1, 1), 1))
      }else{
        proj_dir_prop <- rps(proj_dir_mn[[n_act_prop]], 0) # Propose direction
      }
      cont <- any(feat_type[feat_prop] == 'cont') # Are any proposed features continuous?
      qb <- propose_quant_basis(feat_prop, proj_dir_prop, cont)
      knots_prop <- qb$knots
      ridge_basis_prop <- qb$basis
      if(cont){
        ridge_type_prop <- 'cont'
        n_basis_prop <- df_spline
      }else{ # The proposed features are a mix of categorical and discrete quantitative
        ridge_type_prop <- 'disc'
        n_basis_prop <- 1
      }
    }
    rb <- centre(ridge_basis_prop)
    ridge_basis_prop <- rb$B

    PtP <- t(ridge_basis_prop) %*% ridge_basis_prop # inner product of proposed new basis functions
    BtP <- t(basis_mat) %*% ridge_basis_prop
    Pty <- t(ridge_basis_prop) %*% y

    j_birth <- n_ridge + 1
    idx_K <- seq_len(n_basis_total)
    basis_idx_prop <- (n_basis_total + 1):(n_basis_total + n_basis_prop)
    BtB[idx_K, basis_idx_prop] <<- BtP
    BtB[basis_idx_prop, idx_K] <<- t(BtP)
    BtB[basis_idx_prop, basis_idx_prop] <<- PtP
    Bty[basis_idx_prop] <<- Pty

    idx_new <- seq_len(n_basis_total + n_basis_prop)
    qf_info_prop <- get_qf_info(BtB[idx_new, idx_new, drop = FALSE], Bty[idx_new], intercept)
    log_mh_bd_prop <- get_log_mh_bd(n_ridge_prop, n_quant_prop, n_ridge_max)

    if(!is.null(qf_info_prop)){
      if(qf_info_prop$qf < ssy){
        sse_prop <- ssy - c_var_coefs * qf_info_prop$qf

        # Compute the acceptance probability
        log_mh <- log_mh_bd - log_mh_bd_prop + log_mh_act_feat + # Adjustment for probability of birth proposal
          log_ml_ratio(sse_prop, sse) + # Part of the marginal likelihood
          log(n_ridge_mean/(n_ridge + 1)) # Prior and proposal distribution
        if(prior_coefs == 'zs'){
          log_mh <- log_mh - n_basis_prop/2 * log(var_coefs + 1) # The rest of the marginal likelihood for Zellner-Siow prior
        }else{
          log_mh <- log_mh + log(10e-6) # The rest of the marginal likelihood for flat prior
        }

        if(log(runif(1)) < log_mh){ # Accept the proposal
          n_ridge <<- j_birth
          n_act[j_birth] <<- n_act_prop
          feat[[j_birth]] <<- feat_prop
          knots[[j_birth]] <<- knots_prop
          proj_dir[[j_birth]] <<- proj_dir_prop
          if(!intercept) basis_mn[[j_birth]] <<- rb$mn

          if(ridge_type_prop != 'cat'){
            j_quant <<- c(j_quant, j_birth)
            n_quant <<- n_quant_prop
          }
          ridge_type[j_birth] <<- ridge_type_prop

          basis_idx[[j_birth + 1]] <<- basis_idx_prop
          n_basis_ridge[j_birth + 1] <<- n_basis_prop
          n_basis_total <<- n_basis_total + n_basis_prop
          basis_mat <<- cbind(basis_mat, ridge_basis_prop)

          # Update weights
          if(adapt_act_feat){
            w_n_act[n_act_prop] <<- w_n_act[n_act_prop] + 1
            w_feat[feat_prop] <<- w_feat[feat_prop] + 1
          }

          qf_info <<- qf_info_prop
          sse <<- sse_prop
          log_mh_bd <<- log_mh_bd_prop
          return(TRUE)
        }
      }
    }
    FALSE
  }

  death <- function(){
    j_death <- sample(n_ridge, 1) # Choose random index to delete

    n_act_prop <- n_act[j_death]
    feat_prop <- feat[[j_death]]
    if(adapt_act_feat){
      log_mh_act_feat <<- log(n_act_max) + log((w_n_act[n_act_prop] - 1)/(sum(w_n_act) - 1))
      w_feat_prop <- w_feat
      w_feat_prop[feat_prop] <- w_feat_prop[feat_prop] - 1
      if(n_act_prop > 1){
        log_mh_act_feat <<- log_mh_act_feat + lchoose(p_use, n_act_prop) + log(dwallenius(w_feat_prop, feat_prop)) # Nott, Kuk, and Duc
      }
    }

    basis_idx_vec_prop <- unlist(basis_idx)[-basis_idx[[j_death + 1]]]
    qf_info_prop <- get_qf_info(BtB[basis_idx_vec_prop, basis_idx_vec_prop, drop = FALSE], Bty[basis_idx_vec_prop], intercept)

    n_ridge_prop <- n_ridge - 1
    if(ridge_type[j_death] == 'cat'){
      n_quant_prop <- n_quant
    }else{
      n_quant_prop <- n_quant - 1
    }

    log_mh_bd_prop <- get_log_mh_bd(n_ridge_prop, n_quant_prop, n_ridge_max)

    if(!is.null(qf_info_prop)){
      if(qf_info_prop$qf < ssy){
        sse_prop <- ssy - c_var_coefs * qf_info_prop$qf
        n_basis_prop <- n_basis_ridge[j_death + 1]

        # Compute acceptance probability
        log_mh <- log_mh_bd - log_mh_bd_prop + log_mh_act_feat + # Adjustment for probability of death proposal
          log_ml_ratio(sse_prop, sse) + # Part of the marginal likelihood
          log(n_ridge/n_ridge_mean) # Prior and proposal distribution
        if(prior_coefs == 'zs'){
          log_mh <- log_mh + n_basis_prop/2 * log(var_coefs + 1) # The rest of the marginal likelihood for Zellner-Siow prior
        }else{
          log_mh <- log_mh - log(10e-6) # The rest of the marginal likelihood for flat prior
        }

        if(log(runif(1)) < log_mh){ # Accept the proposal
          if(ridge_type[j_death] != 'cat'){
            n_quant <<- n_quant - 1
            j_quant <<- j_quant[j_quant != j_death]
          }
          j_quant[j_quant > j_death] <<- j_quant[j_quant > j_death] - 1
          ridge_type <<- ridge_type[-j_death]

          idx_left <- seq_len(n_basis_total - n_basis_prop)
          basis_mat <<- basis_mat[, -basis_idx[[j_death + 1]], drop = FALSE]
          BtB[idx_left, idx_left] <<- BtB[basis_idx_vec_prop, basis_idx_vec_prop, drop = FALSE]
          Bty[idx_left] <<- Bty[basis_idx_vec_prop]
          if(j_death < n_ridge){
            for(j in (j_death + 1):n_ridge){
              basis_idx[[j + 1]] <<- basis_idx[[j + 1]] - n_basis_ridge[j_death + 1]
            }
          }
          basis_idx <<- basis_idx[-(j_death + 1)]

          n_basis_ridge <<- n_basis_ridge[-(j_death + 1)]
          n_basis_total <<- n_basis_total - n_basis_prop
          n_ridge <<- n_ridge - 1
          n_act <<- n_act[-j_death]
          feat <<- feat[-j_death]
          knots <<- knots[-j_death]
          proj_dir <<- proj_dir[-j_death]
          if(!intercept) basis_mn <<- basis_mn[-j_death]

          # Update weights
          if(adapt_act_feat){
            w_feat <<- w_feat_prop
            w_n_act[n_act_prop] <<- w_n_act[n_act_prop] - 1
          }

          qf_info <<- qf_info_prop
          sse <<- sse_prop
          log_mh_bd <<- log_mh_bd_prop
          return(TRUE)
        }
      }
    }
    FALSE
  }

  change <- function(){
    if(n_quant == 1){
      j_change <- j_quant
    }else{
      j_change <- sample(j_quant, 1) # Which ridge function should we change?
    }

    if(n_act[j_change] == 1){
      proj_dir_prop <- matrix(sample(c(-1, 1), 1))
    }else{
      proj_dir_prop <- rps(proj_dir[[j_change]], proj_dir_prop_prec) # Get proposed direction
    }

    qb <- propose_quant_basis(feat[[j_change]], proj_dir_prop, ridge_type[j_change] == 'cont')
    knots_prop <- qb$knots
    rb <- centre(qb$basis)
    ridge_basis_prop <- rb$B

    PtP <- t(ridge_basis_prop) %*% ridge_basis_prop # inner product of proposed new basis functions
    BtP <- t(basis_mat) %*% ridge_basis_prop
    Pty <- t(ridge_basis_prop) %*% y

    idx_K <- seq_len(n_basis_total)
    idx_j <- basis_idx[[j_change + 1]]
    BtB_prop <- BtB[idx_K, idx_K, drop = FALSE]
    BtB_prop[idx_j, ] <- t(BtP)
    BtB_prop[, idx_j] <- BtP
    BtB_prop[idx_j, idx_j] <- PtP

    Bty_prop <- Bty[idx_K]
    Bty_prop[idx_j] <- Pty

    qf_info_prop <- get_qf_info(BtB_prop, Bty_prop, intercept)

    if(!is.null(qf_info_prop)){
      if(qf_info_prop$qf < ssy){
        sse_prop <- ssy - c_var_coefs * qf_info_prop$qf

        # Compute the acceptance probability
        log_mh <- log_ml_ratio(sse_prop, sse) # Marginal Likelihood

        if(log(runif(1)) < log_mh){ # Accept the proposal
          knots[[j_change]] <<- knots_prop
          proj_dir[[j_change]] <<- proj_dir_prop
          if(!intercept) basis_mn[[j_change]] <<- rb$mn

          BtB[idx_K, idx_K] <<- BtB_prop
          Bty[idx_K] <<- Bty_prop
          basis_mat[, idx_j] <<- ridge_basis_prop

          qf_info <<- qf_info_prop
          sse <<- sse_prop
          return(TRUE)
        }
      }
    }
    FALSE
  }

  # Gibbs draws of the residual variance, coefficients and their prior variance
  sample_params <- function(){
    if(is.null(qf_info$inv_chol)){
      qf_info <<- append_qf_inv_chol(qf_info, dim = n_basis_total)
    }

    if(is.null(s2_fixed)){
      sd_resid <<- sqrt(1/rgamma(1, s2_shape + n/2, s2_rate + sse/2))
    }

    # Draw coefs
    if(n_basis_total > 0){
      coefs <<- c_var_coefs * qf_info$ls_est +
        sqrt(c_var_coefs) * sd_resid * qf_info$inv_chol %*% rnorm(n_basis_total) # Draw coefs
      temp <- qf_info$chol %*% coefs
      qf2 <- c(t(temp) %*% temp)
    }else{
      coefs <<- numeric(0)
      qf2 <- 0
    }

    if(prior_coefs == 'zs'){
      var_coefs <<- 1/rgamma(1,
                             shape_var_coefs + n_basis_total/2,
                             scale_var_coefs + qf2/(2*sd_resid^2))
      c_var_coefs <<- var_coefs / (var_coefs + 1)
      sse <<- ssy - c_var_coefs * qf_info$qf
    }
  }

  iterate <- function(sample_params = TRUE){
    move_type <- get_move_type(n_ridge, n_quant, n_ridge_max)
    acc <- switch(move_type, birth = birth(), death = death(), change())
    n_prop <<- n_prop + 1
    n_acc <<- n_acc + acc
    if(sample_params) self$sample_params()
    invisible(self)
  }

  step <- function(k = 1L, sample_params = TRUE){
    for(i in seq_len(k)) iterate(sample_params)
    invisible(self)
  }

  fitted <- function(){
    if(n_basis_total == 0) return(numeric(n))
    as.vector(basis_mat %*% coefs)
  }

  s2 <- function() if(is.null(s2_fixed)) sd_resid^2 else s2_fixed

  accept_rate <- function() n_acc / max(n_prop, 1)

  current <- function(){
    list(n_ridge = n_ridge, n_act = n_act, feat = feat, proj_dir = proj_dir, knots = knots,
         coefs = coefs, sd_resid = sd_resid, var_coefs = if(prior_coefs == 'zs') var_coefs,
         w_n_act = w_n_act, w_feat = w_feat)
  }

  save <- function(){
    k <- length(keep$n_ridge) + 1
    keep$n_ridge[k] <<- n_ridge
    keep$n_act[[k]] <<- n_act
    keep$feat[[k]] <<- feat
    keep$proj_dir[[k]] <<- proj_dir
    keep$knots[[k]] <<- knots
    keep$coefs[[k]] <<- coefs
    keep$sd_resid[k] <<- sd_resid
    if(prior_coefs == 'zs') keep$var_coefs[k] <<- var_coefs
    if(!intercept) keep$basis_mn[[k]] <<- basis_mn
    invisible(self)
  }

  n_saved <- function() length(keep$n_ridge)

  clear_saved <- function(){
    keep <<- lapply(keep, function(x) x[0])
    invisible(self)
  }

  as_bppr <- function(...){
    out <- list(n_keep = as.numeric(length(keep$n_ridge)), n_ridge = keep$n_ridge,
                n_act = keep$n_act, feat = keep$feat,
                proj_dir = keep$proj_dir, knots = keep$knots,
                coefs = keep$coefs, sd_resid = keep$sd_resid,
                var_coefs = if(prior_coefs == 'zs') keep$var_coefs else NA,
                mn_X = mn_X, sd_X = sd_X, X = X, y = y,
                n_ridge_mean = n_ridge_mean, n_ridge_max = n_ridge_max,
                n_act_max = n_act_max, df_spline = df_spline, prob_relu = prob_relu,
                prior_coefs = prior_coefs, shape_var_coefs = shape_var_coefs,
                scale_var_coefs = scale_var_coefs, n_dat_min = n_dat_min,
                scale_proj_dir_prop = scale_proj_dir_prop,
                adapt_act_feat = adapt_act_feat, w_n_act = w_n_act,
                w_feat = w_feat)
    out <- c(out, list(...))
    # Fields of the extended model, only when they differ from bppr()'s defaults
    if(!intercept){
      out$intercept <- FALSE
      out$basis_mn <- keep$basis_mn
    }
    if(!jeffreys){
      out$s2_shape <- s2_shape
      out$s2_rate <- s2_rate
    }
    if(!is.null(s2_fixed)) out$s2_fixed <- s2_fixed
    structure(out, class = 'bppr')
  }

  class(self) <- 'bppr_sampler'
  self
}

#' @export
print.bppr_sampler <- function(x, ...){
  cat(sprintf('BayesPPR sampler: n = %d, p = %d, %d ridge functions, %d stored draws, acceptance %.3f\n',
              x$n, x$p, as.integer(x$n_ridge), x$n_saved(), x$accept_rate()))
  invisible(x)
}
