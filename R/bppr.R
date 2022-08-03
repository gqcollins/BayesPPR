########################################################################
## main BayesPPR function
########################################################################

#' @title Bayesian Projection Pursuit Regression
#'
#' @description Fits a BayesPPR model using RJMCMC. Can handle categorical features.
#' @param X a data frame or matrix of predictors. Categorical features should be coded as numeric.
#' @param y a numeric response vector.
#' @param n_ridge_mean mean for Poisson prior on the number of ridge functions.
#' @param n_ridge_max maximum number of ridge functions allowed in the model. Used to avoid memory overload. Defaults to 150 unless the number of observed responses is small.
#' @param n_act_max maximum number of active variables in any given ridge function. Defaults to 3 unless categorical features are detected, in which case the default is larger.
#' @param df_spline degrees of freedom for spline basis. Stability should be examined for anything other than 4.
#' @param prob_relu prior probability that any given ridge function uses a relu transformation.
#' @param prior_coefs form of the prior distribution for the basis coefficients. Default is \code{"zs"} for the Zellner-Siow prior. The other option is \code{"flat"}, which is an improper prior.
#' @param shape_var_coefs shape for IG prior on the variance of the basis function coefficients. Default is for the Zellner-Siow prior. For the flat, improper prior, \code{shape_var_coefs} is ignored.
#' @param rate_var_coefs rate for IG prior on the variance of the basis function coefficients. Default is for the Zellner-Siow prior. For the flat, improper prior, \code{rate_var_coefs} is ignored.
#' @param n_dat_min minimum number of observed non-zero datapoints in a ridge function. Defaults to 20 or 0.1 times the number of observations, whichever is smaller.
#' @param scale_proj_dir_prop scale parameter for generating proposed projection directions. Should be in (0, 1); default is about 0.002.
#' @param adapt_act_feat logical; if \code{TRUE}, use adaptive proposal for feature index sets and number of active features.
#' @param w_n_act vector of weights for number of active variables in a ridge function, used in generating proposed basis functions. If \code{adapt_act_feat == FALSE}, it is also used for the prior distribution. Default is \code{rep(1, n_act_max)}.
#' @param w_feat vector of weights for feature indices used in a ridge function, used in generating proposed basis functions. If \code{adapt_act_feat == FALSE}, it is also used for the prior distribution. Default is \code{rep(1, ncol(X))}.
#' @param n_post number of posterior draws to obtain from the Markov chain after burn-in.
#' @param n_burn number of draws to burn before obtaining \code{n_post} draws for inference. If \code{prior_coefs == "flat"} then these are absorbed into the adapt phase.
#' @param n_adapt number of adaptive MCMC iterations to perform before burn-in. Skips sampling basis coefficients and residual variance to save time.
#' @param n_thin keep every n_thin posterior draws after burn-in.
#' @param print_every print the iteration number every print_every iterations. Use \code{print_every = 0} to silence.
#' @param bppr_init list of initial values for the Markov chain. Used by \link{bppr_resume}.
#' @details Explores BayesPPR model space using RJMCMC. The BayesPPR model has \deqn{y = f(x) + \epsilon,  ~~\epsilon \sim N(0,\sigma^2)} \deqn{f(x) = \beta_0 + \sum_{j=1}^M \beta_j B_j(x)} and \eqn{B_j(x)} is a natural spline basis expansion. We use priors \deqn{\beta \sim N(0,\sigma^2/\tau (B'B)^{-1})} \deqn{M \sim Poisson(\lambda)} as well as the hyper-prior on the variance \eqn{\tau} of the coefficients \eqn{\beta} mentioned in the arguments above.
#' @return An object of class \code{"bppr"}. Predictions can be obtained by passing the entire object to the \code{predict.bppr} function.
#' @keywords nonparametric projection pursuit regression splines
#' @seealso \link{predict.bppr} for prediction.
#' @export
#' @import stats
#' @import utils
#' @example inst/examples.R
#'
bppr <- function(X, y, n_ridge_mean = 10, n_ridge_max = NULL, n_act_max = NULL, df_spline = 4, prob_relu = 2/3, prior_coefs = "zs", shape_var_coefs = NULL, rate_var_coefs = NULL, n_dat_min = NULL, scale_proj_dir_prop = NULL, adapt_act_feat = TRUE, w_n_act = NULL, w_feat = NULL, n_post = 1000, n_burn = 9000, n_adapt = 0, n_thin = 1, print_every = 1000, bppr_init = NULL){
  # Manage posterior draws
  if(n_thin > n_post){
    stop('n_thin > n_post. No posterior samples will be obtained.')
  }
  n_post <- n_post - n_post %% n_thin
  n_keep <- n_post/n_thin
  n_pre <- n_adapt + n_burn
  n_draws <- n_pre + n_post
  idx <- c(rep(1, n_pre), rep(1:n_keep, each = n_thin))
  if(prior_coefs == 'flat'){
    n_adapt <- n_pre
    n_burn <- 0
  }

  # Pre-processing
  n <- length(y)
  p <- ncol(X)

  if(is.null(w_feat)){
    w_feat <- rep(1, p)
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

  if(is.null(n_act_max)){
    n_cat <- sum(feat_type == 'cat')
    n_act_max <- min(3, p - n_cat) + min(3, ceiling(n_cat/2))
  }

  if(is.null(w_n_act)){
    w_n_act <- rep(1, n_act_max)
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
    n_ridge_max <- min(150, floor(length(y)/df_spline) - 2)
  }
  if(n_ridge_max <= 0){
    stop('n_ridge_max <= 0. If n_ridge_max was set by default, df_spline is too large compared to the sample size.')
  }

  knot_quants <- seq(0, 1, length.out = df_spline + 1) # Quantiles for knot locations (except initial knot)

  if(prior_coefs == 'zs'){
    if(is.null(shape_var_coefs)){
      shape_var_coefs <- 0.5
    }
    if(is.null(rate_var_coefs)){
      rate_var_coefs <- n/2
    }
    var_coefs <- numeric(n_keep)
  }else if(prior_coefs == 'flat'){
    var_coefs <- NULL
    c_var_coefs <- 1
  }else{
    stop("prior_coefs must be either 'zs' or 'flat'")
  }
  sd_resid <- numeric(n_keep) # Error standard deviation
  coefs <- vector('list', n_keep) # Basis coefficients
  n_ridge <- numeric(n_keep) # Number of ridge functions
  n_act <- vector('list', n_keep) # Number of active features for jth ridge function
  feat <- vector('list', n_keep) # Features being used in jth ridge function
  knots <- vector('list', n_keep) # Location of knots for nsplines
  proj_dir <- vector('list', n_keep) # Ridge directions
  proj_dir_mn <- lapply(1:n_act_max, function(a) rep(1/sqrt(a), a)) # prior mean for proj_dir (arbitrary, since precision is zero)

  n_basis_ridge <- 1 # Number of basis functions in each ridge function
  ridge_type <- character()
  j_quant <- integer()
  n_quant <- 0
  basis_mat <- matrix(rep(1, n)) # Current basis matrix
  basis_idx <- list(1) # Current indices of segments of basis functions

  # Initialization
  if(is.null(bppr_init)){
    if(prior_coefs == 'zs'){
      var_coefs[1] <- rate_var_coefs / shape_var_coefs
      c_var_coefs <- var_coefs[1] / (var_coefs[1] + 1)
    }
    sd_resid[1] <- 1
    coefs[[1]] <- mean(y)
    n_act[[1]] <- numeric() # Number of active features for jth ridge function
    feat[[1]] <- list() # Features being used in jth ridge function
    knots[[1]] <- list() # Location of knots for nsplines
    proj_dir[[1]] <- list() # Ridge directions
  }else{
    if(prior_coefs == 'zs'){
      var_coefs[1] <- bppr_init$var_coefs
      c_var_coefs <- var_coefs[1] / (var_coefs[1] + 1)
    }
    sd_resid[1] <- bppr_init$sd_resid
    coefs[[1]] <- bppr_init$coefs
    n_ridge[1] <- bppr_init$n_ridge
    n_act[[1]] <- bppr_init$n_act
    feat[[1]] <- bppr_init$feat
    knots[[1]] <- bppr_init$knots
    proj_dir[[1]] <- bppr_init$proj_dir
    if(n_ridge[1] > 0){
      basis_idx_start <- 2
      for(j in 1:n_ridge[1]){
        if(any(feat_type[feat[[1]][[j]]] == 'cont')){
          ridge_type[j] <- 'cont'
          n_basis_ridge[j + 1] <- df_spline
          j_quant <- c(j_quant, j)
          proj <- X_st[, feat[[1]][[j]], drop = FALSE] %*% proj_dir[[1]][[j]]
          basis_mat <- cbind(basis_mat, get_mns_basis(proj, knots[[1]][[j]])) # Get basis function
          basis_idx[[j + 1]] <- basis_idx_start:(basis_idx_start + df_spline - 1)
          basis_idx_start <- basis_idx_start + df_spline
        }else{
          n_basis_ridge[j + 1] <- 1
          if(any(feat_type[feat[[1]][[j]]] == 'disc')){
            ridge_type[j] <- 'disc'
            n_basis_ridge[j + 1] <- 1
            j_quant <- c(j_quant, j)
            basis_mat <- cbind(basis_mat, X_st[, feat[[1]][[j]], drop = FALSE] %*% proj_dir[[1]][[j]])
            basis_idx[[j + 1]] <- basis_idx_start
            basis_idx_start <- basis_idx_start + 1
          }else{
            ridge_type[j] <- 'cat'
            n_basis_ridge[j] <- 1
            basis_mat <- cbind(basis_mat, get_cat_basis(X_st[, feat[[1]][[j]], drop = FALSE]))
            basis_idx[[j + 1]] <- basis_idx_start
            basis_idx_start <- basis_idx_start + 1
          }
        }
      }
      n_quant <- length(j_quant)
    }
    w_n_act <- bppr_init$w_n_act
    w_feat <- bppr_init$w_feat
  }

  n_basis_total <- sum(n_basis_ridge)
  qf_info <- get_qf_info(basis_mat, y)
  if(n_adapt == 0) qf_info <- append_qf_inv_chol(qf_info, dim = n_basis_total)

  ssy <- c(t(y) %*% y) # Keep track of overall sse
  sse <- ssy - c_var_coefs * qf_info$qf
  log_mh_bd <- 0
  log_mh_act_feat <- 0

  if(n_adapt > 0){
    phase <- 'adapt'
  }else if(n_burn > 0){
    phase <- 'burn'
  }else{
    phase <- 'post-burn'
  }
  if(print_every > 0){
    start_time <- Sys.time()
    cat(paste0('MCMC iteration 1/', n_draws, ' (', phase, ') ',
               myTimestamp(), ' n ridge: ', n_ridge[1], '\n'))
    silent <- FALSE
  }else{
    print_every <- n_draws + 2
    silent <- TRUE
  }

  # Run MCMC
  if(n_draws > 1){
    for(it in 2:n_draws){
      # Set current it values to last it values (these will change during the iteration)
      if(idx[it] > idx[it - 1]){
        if(prior_coefs == 'zs') var_coefs[idx[it]] <- var_coefs[idx[it - 1]]
        sd_resid[idx[it]] <- sd_resid[idx[it - 1]]
        coefs[[idx[it]]] <- coefs[[idx[it - 1]]]
        n_ridge[idx[it]] <- n_ridge[idx[it - 1]]
        n_act[[idx[it]]] <- n_act[[idx[it - 1]]]
        feat[[idx[it]]] <- feat[[idx[it - 1]]]
        knots[[idx[it]]] <- knots[[idx[it - 1]]]
        proj_dir[[idx[it]]] <- proj_dir[[idx[it - 1]]]
      }

      if(it == n_adapt + 1){
        if(n_burn > 0){
          phase <- 'burn'
        }else{
          phase <- 'post-burn'
        }
      }
      if(it == n_pre + 1) phase <- 'post-burn'
      if((it - 1) %% print_every == 0  ||  ((it == n_adapt + 1  ||  it == n_pre + 1)  &&  !silent)){
        pr <- paste0('MCMC iteration ', it, '/', n_draws, ' (', phase, ') ',
                     myTimestamp(start_time), ' n ridge: ', n_ridge[idx[it]])
        cat(pr, '\n')
      }

      move_type <- get_move_type(n_ridge[idx[it]], n_quant, n_ridge_max)

      if(move_type == 'birth'){ # Birth step
        n_ridge_prop <- n_ridge[idx[it]] + 1
        n_act_prop <- sample(n_act_max, 1, prob = w_n_act) # Propose number of active features
        if(adapt_act_feat){
          log_mh_act_feat <- -(log(n_act_max) + log(w_n_act[n_act_prop]/sum(w_n_act))) # Nott, Kuk, Duc for n_act
          if(n_act_prop == 1){
            feat_prop <- sample(p, 1)
          }else{
            feat_prop <- sample(p, n_act_prop, prob = w_feat) # Propose features to include
            log_mh_act_feat <- log_mh_act_feat - (lchoose(p, n_act_prop) + log(dwallenius(w_feat, feat_prop))) # Nott, Kuk, Duc for feat
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
          proj_prop <- X_st[, feat_prop, drop = FALSE] %*% proj_dir_prop # Get proposed projection

          if(any(feat_type[feat_prop] == 'cont')){ # Are any proposed features continuous?
            ridge_type_prop <- 'cont'
            max_knot0 <- quantile(proj_prop, p_dat_max)
            rg_knot0 <- (max_knot0 - min(proj_prop)) / prob_relu
            knot0_prop <- max_knot0 - rg_knot0 * runif(1)
            knots_prop <- c(knot0_prop, quantile(proj_prop[proj_prop > knot0_prop], knot_quants)) # Get proposed knots
            ridge_basis_prop <- get_mns_basis(proj_prop, knots_prop) # Get proposed basis functions
            n_basis_prop <- df_spline
          }else{ # The proposed features are a mix of categorical and discrete quantitative
            ridge_type_prop <- 'disc'
            knots_prop <- NA
            ridge_basis_prop <- proj_prop
            n_basis_prop <- 1
          }
        }

        basis_mat_prop <- cbind(basis_mat, ridge_basis_prop)
        qf_info_prop <- get_qf_info(basis_mat_prop, y)
        log_mh_bd_prop <- get_log_mh_bd(n_ridge_prop, n_quant_prop, n_ridge_max)

        if(!is.null(qf_info_prop)){
          if(qf_info_prop$qf < ssy){
            sse_prop <- ssy - c_var_coefs * qf_info_prop$qf

            # Compute the acceptance probability
            log_mh <- log_mh_bd - log_mh_bd_prop + log_mh_act_feat + # Adjustment for probability of birth proposal
              -n/2 * (log(sse_prop) - log(sse)) + # Part of the marginal likelihood
              log(n_ridge_mean/(n_ridge[idx[it]] + 1)) # Prior and proposal distribution
            if(prior_coefs == 'zs'){
              log_mh <- log_mh - n_basis_prop/2 * log(var_coefs[idx[it]] + 1) # The rest of the marginal likelihood for Zellner-Siow prior
            }else{
              log_mh <- log_mh + log(10e-6) # The rest of the marginal likelihood for flat prior
            }

            if(log(runif(1)) < log_mh){ # Accept the proposal
              n_ridge[idx[it]] <- j_birth <- n_ridge[idx[it]] + 1
              n_act[[idx[it]]][j_birth] <- n_act_prop
              feat[[idx[it]]][[j_birth]] <- feat_prop
              knots[[idx[it]]][[j_birth]] <- knots_prop
              proj_dir[[idx[it]]][[j_birth]] <- proj_dir_prop

              if(ridge_type_prop != 'cat'){
                j_quant <- c(j_quant, j_birth)
                n_quant <- n_quant_prop
              }
              ridge_type[j_birth] <- ridge_type_prop

              basis_idx_start <- basis_idx[[j_birth]][n_basis_ridge[j_birth]]
              basis_idx[[j_birth + 1]] <- (basis_idx_start + 1):(basis_idx_start + n_basis_prop)
              n_basis_ridge[j_birth + 1] <- n_basis_prop
              n_basis_total <- n_basis_total + n_basis_prop
              basis_mat <- basis_mat_prop

              # Update weights
              if(adapt_act_feat){
                w_n_act[n_act_prop] <- w_n_act[n_act_prop] + 1
                w_feat[feat_prop] <- w_feat[feat_prop] + 1
              }

              qf_info <- qf_info_prop
              if(it > n_adapt) qf_info <- append_qf_inv_chol(qf_info, dim = n_basis_total)
              sse <- sse_prop
              log_mh_bd <- log_mh_bd_prop
            }
          }
        }
      }else if(move_type == 'death'){ # Death step
        j_death <- sample(n_ridge[idx[it]], 1) # Choose random index to delete

        n_act_prop <- n_act[[idx[it]]][j_death]
        feat_prop <- feat[[idx[it]]][[j_death]]
        if(adapt_act_feat){
          log_mh_act_feat <- log(n_act_max) + log((w_n_act[n_act_prop] - 1)/(sum(w_n_act) - 1))
          w_feat_prop <- w_feat
          w_feat_prop[feat_prop] <- w_feat_prop[feat_prop] - 1
          if(n_act_prop > 1){
            log_mh_act_feat <- log_mh_act_feat + lchoose(p, n_act_prop) + log(dwallenius(w_feat_prop, feat_prop)) # Nott, Kuk, and Duc
          }
        }

        basis_mat_prop <- basis_mat[, -basis_idx[[j_death + 1]], drop = FALSE]
        qf_info_prop <- get_qf_info(basis_mat_prop, y)
        n_ridge_prop <- n_ridge[idx[it]] - 1
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
              -n/2 * (log(sse_prop) - log(sse)) + # Part of the marginal likelihood
              log(n_ridge[idx[it]]/n_ridge_mean) # Prior and proposal distribution
            if(prior_coefs == 'zs'){
              log_mh <- log_mh + n_basis_prop/2 * log(var_coefs[idx[it]] + 1) # The rest of the marginal likelihood for Zellner-Siow prior
            }else{
              log_mh <- log_mh - log(10e-6) # The rest of the marginal likelihood for flat prior
            }

            if(log(runif(1)) < log_mh){ # Accept the proposal
              if(j_death < n_ridge[idx[it]]){
                for(j in (j_death + 1):n_ridge[idx[it]]){
                  basis_idx[[j + 1]] <- basis_idx[[j + 1]] - n_basis_ridge[j_death + 1]
                }
              }

              if(ridge_type[j_death] != 'cat'){
                n_quant <- n_quant - 1
                j_quant <- j_quant[j_quant != j_death]
              }
              j_quant[j_quant > j_death] <- j_quant[j_quant > j_death] - 1
              ridge_type <- ridge_type[-j_death]

              basis_idx <- basis_idx[-(j_death + 1)]
              n_basis_ridge <- n_basis_ridge[-(j_death + 1)]
              n_basis_total <- n_basis_total - n_basis_prop
              n_ridge[idx[it]] <- n_ridge[idx[it]] - 1
              n_act[[idx[it]]] <- n_act[[idx[it]]][-j_death]
              feat[[idx[it]]] <- feat[[idx[it]]][-j_death]
              knots[[idx[it]]] <- knots[[idx[it]]][-j_death]
              proj_dir[[idx[it]]] <- proj_dir[[idx[it]]][-j_death]

              basis_mat <- basis_mat_prop

              # Update weights
              if(adapt_act_feat){
                w_feat <- w_feat_prop
                w_n_act[n_act_prop] <- w_n_act[n_act_prop] - 1
              }

              qf_info <- qf_info_prop
              if(it > n_adapt) qf_info <- append_qf_inv_chol(qf_info, dim = n_basis_total)
              sse <- sse_prop
              log_mh_bd <- log_mh_bd_prop
            }
          }
        }
      }else{ # Change Step
        if(n_quant  == 1){
          j_change <- j_quant
        }else{
          j_change <- sample(j_quant, 1) # Which ridge function should we change?
        }

        if(n_act[[idx[it]]][j_change] == 1){
          proj_dir_prop <- matrix(sample(c(-1, 1), 1))
        }else{
          proj_dir_prop <- rps(proj_dir[[idx[it]]][[j_change]], proj_dir_prop_prec) # Get proposed direction
        }

        proj_prop <- X_st[, feat[[idx[it]]][[j_change]], drop = FALSE] %*% proj_dir_prop # Get proposed projection

        if(ridge_type[j_change] == 'cont'){ # Are any variables continuous for this ridge function?
          max_knot0 <- quantile(proj_prop, p_dat_max)
          rg_knot0 <- (max_knot0 -  min(proj_prop)) / prob_relu
          knot0_prop <- max_knot0 - rg_knot0 * runif(1)
          knots_prop <- c(knot0_prop, quantile(proj_prop[proj_prop > knot0_prop], knot_quants)) # Get proposed knots
          ridge_basis_prop <- get_mns_basis(proj_prop, knots_prop) # Get proposed basis function
        }else{
          knots_prop <- NA
          ridge_basis_prop <- proj_prop
        }

        basis_mat_prop <- basis_mat
        basis_mat_prop[, basis_idx[[j_change + 1]]] <- ridge_basis_prop
        qf_info_prop <- get_qf_info(basis_mat_prop, y)

        if(!is.null(qf_info_prop)){
          if(qf_info_prop$qf < ssy){
            sse_prop <- ssy - c_var_coefs * qf_info_prop$qf

            # Compute the acceptance probability
            log_mh <- -n/2 * (log(sse_prop) - log(sse)) # Marginal Likelihood

            if(log(runif(1)) < log_mh){ # Accept the proposal
              knots[[idx[it]]][[j_change]] <- knots_prop
              proj_dir[[idx[it]]][[j_change]] <- proj_dir_prop

              basis_mat <- basis_mat_prop

              qf_info <- qf_info_prop
              if(it > n_adapt) qf_info <- append_qf_inv_chol(qf_info, dim = n_basis_total)
              sse <- sse_prop
            }
          }
        }
      }

      if(it > n_adapt){
        if(it == n_adapt + 1  &&  is.null(qf_info$inv_chol)){
          qf_info <- append_qf_inv_chol(qf_info, dim = n_basis_total)
        }

        # Draw coefs
        coefs[[idx[it]]] <- c_var_coefs * qf_info$ls_est +
          sqrt(c_var_coefs) * sd_resid[idx[it]] * qf_info$inv_chol %*% rnorm(n_basis_total) # Draw coefs
        preds <- basis_mat %*% coefs[[idx[it]]] # Current predictions of y
        resid <- y - preds # current residuals

        sd_resid[idx[it]] <- sqrt(1/rgamma(1, n/2, c(t(resid) %*% resid)/2))

        if(prior_coefs == 'zs'){
          var_coefs[idx[it]] <- 1/rgamma(1,
                                         shape_var_coefs + n_basis_total/2,
                                         rate_var_coefs + c(t(preds) %*% preds)/(2*sd_resid[idx[it]]^2))
          c_var_coefs <- var_coefs[idx[it]] / (var_coefs[idx[it]] + 1)
          sse <- ssy - c_var_coefs * qf_info$qf
        }
      }
    }
  }

  out <- list(n_keep = n_keep, n_ridge = n_ridge,
              n_act = n_act, feat = feat,
              proj_dir = proj_dir, knots = knots,
              coefs = coefs, sd_resid = sd_resid, var_coefs = var_coefs,
              mn_X = mn_X, sd_X = sd_X, X = X, y = y,
              n_ridge_mean = n_ridge_mean, n_ridge_max = n_ridge_max,
              n_act_max = n_act_max, df_spline = df_spline, prob_relu = prob_relu,
              prior_coefs = prior_coefs, shape_var_coefs = shape_var_coefs,
              rate_var_coefs = rate_var_coefs, n_dat_min = n_dat_min,
              scale_proj_dir_prop = scale_proj_dir_prop,
              adapt_act_feat = adapt_act_feat, w_n_act = w_n_act,
              w_feat = w_feat, n_post = n_post, n_burn = n_burn,
              n_adapt = n_adapt, n_thin = n_thin, print_every = print_every,
              bppr_init = bppr_init, call = match.call())

  if(!silent){
    cat(paste0('MCMC iteration ', n_draws, '/', n_draws, ' (', phase, ') ',
               myTimestamp(start_time), ' n ridge: ', n_ridge[idx[n_draws]], '\n'))
  }

  return(structure(out, class = 'bppr'))
}
