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
#' @param w_n_act_init vector of initial weights for number of active variables in a ridge function, used in generating proposed basis functions. Default is \code{rep(1, n_act_max)}.
#' @param w_feat_init vector of initial weights for features to be used in generating proposed basis functions. Default is \code{rep(1, ncol(X))}.
#' @param n_post number of posterior draws to obtain from the Markov chain after burn-in.
#' @param n_burn number of draws to burn before obtaining \code{n_post} draws for inference.
#' @param n_thin keep every n_thin posterior draws after burn-in.
#' @param print_every print the iteration number every print_every iterations. Use \code{print_every = 0} to silence.
#' @param model "bppr" is the only valid option as of now.
#' @details Explores BayesPPR model space using RJMCMC. The BayesPPR model has \deqn{y = f(x) + \epsilon,  ~~\epsilon \sim N(0,\sigma^2)} \deqn{f(x) = \beta_0 + \sum_{j=1}^M \beta_j B_j(x)} and \eqn{B_j(x)} is a natural spline basis expansion. We use priors \deqn{\beta \sim N(0,\sigma^2/\tau (B'B)^{-1})} \deqn{M \sim Poisson(\lambda)} as well as the hyper-prior on the variance \eqn{\tau} of the coefficients \eqn{\beta} mentioned in the arguments above.
#' @return An object of class \code{"bppr"}. Predictions can be obtained by passing the entire object to the \code{predict.bppr} function.
#' @keywords nonparametric projection pursuit regression splines
#' @seealso \link{predict.bppr} for prediction.
#' @export
#' @import stats
#' @import utils
#' @example inst/examples.R
#'
bppr <- function(X, y, n_ridge_mean = 10, n_ridge_max = NULL, n_act_max = NULL, df_spline = 4, prob_relu = 2/3, prior_coefs = "zs", shape_var_coefs = NULL, rate_var_coefs = NULL, n_dat_min = NULL, scale_proj_dir_prop = NULL, w_n_act_init = NULL, w_feat_init = NULL, n_post = 1000, n_burn = 9000, n_thin = 1, print_every = 1000, model = 'bppr'){
  # Manage posterior draws
  if(n_thin > n_post){
    stop('n_thin > n_post. No posterior samples will be obtained.')
  }
  n_post <- n_post - n_post %% n_thin
  n_draws <- n_burn + n_post
  n_keep <- n_post/n_thin
  idx <- c(rep(1, n_burn), rep(1:n_keep, each = n_thin))

  # Pre-processing
  n <- length(y)
  p <- ncol(X)

  if(is.null(w_feat_init)){
    w_feat_init <- rep(1, p)
  }
  w_feat <- w_feat_init

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
        w_feat_init[j] <- 0
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

  if(is.null(scale_proj_dir_prop)){
    proj_dir_prop_prec <- 1000 # scale_proj_dir_prop = 0.002
  }else if(scale_proj_dir_prop > 1  ||  scale_proj_dir_prop <= 0){
    stop("scale_proj_dir_prop must be in (0, 1]")
  }else{
    proj_dir_prop_prec <- 1/scale_proj_dir_prop
    proj_dir_prop_prec <- (proj_dir_prop_prec - 1) + sqrt(proj_dir_prop_prec * (proj_dir_prop_prec - 1))
  }

  if(is.null(w_n_act_init)){
    w_n_act_init <- rep(1, n_act_max)
  }
  w_n_act <- w_n_act_init

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

  # Initialization
  if(prior_coefs == 'zs'){
    if(is.null(shape_var_coefs)){
      shape_var_coefs <- 0.5
    }
    if(is.null(rate_var_coefs)){
      rate_var_coefs <- n/2
    }
    var_coefs <- numeric(n_keep)
    var_coefs[1] <- 1/rgamma(1, shape_var_coefs, rate_var_coefs)
    c_var_coefs <- var_coefs[1] / (var_coefs[1] + 1)
  }else if(prior_coefs == 'flat'){
    c_var_coefs <- 1
  }else{
    stop("prior_coefs must be either 'zs' or 'flat'")
  }

  sd_resid <- numeric(n_keep) # Error standard deviation
  sd_resid[1] <- 1

  coefs <- lapply(1:n_keep, function(it) numeric(1)) # Basis coefficients
  coefs[[1]] <- mean(y)

  n_ridge <- numeric(n_keep) # Number of ridge functions
  n_act <- lapply(1:n_keep, function(it) numeric()) # Number of active features for jth ridge function
  feat <- lapply(1:n_keep, function(it) list()) # Features being used in jth ridge function
  knots <- lapply(1:n_keep, function(it) list()) # Location of knots for nsplines
  proj_dir <- lapply(1:n_keep, function(it) list()) # Ridge directions
  proj_dir_mn <- lapply(1:n_act_max, function(a) rep(1/sqrt(a), a)) # prior mean for proj_dir (arbitrary, since precision is zero)
  n_basis_ridge <- 1 # Number of basis functions in each ridge function
  n_basis_total <- sum(n_basis_ridge)

  basis_mat <- matrix(rep(1, n)) # Current basis matrix
  qf_info <- get_qf_info(basis_mat, y)
  qf_info <- append_qf_inv_chol(qf_info, dim = n_basis_total)
  basis_idx <- list(1) # Current indices of segments of basis functions

  ssy <- c(t(y) %*% y) # Keep track of overall sse
  sse <- ssy - c_var_coefs * qf_info$qf

  if(n_burn > 0){
    phase <- 'burn'
  }else{
    phase <- 'post-burn'
  }
  if(print_every > 0){
    start_time <- Sys.time()
    cat(paste0('MCMC iteration 1/', n_draws, ' (', phase, ') ',
               myTimestamp(), ' n ridge: 0', '\n'))
    silent <- FALSE
  }else{
    print_every <- n_draws + 2
    silent <- TRUE
  }


  # Run MCMC
  for(it in 2:n_draws){
    if(it == n_burn + 1) phase <- 'post-burn'
    if((it - 1) %% print_every == 0  ||  (it == n_burn + 1  &&  !silent)){
      pr <- paste0('MCMC iteration ', it, '/', n_draws, ' (', phase, ') ',
                   myTimestamp(start_time), ' n ridge: ', n_ridge[idx[it]])
      cat(pr, '\n')
    }

    # Set current it values to last it values (these will change during the iteration)
    if(prior_coefs == 'zs') var_coefs[idx[it]] <- var_coefs[idx[it - 1]]
    sd_resid[idx[it]] <- sd_resid[idx[it - 1]]
    coefs[[idx[it]]] <- coefs[[idx[it - 1]]]
    n_ridge[idx[it]] <- n_ridge[idx[it - 1]]
    n_act[[idx[it]]] <- n_act[[idx[it - 1]]]
    feat[[idx[it]]] <- feat[[idx[it - 1]]]
    knots[[idx[it]]] <- knots[[idx[it - 1]]]
    proj_dir[[idx[it]]] <- proj_dir[[idx[it - 1]]]

    # Perform Reversible Jump Step (birth, death, change)
    if(n_ridge[idx[it]] == 0){
      move_type <- 'birth'
      alpha0 <- log(1/3)
    }else if(n_ridge[idx[it]] == n_ridge_max){
      move_type <- sample(c('death', 'change'), 1)
      alpha0 <- log(2/3)
    }else{
      move_type <- sample(c('birth', 'death', 'change'), 1)
      if(n_ridge[idx[it]] == 1  &&  move_type == 'death'){
        alpha0 <- log(3)
      }else if(n_ridge[idx[it]] == (n_ridge_max - 1)  &&  move_type == 'birth'){
        alpha0 <- log(3/2)
      }else{
        alpha0 <- 0
      }
    }

    if(move_type == 'birth'){ # Birth step
      n_act_prop <- sample(n_act_max, 1, prob = w_n_act) # Propose number of active features
      alpha0 <- alpha0 - (log(n_act_max) + log(w_n_act[n_act_prop]/sum(w_n_act))) # Nott, Kuk, Duc for n_act
      if(n_act_prop == 1){
        feat_prop <- sample(p, 1)
      }else{
        feat_prop <- sample(p, n_act_prop, prob = w_feat) # Propose features to include
        alpha0 <- alpha0 - (lchoose(p, n_act_prop) + log(dwallenius(w_feat, feat_prop))) # Nott, Kuk, Duc for feat
      }

      if(all(feat_type[feat_prop] == 'cat')){ # Are all of the proposed features categorical?
        proj_dir_prop <- knots_prop <- NA
        ridge_basis_prop <- get_cat_basis(X_st[, feat_prop, drop = FALSE])
        n_basis_prop <- 1
      }else{
        if(n_act_prop == 1){
          proj_dir_prop <- matrix(sample(c(-1, 1), 1))
        }else{
          proj_dir_prop <- rps(proj_dir_mn[[n_act_prop]], 0) # Propose direction
        }
        proj_prop <- X_st[, feat_prop, drop = FALSE] %*% proj_dir_prop # Get proposed projection

        if(any(feat_type[feat_prop] == 'cont')){ # Are any proposed features continuous?
          max_knot0 <- quantile(proj_prop, p_dat_max)
          rg_knot0 <- (max_knot0 - min(proj_prop)) / prob_relu
          knot0_prop <- max_knot0 - rg_knot0 * runif(1)
          knots_prop <- c(knot0_prop, quantile(proj_prop[proj_prop > knot0_prop], knot_quants)) # Get proposed knots
          ridge_basis_prop <- get_mns_basis(proj_prop, knots_prop) # Get proposed basis functions
          n_basis_prop <- df_spline
        }else{ # The proposed features are a mix of categorical and discrete quantitative
          knots_prop <- NA
          ridge_basis_prop <- proj_prop
          n_basis_prop <- 1
        }
      }

      basis_mat_prop <- cbind(basis_mat, ridge_basis_prop)
      qf_info_prop <- get_qf_info(basis_mat_prop, y)

      if(!is.null(qf_info_prop)  &&  (sse_prop <- ssy - c_var_coefs * qf_info_prop$qf) > 0){
        # Compute the acceptance probability
        alpha <- alpha0 + # Adjustment for probability of birth proposal
          -n/2 * (log(sse_prop) - log(sse)) + # Part of the marginal likelihood
          log(n_ridge_mean/(n_ridge[idx[it]] + 1)) # Prior and proposal distribution
        if(prior_coefs == 'zs'){
          alpha <- alpha - n_basis_prop/2 * log(var_coefs[idx[it]] + 1) # The rest of the marginal likelihood for Zellner-Siow prior
        }else{
          alpha <- alpha + log(10e-6) # The rest of the marginal likelihood for flat prior
        }

        if(log(runif(1)) < alpha){ # Accept the proposal
          n_ridge[idx[it]] <- j_birth <- n_ridge[idx[it]] + 1
          n_act[[idx[it]]][j_birth] <- n_act_prop
          feat[[idx[it]]][[j_birth]] <- feat_prop
          knots[[idx[it]]][[j_birth]] <- knots_prop
          proj_dir[[idx[it]]][[j_birth]] <- proj_dir_prop

          basis_idx_start <- basis_idx[[j_birth]][n_basis_ridge[j_birth]]
          basis_idx[[j_birth + 1]] <- (basis_idx_start + 1):(basis_idx_start + n_basis_prop)
          n_basis_ridge[j_birth + 1] <- n_basis_prop
          n_basis_total <- n_basis_total + n_basis_prop
          basis_mat <- basis_mat_prop

          #Update weights
          w_n_act[n_act_prop] <- w_n_act[n_act_prop] + 1
          w_feat[feat_prop] <- w_feat[feat_prop] + 1

          qf_info <- qf_info_prop
          qf_info <- append_qf_inv_chol(qf_info, dim = n_basis_total)
          sse <- sse_prop
        }
      }
    }else if(move_type == 'death'){ # Death step
      j_death <- sample(n_ridge[idx[it]], 1) # Choose random index to delete

      n_act_prop <- n_act[[idx[it]]][j_death]
      alpha0 <- alpha0 + log(n_act_max) + log((w_n_act[n_act_prop] - 1)/(sum(w_n_act) - 1))

      feat_prop <- feat[[idx[it]]][[j_death]]
      w_feat_prop <- w_feat;
      w_feat_prop[feat_prop] <- w_feat_prop[feat_prop] - 1
      if(n_act_prop > 1){
        alpha0 <- alpha0 + lchoose(p, n_act_prop) + log(dwallenius(w_feat_prop, feat_prop)) # Nott, Kuk, and Duc
      }

      basis_mat_prop <- basis_mat[, -basis_idx[[j_death + 1]], drop = FALSE]
      qf_info_prop <- get_qf_info(basis_mat_prop, y)

      if(!is.null(qf_info_prop)  &&  (sse_prop <- ssy - c_var_coefs * qf_info_prop$qf) > 0){
        n_basis_prop <- n_basis_ridge[j_death + 1]

        # Compute acceptance probability
        alpha <- alpha0 + # Adjustment for probability of death proposal
          -n/2 * (log(sse_prop) - log(sse)) + # Part of the marginal likelihood
          log(n_ridge[idx[it]]/n_ridge_mean) # Prior and proposal distribution
        if(prior_coefs == 'zs'){
          alpha <- alpha + n_basis_prop/2 * log(var_coefs[idx[it]] + 1) # The rest of the marginal likelihood for Zellner-Siow prior
        }else{
          alpha <- alpha - log(10e-6) # The rest of the marginal likelihood for flat prior
        }

        if(log(runif(1)) < alpha){ # Accept the proposal
          if(j_death < n_ridge[idx[it]]){
            for(j in (j_death + 1):n_ridge[idx[it]]){
              basis_idx[[j + 1]] <- basis_idx[[j + 1]] - n_basis_ridge[j_death + 1]
            }
          }
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
          w_feat <- w_feat_prop
          w_n_act[n_act_prop] <- w_n_act[n_act_prop] - 1

          qf_info <- qf_info_prop
          qf_info <- append_qf_inv_chol(qf_info, dim = n_basis_total)
          sse <- sse_prop
        }
      }
    }else{ # Change Step
      j_change <- sample(n_ridge[idx[it]], 1) # Which ridge function should we change?
      if(!is.na(proj_dir[[idx[it]]][[j_change]][1])){ # Are any variables quantitative for this ridge function?
        if(n_act[[idx[it]]][j_change] == 1){
          proj_dir_prop <- matrix(sample(c(-1, 1), 1))
        }else{
          proj_dir_prop <- rps(proj_dir[[idx[it]]][[j_change]], proj_dir_prop_prec) # Get proposed direction
        }
        proj_prop <- X_st[, feat[[idx[it]]][[j_change]], drop = FALSE] %*% proj_dir_prop # Get proposed projection

        if(!is.na(knots[[idx[it]]][[j_change]][1])){ # Are any variables continuous for this ridge function?
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

        if(!is.null(qf_info_prop)  &&  (sse_prop <- ssy - c_var_coefs * qf_info_prop$qf) > 0){
          # Compute the acceptance probability
          alpha <- -n/2 * (log(sse_prop) - log(sse)) # Marginal Likelihood

          if(log(runif(1)) < alpha){ # Accept the proposal
            knots[[idx[it]]][[j_change]] <- knots_prop
            proj_dir[[idx[it]]][[j_change]] <- proj_dir_prop

            basis_mat <- basis_mat_prop

            qf_info <- qf_info_prop
            qf_info <- append_qf_inv_chol(qf_info, dim = n_basis_total)
            sse <- sse_prop
          }
        }
      }
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

  out <- list(n_keep = n_keep, n_ridge = n_ridge,
              n_act = n_act, feat = feat,
              proj_dir = proj_dir, knots = knots,
              coefs = coefs, sd_resid = sd_resid,
              w_n_act = w_n_act, w_feat = w_feat,
              mn_X = mn_X, sd_X = sd_X,
              df_spline = df_spline,
              X = X, y = y, call = match.call())
  if(prior_coefs == 'zs'){
    out$var_coefs <- var_coefs
  }

  if(!silent){
    cat(paste0('MCMC iteration ', n_draws, '/', n_draws, ' (', phase, ') ',
               myTimestamp(start_time), ' n ridge: ', n_ridge[idx[it]], '\n'))
  }

  structure(out, class = 'bppr')
}
