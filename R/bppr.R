########################################################################
## main BayesPPR function
########################################################################

#' @title Bayesian Projection Pursuit Regression (BayesPPR)
#'
#' @description Fits a BayesPPR model using RJMCMC. Can handle categorical features.
#' @param X a data frame or matrix of predictors. Categorical features should be coded as numeric.
#' @param y a numeric response vector.
#' @param n_ridge_mean mean for Poisson prior on the number of ridge functions.
#' @param n_ridge_max maximum number of ridge functions allowed in the model. Used to avoid memory overload. Defaults to 150 unless the number of observed responses is small.
#' @param n_act_max maximum number of active variables in any given ridge function. Defaults to 3 unless categorical features are detected, in which case the default is larger.
#' @param df_spline degrees of freedom for spline basis. Stability should be examined for anything other than 4.
#' @param prob_relu prior probability that any given ridge function uses a relu transformation.
#' @param var_coefs_shape shape for IG prior on the variance of the basis function coefficients. Default is for the Zellner-Siow prior.
#' @param var_coefs_rate rate for IG prior on the variance of the basis function coefficients. Default is for the Zellner-Siow prior.
#' @param n_dat_act_min minimum number of observed non-zero values in a ridge function. Defaults to 20 or 0.1 times the number of observations, whichever is smaller.
#' @param proj_dir_prop_scale scale parameter for generating proposed projection directions. Should be in (0, 1); default is about 0.002.
#' @param n_act_w_init vector of initial weights for number of active variables in a ridge function, used in generating proposed basis functions. Default is \code{rep(1, n_act_max)}.
#' @param feat_w_init vector of initial weights for features to be used in generating proposed basis functions. Default is \code{rep(1, ncol(X))}.
#' @param n_draws number of draws to obtain from the Markov chain.
#' @param model "bppr" is the only valid option as of now.
#' @details Explores BayesPPR model space using RJMCMC. The BayesPPR model has \deqn{y = f(x) + \epsilon,  ~~\epsilon \sim N(0,\sigma^2)} \deqn{f(x) = \beta_0 + \sum_{j=1}^M \beta_j B_j(x)} and \eqn{B_j(x)} is a natural spline basis expansion. We use priors \deqn{\beta \sim N(0,\sigma^2/\tau (B'B)^{-1})} \deqn{M \sim Poisson(\lambda)} as well as the hyper-prior on the variance \eqn{\tau} of the coefficients \eqn{\beta} mentioned in the arguments above.
#' @return An object of class 'bppr'. Predictions can be obtained by passing the entire object to the predict.bppr function.
#' @keywords nonparametric projection pursuit regression splines
#' @seealso \link{predict.bppr} for prediction.
#' @export
#' @import stats
#' @import utils
#' @example inst/examples.R
#'
bppr <- function(X, y, n_ridge_mean = 10, n_ridge_max = NULL, n_act_max = NULL, df_spline = 4, prob_relu = 2/3, var_coefs_shape = 0.5, var_coefs_rate = length(y)/2, n_dat_act_min = NULL, proj_dir_prop_scale = NULL, n_act_w_init = NULL, feat_w_init = NULL, n_draws = 10000, model = 'bppr'){
  # Pre-processing
  n <- length(y)
  p <- ncol(X)

  if(is.null(feat_w_init)){
    feat_w_init <- rep(1, p)
  }
  feat_w <- feat_w_init

  mn_X <- sd_X <- numeric(p)
  feat_type <- character(p)
  for(j in 1:p){
    n_unique <- length(unique(X[, j]))
    if(n_unique <= 2){
      mn_X[j] <- 0
      sd_X[j] <- 1
      if(n_unique == 1){
        feat_type[j] <- ''
        feat_w_init[j] <- 0
      }else{
        feat_type[j] <- 'cat'
      }
    }else{
      mn_X[j] <- mean(X[, j])
      sd_X[j] <- sd(X[, j])
      X[, j] <- (X[, j] - mn_X[j]) / sd_X[j]
      if(n_unique <= df_spline){
        feat_type[j] <- 'disc'
      }else{
        feat_type[j] <- 'cont'
      }
    }
  }

  if(is.null(n_act_max)){
    p_cat <- sum(feat_type == 'cat')
    n_act_max <- min(3, p - p_cat) + min(3, ceiling(p_cat/2))
  }

  if(is.null(proj_dir_prop_scale)){
    proj_dir_prop_prec <- 1000 # proj_dir_prop_scale = 0.002
  }else if(proj_dir_prop_scale > 1  ||  proj_dir_prop_scale <= 0){
    stop("scale prop must be in (0, 1]")
  }else{
    proj_dir_prop_prec <- 1/proj_dir_prop_scale
    proj_dir_prop_prec <- (proj_dir_prop_prec - 1) + sqrt(proj_dir_prop_prec * (proj_dir_prop_prec - 1))
  }

  if(is.null(n_act_w_init)){
    n_act_w_init <- rep(1, n_act_max)
  }
  n_act_w <- n_act_w_init

  if(is.null(n_dat_act_min)){
    n_dat_act_min <- min(20, 0.1 * n)
  }
  p_dat_inact_max <- 1 - n_dat_act_min / n # Maximum proportion of inactive datapoints in each ridge runction

  if(is.null(n_ridge_max)){
    n_ridge_max <- min(150, floor(length(y)/df_spline) - 1)
  }

  knot_quants <- seq(0, 1, length.out = df_spline + 1) # Quantiles for knot locations

  prob_no_relu <- 1 - prob_relu

  # Initialization
  sd_resid <- numeric(n_draws) # Error standard deviation
  sd_resid[1] <- 1

  coefs <- lapply(1:n_draws, function(it) numeric(1)) # Ridge coefficients
  coefs[[1]] <- mean(y)

  var_coefs <- numeric(n_draws)
  var_coefs[1] <- 1/rgamma(1, var_coefs_shape, var_coefs_rate)

  n_ridge <- numeric(n_draws) # Number of ridge functions
  n_act <- lapply(1:n_draws, function(it) numeric()) # Number of active features for jth ridge function
  feat <- lapply(1:n_draws, function(it) list()) # Features being used in jth ridge function
  knots <- lapply(1:n_draws, function(it) list()) # Location of knots for nsplines
  bias <- lapply(1:n_draws, function(it) numeric()) # Bias term for each ridge function
  proj_dir <- lapply(1:n_draws, function(it) list()) # Ridge directions
  proj_dir_mn <- lapply(1:n_act_max, function(a) rep(1/sqrt(a), a)) # prior mean for proj_dir (arbitrary, since precision is zero)
  n_basis_ridge <- 1 # Number of basis functions in each ridge function
  n_basis_total <- sum(n_basis_ridge)

  basis_mat <- matrix(rep(1, n)) # Current basis matrix
  qf_info <- get_qf_info(basis_mat, y)
  qf_info <- append_qf_inv_chol(qf_info, dim = n_basis_total)
  basis_idx <- list(1) # Current indices of segments of basis functions

  ssy <- c(t(y) %*% y) # Keep track of overall sse
  sse <- ssy - var_coefs[1]/(var_coefs[1] + 1) * qf_info$qf

  # Run MCMC
  for(it in 2:n_draws){
    # Set current it values to last it values (these will change during the iteration)
    sd_resid[it] <- sd_resid[it - 1]
    coefs[[it]] <- coefs[[it - 1]]
    var_coefs[it] <- var_coefs[it - 1]
    n_ridge[it] <- n_ridge[it - 1]
    n_act[[it]] <- n_act[[it - 1]]
    feat[[it]] <- feat[[it - 1]]
    knots[[it]] <- knots[[it - 1]]
    bias[[it]] <- bias[[it - 1]]
    proj_dir[[it]] <- proj_dir[[it - 1]]

    qf_info <- get_qf_info(basis_mat, y)
    qf_info <- append_qf_inv_chol(qf_info, dim = n_basis_total)

    # Perform Reversible Jump Step (birth, death, change)
    if(n_ridge[it] == 0){
      move_type <- 'birth'
      alpha0 <- log(1/3)
    }else if(n_ridge[it] == n_ridge_max){
      move_type <- sample(c('death', 'change'), 1)
      alpha0 <- log(2/3)
    }else{
      move_type <- sample(c('birth', 'death', 'change'), 1)
      if(n_ridge[it] == 1  &&  move_type == 'death'){
        alpha0 <- log(3)
      }else if(n_ridge[it] == (n_ridge_max - 1)  &&  move_type == 'birth'){
        alpha0 <- log(3/2)
      }else{
        alpha0 <- 0
      }
    }

    if(move_type == 'birth'){ # Birth step
      n_act_prop <- sample(n_act_max, 1, prob = n_act_w) # Propose number of active features
      alpha0 <- alpha0 - (log(n_act_max) + log(n_act_w[n_act_prop]/sum(n_act_w))) # Nott, Kuk, Duc for n_act
      if(n_act_prop == 1){
        feat_prop <- sample(p, 1)
      }else{
        feat_prop <- sample(p, n_act_prop, prob = feat_w) # Propose features to include
        alpha0 <- alpha0 - (lchoose(p, n_act_prop) + log(dwallenius(feat_w, feat_prop))) # Nott, Kuk, Duc for feat
      }

      if(all(feat_type[feat_prop] == 'cat')){ # Are all of the proposed features categorical?
        proj_dir_prop <- bias_prop <- knots_prop <- NA
        ridge_basis_prop <- get_cat_basis(X[, feat_prop, drop = FALSE])
        n_basis_prop <- 1
      }else{
        if(n_act_prop == 1){
          proj_dir_prop <- matrix(sample(c(-1, 1), 1))
        }else{
          proj_dir_prop <- rps(proj_dir_mn[[n_act_prop]], 0) # Propose direction
        }
        proj_prop <- X[, feat_prop, drop = FALSE] %*% proj_dir_prop # Get proposed projection

        if(any(feat_type[feat_prop] == 'cont')){ # Are any proposed features continuous?
          if(runif(1) < prob_no_relu){ # Try not using relu with some probability
            bias_prop <- NA
            knots_prop <- quantile(proj_prop, knot_quants) # Get proposed knots
            ridge_basis_prop <- get_ns_basis(proj_prop, knots_prop) # Get proposed basis functions
          }else{
            min_X_proj_dir <- min(proj_prop)
            rg <- quantile(proj_prop, p_dat_inact_max) - min_X_proj_dir
            bias_prop <- -(rg * runif(1) + min_X_proj_dir)
            proj_prop_trans <- relu(bias_prop + proj_prop) # Get transformation of projection
            knots_prop <- quantile(proj_prop_trans[proj_prop_trans > 0], knot_quants) # Get proposed knots
            ridge_basis_prop <- get_ns_basis(proj_prop_trans, knots_prop) # Get proposed basis functions
          }
          n_basis_prop <- df_spline
        }else{ # The proposed features are a mix of categorical and discrete quantitative
          bias_prop <- knots_prop <- NA
          ridge_basis_prop <- proj_prop
          n_basis_prop <- 1
        }
      }

      basis_mat_prop <- cbind(basis_mat, ridge_basis_prop)
      qf_info_prop <- get_qf_info(basis_mat_prop, y)

      if(!is.null(qf_info_prop)  &&  (sse_prop <- ssy - var_coefs[it]/(var_coefs[it] + 1) * qf_info_prop$qf) > 0){
        # Compute the acceptance probability
        alpha <- alpha0 + # Adjustment for probability of birth proposal
          -n/2 * (log(sse_prop) - log(sse)) - n_basis_prop/2 * log(var_coefs[it] + 1) + # Marginal likelihood
          log(n_ridge_mean/(n_ridge[it] + 1)) # Prior and proposal distribution

        if(log(runif(1)) < alpha){ # Accept the proposal
          n_ridge[it] <- j_birth <- n_ridge[it] + 1
          n_act[[it]][j_birth] <- n_act_prop
          feat[[it]][[j_birth]] <- feat_prop
          knots[[it]][[j_birth]] <- knots_prop
          bias[[it]][j_birth] <- bias_prop
          proj_dir[[it]][[j_birth]] <- proj_dir_prop

          basis_idx_start <- basis_idx[[j_birth]][n_basis_ridge[j_birth]]
          basis_idx[[j_birth + 1]] <- (basis_idx_start + 1):(basis_idx_start + n_basis_prop)
          n_basis_ridge[j_birth + 1] <- n_basis_prop
          n_basis_total <- n_basis_total + n_basis_prop
          basis_mat <- basis_mat_prop

          #Update weights
          n_act_w[n_act_prop] <- n_act_w[n_act_prop] + 1
          feat_w[feat_prop] <- feat_w[feat_prop] + 1

          qf_info <- qf_info_prop
          qf_info <- append_qf_inv_chol(qf_info, dim = n_basis_total)
          sse <- sse_prop
        }
      }
    }else if(move_type == 'death'){ # Death step
      j_death <- sample(n_ridge[it], 1) # Choose random index to delete

      n_act_prop <- n_act[[it]][j_death]
      alpha0 <- alpha0 + log(n_act_max) + log((n_act_w[n_act_prop] - 1)/(sum(n_act_w) - 1))

      feat_prop <- feat[[it]][[j_death]]
      feat_w_prop <- feat_w;
      feat_w_prop[feat_prop] <- feat_w_prop[feat_prop] - 1
      if(n_act_prop > 1){
        alpha0 <- alpha0 + lchoose(p, n_act_prop) + log(dwallenius(feat_w_prop, feat_prop)) # Nott, Kuk, and Duc
      }

      basis_mat_prop <- basis_mat[, -basis_idx[[j_death + 1]], drop = FALSE]
      qf_info_prop <- get_qf_info(basis_mat_prop, y)

      if(!is.null(qf_info_prop)  &&  (sse_prop <- ssy - var_coefs[it]/(var_coefs[it] + 1) * qf_info_prop$qf) > 0){
        n_basis_prop <- n_basis_ridge[j_death + 1]

        # Compute acceptance probability
        alpha <- alpha0 + # Adjustment for probability of death proposal
          -n/2 * (log(sse_prop) - log(sse)) + n_basis_prop/2 * log(var_coefs[it] + 1) + # Marginal likelihood
          log(n_ridge[it]/n_ridge_mean) # Prior and proposal distribution

        if(log(runif(1)) < alpha){ # Accept the proposal
          if(j_death < n_ridge[it]){
            for(j in (j_death + 1):n_ridge[it]){
              basis_idx[[j + 1]] <- basis_idx[[j + 1]] - n_basis_ridge[j_death + 1]
            }
          }
          basis_idx <- basis_idx[-(j_death + 1)]
          n_basis_ridge <- n_basis_ridge[-(j_death + 1)]
          n_basis_total <- n_basis_total - n_basis_prop
          n_ridge[it] <- n_ridge[it] - 1
          n_act[[it]] <- n_act[[it]][-j_death]
          feat[[it]] <- feat[[it]][-j_death]
          knots[[it]] <- knots[[it]][-j_death]
          bias[[it]] <- bias[[it]][-j_death]
          proj_dir[[it]] <- proj_dir[[it]][-j_death]

          basis_mat <- basis_mat_prop

          # Update weights
          feat_w <- feat_w_prop
          n_act_w[n_act_prop] <- n_act_w[n_act_prop] - 1

          qf_info <- qf_info_prop
          qf_info <- append_qf_inv_chol(qf_info, dim = n_basis_total)
          sse <- sse_prop
        }
      }
    }else{ # Change Step
      j_change <- sample(n_ridge[it], 1) # Which ridge function should we change?
      if(!is.na(proj_dir[[it]][[j_change]][1])){ # Are any variables quantitative for this ridge function?
        if(n_act[[it]][j_change] == 1){
          proj_dir_prop <- matrix(sample(c(-1, 1), 1))
        }else{
          proj_dir_prop <- rps(proj_dir[[it]][[j_change]], proj_dir_prop_prec) # Get proposed direction
        }
        proj_prop <- X[, feat[[it]][[j_change]], drop = FALSE] %*% proj_dir_prop # Get proposed projection

        if(!is.na(knots[[it]][[j_change]][1])){ # Are any variables continuous for this ridge function?
          if(runif(1) < prob_no_relu){
            bias_prop <- NA
            knots_prop <- quantile(proj_prop, knot_quants) # Get proposed knots
            ridge_basis_prop <- get_ns_basis(proj_prop, knots_prop) # Get proposed basis function
          }else{
            min_X_proj_dir <- min(proj_prop)
            rg <- quantile(proj_prop, p_dat_inact_max) - min_X_proj_dir
            bias_prop <- -(rg * runif(1) + min_X_proj_dir)
            proj_prop_trans <- relu(bias_prop + proj_prop) # Get transformation of projection
            knots_prop <- quantile(proj_prop_trans[proj_prop_trans > 0], knot_quants) # Get proposed knots
            ridge_basis_prop <- get_ns_basis(proj_prop_trans, knots_prop) # Get proposed basis function
          }
        }else{
          bias_prop <- knots_prop <- NA
          ridge_basis_prop <- proj_prop
        }

        basis_mat_prop <- basis_mat
        basis_mat_prop[, basis_idx[[j_change + 1]]] <- ridge_basis_prop
        qf_info_prop <- get_qf_info(basis_mat_prop, y)

        if(!is.null(qf_info_prop)  &&  (sse_prop <- ssy - var_coefs[it]/(var_coefs[it] + 1) * qf_info_prop$qf) > 0){
          # Compute the acceptance probability
          alpha <- -n/2 * (log(sse_prop) - log(sse)) # Marginal Likelihood

          if(log(runif(1)) < alpha){ # Accept the proposal
            knots[[it]][[j_change]] <- knots_prop
            bias[[it]][j_change] <- bias_prop
            proj_dir[[it]][[j_change]] <- proj_dir_prop

            basis_mat <- basis_mat_prop

            qf_info <- qf_info_prop
            qf_info <- append_qf_inv_chol(qf_info, dim = n_basis_total)
            sse <- sse_prop
          }
        }
      }
    }

    # Draw coefs
    coefs[[it]] <- var_coefs[it]/(var_coefs[it] + 1) * qf_info$ls_est +
      sqrt(var_coefs[it]/(var_coefs[it] + 1)) * sd_resid[it] * qf_info$inv_chol %*% rnorm(n_basis_total) # Draw coefs
    preds <- basis_mat %*% coefs[[it]] # Current predictions of y
    resid <- y - preds # current residuals

    sd_resid[it] <- sqrt(1/rgamma(1, n/2, c(t(resid) %*% resid)/2))

    var_coefs[it] <- 1/rgamma(1,
                              var_coefs_shape + n_basis_total/2,
                              var_coefs_rate + c(t(preds) %*% preds)/(2*sd_resid[it]^2))

    sse <- ssy - var_coefs[it]/(var_coefs[it] + 1) * qf_info$qf
  }

  structure(list(n_ridge = n_ridge, n_act = n_act, feat = feat, proj_dir = proj_dir, bias = bias, knots = knots,
                 coefs = coefs, var_coefs = var_coefs, sd_resid = sd_resid,
                 n_act_w = n_act_w, feat_w = feat_w,
                 mn_X = mn_X, sd_X = sd_X,
                 df_spline = df_spline, n_ridge_mean = n_ridge_mean, model = model),
            class = 'bppr')
}
