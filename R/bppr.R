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
#' @param scale_var_coefs scale for IG prior on the variance of the basis function coefficients. Default is for the Zellner-Siow prior. For the flat, improper prior, \code{scale_var_coefs} is ignored.
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
#' @param intercept logical; if \code{FALSE}, the model has no intercept and the basis functions of every ridge function are centred at their training means, so every draw of \eqn{f} averages zero over the training inputs.
#' @param s2_shape,s2_rate shape and rate of an inverse-gamma prior on the residual variance \eqn{\sigma^2}. The default \code{0, 0} is the improper prior \eqn{1/\sigma^2}.
#' @param s2_fixed if not \code{NULL}, \eqn{\sigma^2} is fixed at this value and not sampled.
#' @param keep_sampler logical; if \code{TRUE}, the returned object carries the live chain in \code{$sampler} (a \link{bppr_sampler}), so it can be continued, with a new response if desired, by \link{bppr_resume} without rebuilding its state.
#' @details Explores BayesPPR model space using RJMCMC. The BayesPPR model has \deqn{y = f(x) + \epsilon,  ~~\epsilon \sim N(0,\sigma^2)} \deqn{f(x) = \beta_0 + \sum_{j=1}^M \beta_j B_j(x)} and \eqn{B_j(x)} is a natural spline basis expansion. We use priors \deqn{\beta \sim N(0,\sigma^2/\tau (B'B)^{-1})} \deqn{M \sim Poisson(\lambda)} as well as the hyper-prior on the variance \eqn{\tau} of the coefficients \eqn{\beta} mentioned in the arguments above.
#' @return An object of class \code{"bppr"}. Predictions can be obtained by passing the entire object to the \code{predict.bppr} function.
#' @keywords nonparametric projection pursuit regression splines
#' @seealso \link{predict.bppr} for prediction.
#' @export
#' @import stats
#' @import utils
#' @example inst/examples.R
#'
bppr <- function(X, y, n_ridge_mean = 10, n_ridge_max = NULL, n_act_max = NULL, df_spline = 4, prob_relu = 2/3, prior_coefs = "zs", shape_var_coefs = NULL, scale_var_coefs = NULL, n_dat_min = NULL, scale_proj_dir_prop = NULL, adapt_act_feat = TRUE, w_n_act = NULL, w_feat = NULL, n_post = 1000, n_burn = 9000, n_adapt = 0, n_thin = 1, print_every = 1000, bppr_init = NULL, intercept = TRUE, s2_shape = 0, s2_rate = 0, s2_fixed = NULL, keep_sampler = FALSE){
  # Manage posterior draws
  if(n_thin > n_post){
    stop('n_thin > n_post. No posterior samples will be obtained.')
  }
  n_post <- n_post - n_post %% n_thin
  n_pre <- n_adapt + n_burn
  if(prior_coefs == 'flat'){
    n_adapt <- n_pre
    n_burn <- 0
  }

  sampler <- bppr_sampler(X, y, n_ridge_mean = n_ridge_mean, n_ridge_max = n_ridge_max, n_act_max = n_act_max,
                          df_spline = df_spline, prob_relu = prob_relu, prior_coefs = prior_coefs,
                          shape_var_coefs = shape_var_coefs, scale_var_coefs = scale_var_coefs,
                          n_dat_min = n_dat_min, scale_proj_dir_prop = scale_proj_dir_prop,
                          adapt_act_feat = adapt_act_feat, w_n_act = w_n_act, w_feat = w_feat,
                          intercept = intercept, s2_shape = s2_shape, s2_rate = s2_rate,
                          s2_fixed = s2_fixed, bppr_init = bppr_init)
  print_every <- run_bppr_chain(sampler, n_post = n_post, n_burn = n_burn, n_adapt = n_adapt,
                                n_thin = n_thin, print_every = print_every)

  out <- sampler$as_bppr(n_post = n_post, n_burn = n_burn,
                         n_adapt = n_adapt, n_thin = n_thin, print_every = print_every,
                         bppr_init = bppr_init, call = match.call())
  if(keep_sampler) out$sampler <- sampler
  return(out)
}

# Run a chain of n_adapt + n_burn + n_post iterations on a sampler, storing every
# n_thin-th post-burn state. Iteration 1 is the sampler's current state, so
# n_adapt + n_burn + n_post - 1 moves are made. Returns the print_every stored
# in the fit (n_draws + 2 when silent).
run_bppr_chain <- function(sampler, n_post, n_burn, n_adapt, n_thin, print_every){
  n_pre <- n_adapt + n_burn
  n_draws <- n_pre + n_post

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
               myTimestamp(), ' n ridge: ', sampler$n_ridge, '\n'))
    silent <- FALSE
  }else{
    print_every <- n_draws + 2
    silent <- TRUE
  }

  for(it in seq_len(n_draws)){
    if(it > 1){
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
                     myTimestamp(start_time), ' n ridge: ', sampler$n_ridge)
        cat(pr, '\n')
      }
      sampler$iterate(sample_params = it > n_adapt)
    }
    if(it > n_pre  &&  (it - n_pre) %% n_thin == 0) sampler$save()
  }

  if(!silent){
    cat(paste0('MCMC iteration ', n_draws, '/', n_draws, ' (', phase, ') ',
               myTimestamp(start_time), ' n ridge: ', sampler$n_ridge, '\n'))
  }
  print_every
}
