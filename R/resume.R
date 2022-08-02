########################################################################
## Resume the MCMC given a current bppr object
########################################################################

#' @title Resume Sampling of BPPR Parameters
#'
#' @description Resumes fitting a BayesPPR model where the last Markov chain left off, using the same model parameters.
#' @param object an object of class \code{"bppr"} containing the values of the current Markov chain, the latest of which constitute initial values for the new chain.
#' @param append logical; if \code{TRUE}, new samples will be appended to samples in \code{object}. Otherwise, output will only include new samples.
#' @param n_post number of posterior draws to obtain from the Markov chain after burn-in.
#' @param n_burn number of draws to burn before obtaining \code{n_post} draws for inference. If \code{prior_coefs == "flat"} then these are absorbed into the adapt phase.
#' @param n_adapt number of adaptive MCMC iterations to perform before burn-in. Skips sampling basis coefficients and residual variance to save time.
#' @param n_thin keep every n_thin posterior draws after burn-in.
#' @param print_every print the iteration number every print_every iterations. Use \code{print_every = 0} to silence.
#' @details Continues exploring BayesPPR model space using RJMCMC, beginning at the values contained in the last iteration of \code{object}.
#' @return An object of class \code{"bppr"}. Predictions can be obtained by passing the entire object to the \code{predict.bppr} function.
#' @keywords nonparametric projection pursuit regression splines
#' @seealso \link{predict.bppr} for prediction.
#' @export
#' @import stats
#' @import utils
#' @example
#' # See examples in bppr documentation.
#'
bppr_resume <- function(object, append = FALSE, n_post = 1000, n_burn = 9000, n_adapt = 0, n_thin = 1, print_every = 1000){
  n_keep <- object$n_keep
  if(object$prior_coefs == 'zs'){
    var_coefs_init <- object$var_coefs[n_keep]
  }else{
    var_coefs_init <- NULL
  }

  bppr_init <- list(n_ridge = object$n_ridge[[n_keep]],
                    n_act = object$n_act[[n_keep]], feat = object$feat[[n_keep]],
                    proj_dir = object$proj_dir[[n_keep]], knots = object$knots[[n_keep]],
                    coefs = object$coefs[[n_keep]], sd_resid = object$sd_resid[n_keep],
                    var_coefs = var_coefs_init)

  fit <- bppr(X = object$X, y = object$y, n_ridge_mean = object$n_ridge_mean,
              n_ridge_max = object$n_ridge_max, n_act_max = object$n_act_max,
              df_spline = object$df_spline, prob_relu = object$prob_relu,
              prior_coefs = object$prior_coefs, shape_var_coefs = object$shape_var_coefs,
              rate_var_coefs = object$rate_var_coefs, n_dat_min = object$n_dat_min,
              scale_proj_dir_prop = object$scale_proj_dir_prop,
              adapt_act_feat = object$adapt_act_feat, w_n_act = object$w_n_act,
              w_feat = object$w_feat, n_post = n_post, n_burn = n_burn,
              n_adapt = n_adapt, n_thin = n_thin, print_every = print_every,
              bppr_init = bppr_init)

  if(append){
    fit$n_ridge <- c(object$n_ridge, fit$n_ridge)
    fit$n_act <- c(object$n_act, fit$n_act)
    fit$feat <- c(object$feat, fit$feat)
    fit$proj_dir <- c(object$proj_dir, fit$proj_dir)
    fit$knots <- c(object$knots, fit$knots)
    fit$coefs <- c(object$coefs, fit$coefs)
    fit$sd_resid <- c(object$sd_resid, fit$sd_resid)
    fit$var_coefs <- c(object$var_coefs, fit$var_coefs)
    fit$n_keep <- object$n_keep + fit$n_keep
    fit$n_post <- c(object$n_post, fit$n_post)
    fit$n_burn <- c(object$n_burn, fit$n_burn)
    fit$n_adapt <- c(object$n_adapt, fit$n_adapt)
  }

  return(fit)
}
