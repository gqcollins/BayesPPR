########################################################################
## BayesPPR for multivariate response
########################################################################

#' @title BayesPPR for Multivariate Response
#'
#' @description Fits a BayesPPR model to multivariate response using RJMCMC. Can handle categorical features.
#' @param X a data frame or matrix of predictors. Categorical features should be coded as numeric.
#' @param Y a data frame or matrix of responses, with columns representing different dimensions of the response.
#' @param n_pc number of principle components to be used in pca approximation of Y.
#' @param prop_var proportion of variance to be explained by the first \code{n_pc} principle components. Used for automatic selection of \code{n_pc} if \code{is.null(n_pc)}.
#' @param n_cores number of cores the user desired to utilize for parallel computation.
#' @param par_type a character variable dicatating the type of parallel computation to be performed. Supported values include \code{"fork"}, which uses \code{parallel::mclapply()}, and \code{"socket"}, which uses \code{parallel::parLapply()}.
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
#' @param bppr_init_list list of length \code{n_pc}, each element containing a list initial values for the Markov chain. Used by \link{bppr_resume}.
#' @details Explores BayesPPR model space using RJMCMC. The BayesPPR model has \deqn{y = f(x) + \epsilon,  ~~\epsilon \sim N(0,\sigma^2)} \deqn{f(x) = \beta_0 + \sum_{j=1}^M \beta_j B_j(x)} and \eqn{B_j(x)} is a natural spline basis expansion. We use priors \deqn{\beta \sim N(0,\sigma^2/\tau (B'B)^{-1})} \deqn{M \sim Poisson(\lambda)} as well as the hyper-prior on the variance \eqn{\tau} of the coefficients \eqn{\beta} mentioned in the arguments above.
#' @return An object of class \code{"bppr_pca"}. Predictions can be obtained by passing the entire object to the \code{predict.bppr_pca} function.
#' @keywords nonparametric projection pursuit regression splines principle component analysis
#' @seealso \link{predict.bppr_pca} for prediction.
#' @export
#' @import stats
#' @import utils
#' @examples
#' # See examples in bppr documentation.
#'
bppr_pca <- function(X, Y, n_pc = NULL, prop_var = 0.99, n_cores = 1, par_type = 'fork', n_ridge_mean = 10, n_ridge_max = NULL, n_act_max = NULL, df_spline = 4, prob_relu = 2/3, prior_coefs = "zs", shape_var_coefs = NULL, scale_var_coefs = NULL, n_dat_min = NULL, scale_proj_dir_prop = NULL, adapt_act_feat = TRUE, w_n_act = NULL, w_feat = NULL, n_post = 1000, n_burn = 9000, n_adapt = 0, n_thin = 1, print_every = NULL, bppr_init_list = NULL){

  pca_Y <- pca_setup(X, Y, n_pc = n_pc, prop_var = prop_var)
  n_pc <- pca_Y$n_pc

  n_cores_max <- parallel::detectCores()

  if(n_cores > n_cores_max){
    warning(paste0("Specified n_cores = ", n_cores, ". Proceeding with n_cores = min(n_cores, n_pc, detectCores()) = ",
                   min(n_cores, n_pc, n_cores_max)))
  }
  n_cores <- min(n_cores, n_pc, n_cores_max)

  if(is.null(print_every)){
    if(n_cores == 1){
      print_every <- 1000
    }else{
      print_every <- 0
    }
  }

  if(is.null(bppr_init_list)){
    bppr_init_code <- "NULL"
  }else{
    bppr_init_code <- "bppr_init_list[[i]]"
  }

  run_bppr <- parse(text = paste0(
  "bppr(X, pca_Y$Y_new[, i], n_ridge_mean = n_ridge_mean, n_ridge_max = n_ridge_max,
        n_act_max = n_act_max, df_spline = df_spline, prob_relu = prob_relu,
        prior_coefs = prior_coefs, shape_var_coefs = shape_var_coefs,
        scale_var_coefs = scale_var_coefs, n_dat_min = n_dat_min,
        scale_proj_dir_prop = scale_proj_dir_prop, adapt_act_feat = adapt_act_feat,
        w_n_act = w_n_act, w_feat = w_feat, n_post = n_post, n_burn = n_burn,
        n_adapt = n_adapt, n_thin = n_thin, print_every = print_every,
        bppr_init = ", bppr_init_code, ")")
  )

  if(n_cores == 1){
    fit_list <- lapply(1:n_pc, function(i) eval(run_bppr))
  }else if(par_type == "socket"){
    cl <- parallel::makeCluster(n_cores, setup_strategy = "sequential")
    fit_list <- parallel::parLapply(cl, 1:n_pc, function(i) eval(run_bppr))
    parallel::stopCluster(cl)
  }else if(par_type == "fork"){
    fit_list <- parallel::mclapply(1:pca_Y$n_pc, function(i) eval(run_bppr),
      mc.cores = n_cores, mc.preschedule = FALSE)
  }

 structure(list(pca_Y = pca_Y, fit_list = fit_list,
                n_cores = n_cores, par_type = par_type),
           class = 'bppr_pca')
}
