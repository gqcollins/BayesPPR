########################################################################
## prediction function for multivariate response
########################################################################

#' @title BayesPPR Predictions for Multivariate Response
#'
#' @description Predict function for BayesPPR. Outputs the posterior predictive samples for the desired MCMC iterations.
#' @param object a fitted model, output from the \code{bppr} function.
#' @param newdata a matrix of new input values at which to predict.  The columns should correspond to the same variables used in the \code{bppr} function.
#' @param idx_use index of Markov samples to use when generating predictions.
#' @param n_cores number of cores the user desired to utilize for parallel computation.
#' @param par_type a character variable dicatating the type of parallel computation to be performed. Supported values include \code{"fork"}, which uses \code{parallel::mclapply()}, and \code{"socket"}, which uses \code{parallel::parLapply()}.
#' @param ... further arguments passed to or from other methods.
#' @details Bare-bones methods. Could be improved for efficiency.
#' @return An matrix where the first dimension corresponds to MCMC iterations indexed by idx_use, the second dimensions is the same number of rows as \code{newdata}, and the third corresponds to the dimension of the response.
#' @seealso \link{bppr} for model fitting.
#' @export
#' @examples
#' # See examples in bppr documentation.
#'
predict.bppr_pca <-function(object, newdata, idx_use = NULL, n_cores = 1, par_type = 'fork', ...){

  n_keep <- object$fit_list[[1]]$n_keep
  if(is.null(idx_use)){
    idx_use <- 1:n_keep
  }else if(max(idx_use) > n_keep){
    stop("invalid 'idx_use'")
  }
  n_use <- length(idx_use)

  n <- nrow(newdata)
  D <- length(object$pca_Y$mn_Y)

  n_pc <- object$pca_Y$n_pc

  run_predict <- parse(text = 'predict(object$fit_list[[i]], newdata, idx_use = idx_use, ...)')
  run_pca_reverse <- parse(text = 'pca_reverse(preds_Y_new[i, , ], object$pca_Y)')

  if(n_cores == 1){
    preds_Y_new <- array(
      unlist(lapply(1:n_pc, function(i) eval(run_predict))),
      dim = c(n_use, n, n_pc))
    out <- array(
      unlist(lapply(1:n_use, function(i) eval(run_pca_reverse))),
      dim = c(D, n, n_use))
  }else if(par_type == "socket"){
    cl <- parallel::makeCluster(min(n_cores, n_pc, parallel::detectCores()),
                                setup_strategy = "sequential")
    parallel::clusterExport(cl, varlist = c("newdata"), envir = environment())
    preds_Y_new <- array(
      unlist(parallel::parLapply(cl, 1:n_pc, function(i) eval(run_predict))),
      dim = c(n_use, n, n_pc))
    out <- array(
      unlist(parallel::parLapply(cl, 1:n_use, function(i) eval(run_pca_reverse))),
      dim = c(D, n, n_use))
    parallel::stopCluster(cl)
  }else if(par_type == "fork"){
    preds_Y_new <- array(
      unlist(parallel::mclapply(1:n_pc, function(i) eval(run_predict),
                                mc.cores = n_cores)),
      dim = c(n_use, n, n_pc))
    out <- array(
      unlist(parallel::mclapply(1:n_use, function(i) eval(run_pca_reverse),
             mc.cores = n_cores)),
      dim = c(D, n, n_use))
  }

  out <- aperm(out, c(3, 2, 1))
  return(out)
}
