########################################################################
## prediction function
########################################################################

#' @title BayesPPR Predictions
#'
#' @description Predict function for BayesPPR. Outputs the posterior predictive samples for the desired MCMC iterations.
#' @param object a fitted model, output from the \code{bppr} function.
#' @param newdata a matrix of new input values at which to predict.  The columns should correspond to the same variables used in the \code{bppr} function.
#' @param idx_use index of Markov samples to use when generating predictions.
#' @param ... further arguments passed to or from other methods.
#' @details Bare-bones methods. Could be improved for efficiency.
#' @return A matrix with the same number of rows as \code{newdata} and columns corresponding to all MCMC iterations indexed by idx_use. These are samples from the posterior predictive distribution.
#' @seealso \link{bppr} for model fitting.
#' @export
#' @examples
#' # See examples in bppr documentation.
#'
predict.bppr <- function(object, newdata, idx_use = NULL, ...){
  newdata <- as.matrix(newdata)
  n <- nrow(newdata)
  p <- ncol(newdata)
  mn_X <- object$mn_X
  sd_X <- object$sd_X
  for(j in 1:p){
    newdata[, j] <- (newdata[, j] - mn_X[j]) / sd_X[j]
  }
  coefs <- object$coefs
  bias <- object$bias
  proj_dir <- object$proj_dir
  n_ridge <- object$n_ridge
  n_act <- object$n_act
  feat <- object$feat
  knots <- object$knots
  df_spline <- object$df_spline
  n_keep <- object$n_keep

  if(is.null(idx_use)){
    idx_use <- 1:n_keep
  }else if(max(idx_use) > n_keep){
    stop("invalid 'idx_use'")
  }
  n_use <- length(idx_use)

  preds <- matrix(0, nrow = n_use, ncol = n)
  for(i in 1:n_use){
    preds[i, ] <- coefs[[idx_use[i]]][1]
    if(n_ridge[idx_use[i]] > 0){
      basis_idx_start <- 2
      for(j in 1:n_ridge[idx_use[i]]){
        if(is.na(proj_dir[[idx_use[i]]][[j]][1])){ # No quantitative feature in this basis
          ridge_basis <- get_cat_basis(newdata[, feat[[idx_use[i]]][[j]], drop = FALSE])
          basis_idx <- basis_idx_start
          basis_idx_start <- basis_idx_start + 1
        }else if(is.na(knots[[idx_use[i]]][[j]][1])){ # No continuous features in this basis
          ridge_basis <- newdata[, feat[[idx_use[i]]][[j]], drop = FALSE] %*% proj_dir[[idx_use[i]]][[j]]
          basis_idx <- basis_idx_start
          basis_idx_start <- basis_idx_start + 1
        }else{ # At least one continuous variable in this basis
          proj <- relu(bias[[idx_use[i]]][j] + newdata[, feat[[idx_use[i]]][[j]], drop = FALSE] %*% proj_dir[[idx_use[i]]][[j]]) # Get relu of projection
          ridge_basis <- get_ns_basis(proj, knots[[idx_use[i]]][[j]]) # Get basis function
          basis_idx <- basis_idx_start:(basis_idx_start + df_spline - 1)
          basis_idx_start <- basis_idx_start + df_spline
        }
        preds[i, ] <- preds[i, ] + ridge_basis %*% coefs[[idx_use[i]]][basis_idx] # Get predictions for jth basis function
      }
    }
  }

  return(preds)
}
