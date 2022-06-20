########################################################################
## prediction function
########################################################################

#' @title BayesPPR Predictions
#'
#' @description Predict function for BayesPPR  Outputs the posterior predictive samples for the desired MCMC iterations.
#' @param object a fitted model, output from the \code{bppr} function.
#' @param newdata a matrix of new input values at which to predict.  The columns should correspond to the same variables used in the \code{bppr} function.
#' @param idx_use index of Markov samples to use when generating predictions.
#' @param ... further arguments passed to or from other methods.
#' @details bare-bones methods. Could be improved for efficiency.
#' #' @return This returns a matrix with the same number of rows as \code{newdata} and columns corresponding to all MCMC iterations after the first \code{n_ignore}.  These are samples from the posterior predictive distribution.
#' @seealso \link{bppr} for model fitting.
#' @export
#' @examples
#' # See examples in bppr documentation.
#'
predict.bppr <- function(object, newdata, idx_use = NULL,...){
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

  n_keep <- length(proj_dir)
  if(is.null(idx_use)){
    idx_use <- 1:n_keep
    n_use <- n_keep
  }else if(max(idx_use) > n_keep){
    stop("invalid 'idx_use'")
  }else{
    n_use <- length(idx_use)
  }

  preds <- matrix(0, nrow = n_use, ncol = n)
  basis_int <- matrix(rep(1, n))
  for(i in 1:n_use){
    basis_mat <- basis_int
    if(n_ridge[idx_use[i]] > 0){
      for(j in 1:n_ridge[idx_use[i]]){
        if(is.na(proj_dir[[idx_use[i]]][[j]][1])){ # No quantitative feature in this basis
          ridge_basis <- get_cat_basis(newdata[, feat[[idx_use[i]]][[j]], drop = FALSE])
        }else if(is.na(knots[[idx_use[i]]][[j]][1])){ # No continuous features in this basis
          ridge_basis <- newdata[, feat[[idx_use[i]]][[j]], drop = FALSE] %*% proj_dir[[idx_use[i]]][[j]]
        }else{ # At least one continuous variable in this basis
          if(is.na(bias[[idx_use[i]]][j])){ # No knot or relu in this basis
            proj <- newdata[, feat[[idx_use[i]]][[j]], drop = FALSE] %*% proj_dir[[idx_use[i]]][[j]]
          }else{ # Relu present
            proj <- relu(bias[[idx_use[i]]][j] + newdata[, feat[[idx_use[i]]][[j]], drop = FALSE] %*% proj_dir[[idx_use[i]]][[j]]) # Get relu of projection
          }
          ridge_basis <- get_ns_basis(proj, knots[[idx_use[i]]][[j]]) # Get basis function
        }
        basis_mat <- cbind(basis_mat, ridge_basis)
      }
    }
    preds[i, ] <- basis_mat %*% coefs[[idx_use[i]]] # Get predictions for jth basis function
  }

  return(preds)
}
