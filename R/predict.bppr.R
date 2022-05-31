########################################################################
## prediction function
########################################################################

#' @title BayesPPR Predictions
#'
#' @description Predict function for BayesPPR  Outputs the posterior predictive samples for the desired MCMC iterations.
#' @param object a fitted model, output from the \code{bppr} function.
#' @param newdata a matrix of new input values at which to predict.  The columns should correspond to the same variables used in the \code{bppr} function.
#' @param n_ignore number of samples to ignore at the beginning of the Markov chain, when generating predictions.
#' @param ... further arguments passed to or from other methods.
#' @details bare-bones methods. Could be improved for efficiency.
#' #' @return This returns a matrix with the same number of rows as \code{newdata} and columns corresponding to all MCMC iterations after the first \code{n_ignore}.  These are samples from the posterior predictive distribution.
#' @seealso \link{bppr} for model fitting.
#' @export
#' @examples
#' # See examples in bppr documentation.
#'
predict.bppr <- function(object, newdata, n_ignore = 0,...){
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

  n_draws <- length(proj_dir)
  if(n_ignore >= n_draws){
    stop('n_ignore >= n_draws')
  }
  n_keep <- n_draws - n_ignore
  it_keep <- (n_ignore + 1):n_draws
  preds <- matrix(0, nrow = n_keep, ncol = n)
  basis_int <- matrix(rep(1, n))
  for(it in it_keep){
    basis_mat <- basis_int
    if(n_ridge[it] > 0){
      for(j in 1:n_ridge[it]){
        if(is.na(proj_dir[[it]][[j]][1])){ # No quantitative feature in this basis
          ridge_basis <- get_cat_basis(newdata[, feat[[it]][[j]], drop = FALSE])
        }else if(is.na(knots[[it]][[j]][1])){ # No continuous features in this basis
          ridge_basis <- newdata[, feat[[it]][[j]], drop = FALSE] %*% proj_dir[[it]][[j]]
        }else{ # At least one continuous variable in this basis
          if(is.na(bias[[it]][j])){ # No knot or relu in this basis
            proj <- newdata[, feat[[it]][[j]]] %*% proj_dir[[it]][[j]]
          }else{ # Relu present
            proj <- relu(bias[[it]][j] + newdata[, feat[[it]][[j]]] %*% proj_dir[[it]][[j]]) # Get relu of projection
          }
          ridge_basis <- get_ns_basis(proj, knots[[it]][[j]]) # Get basis function
        }
        basis_mat <- cbind(basis_mat, ridge_basis)
      }
    }
    preds[it - n_ignore, ] <- basis_mat %*% coefs[[it]] # Get predictions for jth basis function
  }

  return(preds)
}
