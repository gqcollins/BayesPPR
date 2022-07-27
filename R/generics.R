
#' @title Print BPPR Details
#'
#' @description Print some of the details of a BPPR model.
#' @param x a \code{bppr} object, returned from \code{bppr}.
#' @param ... further arguments passed to or from other methods.
#' @export
#'
print.bppr<-function(x,...){
  cat("\nCall:\n",deparse(x$call),'\n')

  cat("\n              Number of variables: ",ncol(x$X),sep='')
  cat("\n                      Sample size: ",nrow(x$X),sep='')

  cat('\n\n')

}


#' @title Summarize BPPR Details
#'
#' @description Summarize some of the details of a BPPR model.
#' @param object a \code{bppr} object, returned from \code{bppr}.
#' @param ... further arguments passed to or from other methods.
#' @export
#'
summary.bppr<-function(object,...){
  cat("\nCall:\n",deparse(object$call),'\n')

  cat("\n              Number of variables: ",ncol(object$X),sep='')
  cat("\n                      Sample size: ",nrow(object$X),sep='')

  cat("\n\nNumber of ridge functions (range):",range(object$n_ridge))
  cat("\n    Posterior mean error sd:",mean(object$sd_resid),'\n\n')

}


