#' @title BPPR Plot Diagnostics
#'
#' @description Generate diagnostic plots for BPPR model fit.
#' @param x a \code{bppr} object.
#' @param quants quantiles for intervals, if desired.  NULL if not desired.
#' @param pred logical, whether posterior predictions should be plotted (defaults to TRUE).
#' @param ... graphical parameters.
#' @details The first two plots are trace plots for diagnosing convergence.  The third plot is posterior predicted vs observed, with intervals for predictions.  The fourth plot is a histogram of the residuals (of the posterior mean model), with a red curve showing the assumed Normal density (using posterior mean variance). If \code{pred=FALSE} the third and fourth plots are omitted.
#' @export
#' @import graphics
#' @seealso \link{bppr}, \link{predict.bppr}
#' @examples
#' # See examples in bppr documentation.
#'
plot.bppr <- function(x, quants = c(.025, .975), pred = TRUE, ...){
  if(class(x) != 'bppr'){
    stop('x must be an object of class bppr')
  }

  op <- par(no.readonly = TRUE)
  if(pred){
    par(mfrow = c(2, 2))
  }else{
    par(mfrow = c(1, 2))
  }
  plot(x$n_ridge, type = 'l', ylab = 'Number of Ridge Functions', xlab = 'MCMC Iteration (Post-Burn)')
  plot(x$sd_resid, type = 'l', ylab = 'Residual SD', xlab = 'MCMC Iteration (Post-Burn)')
  if(pred){
    margin <- 2
    preds <- predict(x, x$X)
    yhat <- colMeans(preds)
    if(!is.null(quants)){
      qq <- apply(preds + rnorm(prod(dim(preds)), sd = x$sd_resid), margin, quantile, probs = quants)
      ylim <- range(qq)
      ylab <- 'Interval'
    } else{
      ylim <- range(yhat)
      ylab <- 'Mean'
    }
    plot(x$y, yhat, ylim = ylim, ylab = paste('Posterior Predictive', ylab), xlab = 'Observed Response', main = 'Training Fit', type = 'n', ...)
    if(!is.null(quants)){
      segments(x$y, qq[1,], x$y, qq[2,], col='lightgrey')
    }
    points(x$y, yhat)
    abline(a = 0, b = 1, col = 2)

    hist_dat <- hist(x$y - yhat, plot = FALSE)
    xx <- seq(min(hist_dat$breaks), max(hist_dat$breaks), length.out = 200)
    mn_sd <- mean(x$sd_resid)
    plot(hist_dat, freq = FALSE, col = 'lightgrey', main = 'Posterior Mean Residuals',
         xlab = 'Residuals', ylim = c(0, max(hist_dat$density, dnorm(0, sd = mn_sd))))
    lines(xx, dnorm(xx, sd = mn_sd), col = 2)
  }
  mtext('BPPR Diagnostics', 3, -2, outer=TRUE)
  par(op)
}

