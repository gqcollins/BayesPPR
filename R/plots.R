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
plot.bppr<-function(x,quants=c(.025,.975),pred=TRUE,...){
  if(class(x)!='bppr')
    stop('x must be an object of class bppr')

  op<-par(no.readonly=T)
  if(pred)
    par(mfrow=c(2,2))
  else
    par(mfrow=c(1,2))
  plot(x$n_ridge,type='l',ylab='number of ridge functions',xlab='MCMC iteration (post-burn)')
  plot(x$sd_resid,type='l',ylab='residual sd',xlab='MCMC iteration (post-burn)')
  if(pred){
    margin<-2
    preds<-predict(x, x$X)
    yhat<-colMeans(preds)
    if(!is.null(quants)){
      qq<-apply(preds+rnorm(prod(dim(preds)), sd=x$sd_resid),margin,quantile,probs=quants)
      ylim=range(qq)
      ylab='interval'
    } else{
      ylim=range(yhat)
      ylab='mean'
    }
    plot(x$y,yhat,ylim=ylim,ylab=paste('posterior predictive',ylab),xlab='observed',main='Training Fit',type='n',...)
    if(!is.null(quants))
      segments(x$y,qq[1,],x$y,qq[2,],col='lightgrey')
    points(x$y,yhat)
    abline(a=0,b=1,col=2)

    hist(x$y-yhat,freq=F,main='Posterior mean residuals',xlab='residuals')
    curve(dnorm(xx,sd=mean(x$sd_resid)),xname='xx',col=2,add=T)
  }
  mtext('BPPR Diagnostics',3,-2,outer=T)
  par(op)
}

