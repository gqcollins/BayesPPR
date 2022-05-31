get_qf_info <- function(X, y){
  Xt <- t(X)
  XtX <- Xt %*% X
  chol_XtX <- tryCatch(chol(XtX), error = function(e) matrix(FALSE))
  if(chol_XtX[1, 1]){
    d <- diag(chol_XtX)
    if(length(d) > 1){
      if((max(d[-1]) / min(d)) > 1000){
        return(NULL)
      }
    }
    Xty <- Xt %*% y
    bhat <- backsolve(chol_XtX, forwardsolve(chol_XtX, Xty,
                                             transpose = TRUE,
                                             upper.tri = TRUE))
    qf <- t(Xty) %*% bhat
    return(list(chol = chol_XtX, ls_est = bhat, qf = qf))
  }else{
    return(NULL)
  }
}

append_qf_inv_chol <- function(qf_info, dim){
  qf_info$inv_chol <- backsolve(qf_info$chol, diag(dim))
  return(qf_info)
}


