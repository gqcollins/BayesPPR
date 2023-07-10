get_qf_info <- function(BtB, Bty){
  chol_BtB <- tryCatch(chol(BtB), error = function(e) matrix(FALSE))
  if(chol_BtB[1, 1]){
    d <- diag(chol_BtB)
    if(length(d) > 1){
      if((max(d[-1]) / min(d)) > 1000){
        return(NULL)
      }
    }
    bhat <- backsolve(chol_BtB, forwardsolve(chol_BtB, Bty,
                                             transpose = TRUE,
                                             upper.tri = TRUE))
    qf <- t(Bty) %*% bhat
    return(list(chol = chol_BtB, ls_est = bhat, qf = qf))
  }else{
    return(NULL)
  }
}

append_qf_inv_chol <- function(qf_info, dim){
  qf_info$inv_chol <- backsolve(qf_info$chol, diag(dim))
  return(qf_info)
}


