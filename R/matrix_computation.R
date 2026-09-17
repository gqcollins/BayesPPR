get_qf_info <- function(BtB, Bty){
  chol_BtB <- tryCatch(chol(BtB), error = function(e) matrix(FALSE))
  if(chol_BtB[1, 1]  &&  !any(is.na(chol_BtB))){
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

check_idx_use <- function(idx_use, n_keep){ # Validate an index into the kept MCMC draws
  if(!is.numeric(idx_use)  ||  length(idx_use) == 0  ||  anyNA(idx_use)){
    stop("'idx_use' must be a non-empty numeric vector without missing values")
  }
  if(any(idx_use != round(idx_use))){
    stop("'idx_use' must contain whole numbers")
  }
  idx_use <- as.integer(round(idx_use))
  if(any(idx_use < 1)  ||  any(idx_use > n_keep)){
    stop(paste0("'idx_use' must be between 1 and ", n_keep))
  }
  return(idx_use)
}


