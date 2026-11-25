# Cholesky factor of BtB, least-squares estimate and quadratic form Bty' BtB^-1 Bty.
# NULL if BtB is singular or badly conditioned. With an intercept the first pivot
# (the intercept's) is left out of the conditioning check. An empty basis is the
# zero model.
get_qf_info <- function(BtB, Bty, intercept = TRUE){
  if(length(Bty) == 0){
    return(list(chol = matrix(0, 0, 0), ls_est = numeric(0), qf = 0))
  }
  chol_BtB <- tryCatch(chol(BtB), error = function(e) matrix(FALSE))
  if(chol_BtB[1, 1]  &&  !any(is.na(chol_BtB))){
    d <- diag(chol_BtB)
    if(length(d) > 1){
      d_check <- if(intercept) d[-1] else d
      if((max(d_check) / min(d)) > 1000){
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

# Same as get_qf_info() for a new Bty, reusing the Cholesky factor (and its
# inverse, if present), which depend only on the basis
update_qf_info <- function(qf_info, Bty){
  if(length(Bty) == 0){
    qf_info$ls_est <- numeric(0)
    qf_info$qf <- 0
    return(qf_info)
  }
  qf_info$ls_est <- backsolve(qf_info$chol, forwardsolve(qf_info$chol, Bty,
                                                         transpose = TRUE,
                                                         upper.tri = TRUE))
  qf_info$qf <- t(Bty) %*% qf_info$ls_est
  return(qf_info)
}

append_qf_inv_chol <- function(qf_info, dim){
  qf_info$inv_chol <- if(dim == 0) matrix(0, 0, 0) else backsolve(qf_info$chol, diag(dim))
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


