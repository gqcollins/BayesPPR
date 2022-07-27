pca_setup <- function(X, Y, n_pc = NULL, prop_var = 0.99){
  if(prop_var > 1  ||  prop_var < 0){
    stop('prop_var must be between 0 and 1')
  }

  n <- nrow(X)
  Y <- as.matrix(Y)

  if(nrow(Y) == n){
    D <- ncol(Y)
  }else{
    if(ncol(Y) == n){
      Y <- t(Y)
      D <- ncol(Y)
    }else{
      stop('X, Y dimension mismatch')
    }
  }

  if(D == 1){
    stop('univariate Y: use bppr instead of bppr_pca')
  }

  if(n == D){
    warning("Caution: because Y is square, please ensure that each row of X corresponds to a row of Y (and not a column)")
  }

  if(!is.null(n_pc)){
    if(n_pc > max(n, D)){
      n_pc <- max(n, D)
      warning('n_pc too large, using all PCs intead')
    }
  }

  mn_Y <- sd_Y <- numeric(D)
  Y_std <- matrix(nrow = n, ncol = D)
  for(k in 1:D){
    mn_Y[k] <- mean(Y[, k])
    sd_Y[k] <- sd(Y[, k])
    Y_std[, k] <- (Y[, k] - mn_Y[k]) / sd_Y[k]
  }

  svd_Y <- svd(t(Y_std))
  ev_Y <- svd_Y$d^2

  if(is.null(n_pc)){
    n_pc <- which(cumsum(ev_Y / sum(ev_Y)) > prop_var)[1]
  }

  basis_Y <- svd_Y$u[, 1:n_pc, drop=FALSE] %*% diag(svd_Y$d[1:n_pc], nrow = n_pc) # columns are basis functions
  Y_new <- svd_Y$v[, 1:n_pc, drop=FALSE]

  trunc_error <- Y_new %*% t(basis_Y) - Y_std

  out <- list(X = X, Y = Y, n_pc = n_pc, basis_Y = basis_Y, Y_new = Y_new,
              trunc_error = trunc_error, mn_Y = mn_Y, sd_Y = sd_Y, ev_Y = ev_Y)

  return(out)
}

pca_reverse <- function(preds_Y_new, pca_Y){
  pca_Y$basis_Y %*% t(preds_Y_new) * pca_Y$sd_Y + pca_Y$mn_Y
}
