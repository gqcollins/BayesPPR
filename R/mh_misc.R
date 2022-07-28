get_move_type <- function(n_ridge, n_quant, n_ridge_max){
  if(n_ridge == 0){
    move_type <- 'birth'
  }else if(n_ridge == n_ridge_max){
    if(n_quant == 0){
      move_type <- 'death'
    }else{
      move_type <- sample(c('death', 'change'), 1)
    }
  }else{
    if(n_quant == 0){
      move_type <- sample(c('birth', 'death'), 1)
    }else{
      move_type <- sample(c('birth', 'death', 'change'), 1)
    }
  }

  return(move_type)
}

get_log_mh_bd <- function(n_ridge_prop, n_quant_prop, n_ridge_max){
  if(n_ridge_prop == 0){
    alpha0 <- 0
  }else if(n_ridge_prop == n_ridge_max){
    if(n_quant_prop == 0){
      alpha0 <- 0
    }else{
      alpha0 <- log(2)
    }
  }else{
    if(n_quant_prop == 0){
      alpha0 <- log(2)
    }else{
      alpha0 <- log(3)
    }
  }

  return(alpha0)
}

myTimestamp <- function(start_time = NULL){
  if(is.null(start_time)){
    x <- Sys.time()
    return(paste('--', format(x, "%b %d %X"), '--'))
  }else{
    x <- Sys.time() - start_time
    return(paste('--', format(round(x, 2)), '--'))
  }
}
