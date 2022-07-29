relu <- function(x){ # Get relu of input
  (abs(x) + x) / 2
}

get_mns_basis <- function(u, knots){ # Get modified natural spline basis
  n_knots <- length(knots)
  df <- n_knots - 2

  basis <- relu(u - knots[1]) # Initialize basis -- would be basis <- u for ns

  if(df > 1){
    n_internal_knots <- n_knots - 3
    r <- d <- list() # Temporary lists for computing basis
    for(k in 1:(n_knots - 1)){
      r[[k]] <- relu(u - knots[k + 1])^3
    }
    for(k in 1:df){
      d[[k]] <- (r[[k]] - r[[n_knots - 1]]) / (knots[n_knots] - knots[k + 1])
    }
    for(k in 1:n_internal_knots){
      basis <- cbind(basis, d[[k]] - d[[df]])
    }
  }

  basis
}

get_cat_basis <- function(X){
  if(ncol(X) > 1){
    tp <- (1 - X[, 1])
    for(j in 2:ncol(X)){
      tp <- tp * (1 - X[, j])
    }
    basis <- matrix(1 - tp)
  }else{
    basis <- X
  }
  basis
}

