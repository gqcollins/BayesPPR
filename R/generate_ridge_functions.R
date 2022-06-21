relu <- function(x){ # Get relu of input
  (abs(x) + x) / 2
}

get_ns_basis <- function(u, knots){ # Get natural spline basis
  n_knots <- length(knots)
  df <- n_knots - 1

  basis <- u # Initialize basis

  if(df > 1){
    n_internal_knots <- n_knots - 2
    r <- d <- list() # Temporary lists for computing basis
    for(k in 1:n_knots){
      r[[k]] <- relu(u - knots[k])^3
    }
    for(k in 1:df){
      d[[k]] <- (r[[k]] - r[[n_knots]]) / (knots[n_knots] - knots[k])
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

