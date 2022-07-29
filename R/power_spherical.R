rps <- function(mu, kappa){ # Get a random draw from the "power spherical" distribution, which has density f(x) \propto (1 + t(mu) %*% x)^kappa
  d <- length(mu)

  uhat <- -mu
  uhat[1] <- uhat[1] + 1
  u <- uhat / sqrt(sum(uhat^2))

  b <- (d - 1)/2
  a <- b + kappa
  z <- rbeta(1, a, b)
  t <- 2*z - 1

  temp <- rnorm(d - 1)
  v <- temp / sqrt(sum(temp^2))
  y <- c(t, sqrt(1 - t^2) * v)

  uy <- sum(u * y)
  x <- y - 2 * u * uy
  return(matrix(x))
}
