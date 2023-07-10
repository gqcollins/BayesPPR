f <- function(xx)
{
  ##########################################################################
  #
  # WELCH ET AL. (1992) FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it
  # and/or modify it under the terms of the GNU General Public License as
  # published by the Free Software Foundation; version 2.0 of the License.
  # Accordingly, this program is distributed in the hope that it will be
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUT:
  #
  # xx = c(x1, x2, ..., 20)
  #
  ##########################################################################

  x1  <- xx[1] - 0.5
  x2  <- xx[2] - 0.5
  x3  <- xx[3] - 0.5
  x4  <- xx[4] - 0.5
  x5  <- xx[5] - 0.5
  x6  <- xx[6] - 0.5
  x7  <- xx[7] - 0.5
  x8  <- xx[8] - 0.5
  x9  <- xx[9] - 0.5
  x10 <- xx[10] - 0.5
  x11 <- xx[11] - 0.5
  x12 <- xx[12] - 0.5
  x13 <- xx[13] - 0.5
  x14 <- xx[14] - 0.5
  x15 <- xx[15] - 0.5
  x16 <- xx[16] - 0.5
  x17 <- xx[17] - 0.5
  x18 <- xx[18] - 0.5
  x19 <- xx[19] - 0.5
  x20 <- xx[20] - 0.5

  term1 <- 5*x12 / (1+x1)
  term2 <- 5 * (x4-x20)^2
  term3 <- x5 + 40*x19^3 - 5*x19
  term4 <- 0.05*x2 + 0.08*x3 - 0.03*x6
  term5 <- 0.03*x7 - 0.09*x9 - 0.01*x10
  term6 <- -0.07*x11 + 0.25*x13^2 - 0.04*x14
  term7 <- 0.06*x15 - 0.01*x17 - 0.03*x18

  y <- term1 + term2 + term3 + term4 + term5 + term6 + term7
  return(y)
}


set.seed(71023)
n <- 5000
p <- 20
q <- 2 # 2 categorical inputs
X <- cbind(matrix(runif(n*p), n), sample(0:1, n, replace = TRUE), sample(0:1, n, replace = TRUE))
y <- apply(X, 1, f) + X[, 21] - X[, 22] + 0.5 * X[, 21]*X[, 22] + rnorm(n, sd = 0.2)

fit <- bppr(X, y)
plot(fit$sd_resid, type = 'l')

X_test <- cbind(matrix(runif(n*p), n), sample(0:1, n, replace = TRUE), sample(0:1, n, replace = TRUE))
preds <- predict(fit, X_test)
pred_mn <- apply(preds, 2, mean)
plot(apply(X_test, 1, f), pred_mn)
abline(0, 1, col = 2)
sqrt(mean((apply(X_test, 1, f) - pred_mn)^2))


