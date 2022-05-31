dwallenius <- function(feat_w, feat){ # Wallenius non-central t
  feat_w <- feat_w / sum(feat_w)
  alpha <- feat_w[feat] / sum(feat_w[-feat])
  j <- length(alpha)
  ss <- 1 + (-1)^j * 1/(sum(alpha) + 1)
  if(j > 1){
    for(i in 1:(j - 1))
      ss <- ss + (-1)^(i) * sum(1/(colSums(combn(alpha, i)) + 1))
  }
  ss
}
