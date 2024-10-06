####################
# Custom Functions #
####################



# 1. Covariance Inversion Function ---------------------------------------------

mat_inv <- function(mat) {
  # mat = valid covariance matrix
  return(chol2inf(chol(mat)))
}

# 2. Informativeness Function --------------------------------------------------


info_value <- function(Sig, Psi, alpha0 = NULL, numthres = 3, dist = c("pois", "binom")) {
  # Sig: spatial covariance
  # Psi <- diag(non-spatial variances)
  
  var <- 1 / diag(solve(Psi + (Psi + Sig) / numthres)) # conditional variance

  if (dist == "pois") {
    info <- 1 / (exp(var) - 1)
    return(as.vector(info))
  } else if (dist == "binom") {
    info <- (1 + exp(alpha0)) / var - exp(alpha0) / (1 + exp(alpha0))
    return(as.vector(info))
  }
}
