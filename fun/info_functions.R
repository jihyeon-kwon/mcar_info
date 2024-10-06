################################
# Function for Informativeness #
################################

# Sig: spatial covariance
# Psi <- diag(non-spatial variances)

info_value <- function(Sig, Psi, alpha0 = NULL, numthres = 3, dist = c("pois", "binom")) {
  n <- dim(Sig)[1]
  var <- 1 / diag(solve(Psi + (Psi + Sig) / numthres)) # conditional variance
  
  if (dist == "pois") {
    info <- 1 / (exp(var) - 1)
    return(as.vector(info))
    
  } else if (dist == "binom") {
    info <- (1 + exp(alpha0)) / var - exp(alpha0) / (1 + exp(alpha0))
    return(as.vector(info))
    
  }
}


