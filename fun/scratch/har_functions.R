###################################################################################
# Functions that utilize the Heterogeneous AR(1) structure of a covariance matrix #
###################################################################################



# 2. har_cov_mat ----------------------------------------------------------------

# function to produce an har structured covariance matrix

har_cov_mat <- function(sig2, rho) {
  n <- length(sig2)
  S <- diag(sqrt(sig2)) # sd matrix
  R <- toeplitz(rho^(0:(n-1))) # correlation matrix
  G <- S %*% R %*% S # covariance matrix

  return(G)
}


# 3. har_cov_inv ----------------------------------------------------------------

# function to produce an inverse of an har structured covariance matrix
# inv1 uses the sig2 and rho vectors
# inv2 uses the covariance matrix G

har_cov_inv1 <- function(sig2, rho) {
  # get G using the sig2 and rho vectors
  G <- har_cov_mat(sig2, rho)
  dim <- dim(G)[1]

  # get the inverse of G
  G_inv <- matrix(0, nrow = dim, ncol = dim)

  # i = j = 1
  G_inv[1, 1] <- 1 / (sig2[1] * (1 - rho^2))

  # if 1 < i = j < Ng
  for (k in 2:(dim - 1)) {
    G_inv[k, k] <- (1 + rho^2) / (sig2[k] * (1 - rho^2))
  }

  # i = j = Ng
  G_inv[dim, dim] <- 1 / (sig2[dim] * (1 - rho^2))

  # |i - j| = 1
  for (i in 1:(dim - 1)) {
    j <- i + 1
    G_inv[i, j] <- G_inv[j, i] <- -rho / (sqrt(sig2[i] * sig2[j]) * (1 - rho^2))
  }

  return(G_inv)
}


har_cov_inv2 <- function(G) {
  dim <- dim(G)[1]

  # get sig2 and rho
  sig2 <- diag(G)
  rho <- G[1, 2] / sqrt(sig2[1] * sig2[2])


  # get the inverse of G
  G_inv <- matrix(0, nrow = dim, ncol = dim)

  # i = j = 1
  G_inv[1, 1] <- 1 / (sig2[1] * (1 - rho^2))
  
  # if 1 < i = j < Ng
  for (k in 2:(dim - 1)) {
    G_inv[k, k] <- (1 + rho^2) / (sig2[k] * (1 - rho^2))
  }
  
  # i = j = Ng
  G_inv[dim, dim] <- 1 / (sig2[dim] * (1 - rho^2))
  
  # |i - j| = 1
  for (i in 1:(dim - 1)) {
    j <- i + 1
    G_inv[i, j] <- G_inv[j, i] <- -rho / (sqrt(sig2[i] * sig2[j]) * (1 - rho^2))
  }

  return(G_inv)
}


