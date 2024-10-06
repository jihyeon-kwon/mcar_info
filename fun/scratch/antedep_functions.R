#######################################################################
# Functions that utilize the antedep structure of a covariance matrix #
#######################################################################


# 1. antedep_corr_mat ---------------------------------------------------------------

# function to produce an antedep structured correlation matrix
# when rho's for each group is given

antedep_corr_mat <- function(rho) {
  dim <- length(rho)
  total_length <- (dim + 1) * dim / 2
  R <- matrix(1, nrow = dim + 1, ncol = dim + 1)
  upper <- numeric(total_length)
  lower <- numeric(total_length)

  index <- 1
  for (i in 1:dim) {
    for (j in 1:i) {
      upper[index] <- prod(rho[i:j])
      index <- index + 1
    }
  }

  index <- 1
  for (i in 1:dim) {
    for (j in i:dim) {
      lower[index] <- prod(rho[i:j])
      index <- index + 1
    }
  }

  R[upper.tri(R)] <- upper
  R[lower.tri(R)] <- lower
  return(R)
}


# 2. antedep_cov_mat ----------------------------------------------------------------

# function to produce an antedep structured covariance matrix

antedep_cov_mat <- function(sig2, rho) {
  dim <- length(sig2)
  S <- diag(sqrt(sig2)) # sd matrix
  R <- antedep_corr_mat(rho) # correlation matrix
  G <- S %*% R %*% S # covariance matrix

  return(G)
}


# 3. antedep_cov_inv ----------------------------------------------------------------

# function to produce an inverse of an antedep structured covariance matrix
# inv1 uses the sig2 and rho vectors
# inv2 uses the covariance matrix G

antedep_cov_inv1 <- function(sig2, rho) {
  # get G using the sig2 and rho vectors
  G <- antedep_cov_mat(sig2, rho)
  dim <- dim(G)[1]

  # get the inverse of G
  G_inv <- matrix(0, nrow = dim, ncol = dim)

  # i = j = 1
  G_inv[1, 1] <- 1 / (sig2[1] * (1 - rho[1]^2))

  # if 1 < i = j < Ng
  for (k in 2:(dim - 1)) {
    G_inv[k, k] <- (1 - prod(rho[(k - 1):k])^2) / (sig2[k] * (1 - rho[k - 1]^2) * (1 - rho[k]^2))
  }

  # i = j = Ng
  G_inv[dim, dim] <- 1 / (sig2[dim] * (1 - rho[dim - 1]^2))

  # |i - j| = 1
  for (i in 1:(dim - 1)) {
    j <- i + 1
    G_inv[i, j] <- G_inv[j, i] <- -rho[i] / (sqrt(sig2[i] * sig2[j]) * (1 - rho[i]^2))
  }

  return(G_inv)
}

antedep_cov_inv2 <- function(G) {
  dim <- dim(G)[1]

  # get sig2 and rho
  sig2 <- diag(G)
  rho <- numeric(dim - 1)
  for (j in 1:(dim - 1)) {
    rho[j] <- G[j, (j + 1)] / sqrt(prod(sig2[j:(j + 1)]))
  }

  # get the inverse of G
  G_inv <- matrix(0, nrow = dim, ncol = dim)

  # i = j = 1
  G_inv[1, 1] <- 1 / (sig2[1] * (1 - rho[1]^2))

  # if 1 < i = j < Ng
  for (k in 2:(dim - 1)) {
    G_inv[k, k] <- (1 - prod(rho[(k - 1):k])^2) / (sig2[k] * (1 - rho[k - 1]^2) * (1 - rho[k]^2))
  }

  # i = j = Ng
  G_inv[dim, dim] <- 1 / (sig2[dim] * (1 - rho[dim - 1]^2))

  # |i - j| = 1
  for (i in 1:(dim - 1)) {
    j <- i + 1
    G_inv[i, j] <- G_inv[j, i] <- -rho[i] / (sqrt(sig2[i] * sig2[j]) * (1 - rho[i]^2))
  }

  return(G_inv)
}

