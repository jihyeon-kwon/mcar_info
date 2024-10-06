#######################################
# Extract Correlation from Covariance #
#######################################

extract_corr_mat <- function(G) {
  dim <- dim(G)[1]
  C <- matrix(NA, ncol = dim, nrow = dim)
  S_inv <- diag(1 / sqrt(diag(G))) # inverse of sd matrix
  C <- S_inv %*% G %*% S_inv
  return(C)
}

extract_corr <- function(G) {
  dim <- dim(G)[1]
  r <- numeric((dim - 1) * dim / 2)
  C <- extract_corr_mat(G)
  r <- C[lower.tri(C)]
  return(r)
}

get_corr_lower_mat <- function(rho, dim = dim) {
  C <- diag(1, dim)
  C[lower.tri(C)] <- rho
  return(C)
 }
# 
# get_cov_mat <- function(sig2, rho) {
#   dim <- length(sig2)
#   C <- get_corr_mat(rho, dim = dim)
#   S <- diag(sqrt(sig2)) # inverse of sd matrix
#   G <- S %*% C %*% S
#   return(G)
# }


# cov2cor(G)
# extract_corr_mat(G)
# 
# benchmark_result <- benchmark(
#   cov2cor(G),
#   extract_corr_mat(G),
#   replications = 10000,  # Number of times to repeat the test
#   columns = c("test", "replications", "elapsed", "relative", "user.self", "sys.self")
# )
# print(benchmark_result)