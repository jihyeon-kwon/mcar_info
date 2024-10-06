######################################
# Function to create simulation data #
######################################


# function to create MCAR random effect ----------------------------------------

create_mcar_random <- function(
    file = NULL,
    Ns = Ns, Nt = Nt, n_reps = n_reps, # basic info
    D = D, W = W, G_true = G_true # spatial info
    ) {
  if (!require(mvtnorm)) {
    warning("mvtnorm package is not available")
  }

  # create MCAR random effect
  Q <- D - W
  E <- eigen(Q)

  # create u
  u <- replicate(n_reps, rmvnorm(Ns - 1, rep(0, Nt), G_true), simplify = "array")

  # check
  if (FALSE) {
    dim(u)
    round(cov(t(u[5, , ])), digits = 3)
    round(G_true, digits = 3)
  }

  # create z
  constant <- E$vectors[, 1:(Ns - 1)] %*% diag(E$values[1:(Ns - 1)]^(-1 / 2))
  z <- array(dim = c(Ns, Nt, n_reps))
  for (j in 1:Nt) {
    z[, j, ] <- constant %*% u[, j, ]
  }
  if (!is.null(file)) {
    save(z, D, W, G_true, file = file)
  }
  # add z to the global environment
  assign("z", z, envir = .GlobalEnv)
  return(z)
}


# create mcar data -------------------------------------------------------------

# this one needs an input of data generated from create_mcar_random function.

create_mcar_data <- function(
    file = file,
    n = n, beta0_true = beta0_true, tau2_true = tau2_true,
    model = c("pois", "binom"),
    z = z # random effect
    ) {
  Ns <- dim(z)[[1]]
  Nt <- dim(z)[[2]]
  n_reps <- dim(z)[[3]]

  # check dimension of n and beta0
  if (any(
    !identical(dim(n), c(Ns, Nt)),
    length(beta0_true) != Nt,
    length(tau2_true) != Nt
  )) {
    cat(paste0(
      "n must be ", Ns, " by ", Nt, " matrix \n",
      "beta0 and tau2_true must be length ", Nt
    ))
    return()
  }

  # empty array
  theta <- lambda <- array(dim = c(Ns, Nt, n_reps)) 
  Y <- array(dim = c(Ns, Nt, n_reps))
  
  if (model == "pois") {
    for (r in 1:n_reps) {
      for (i in 1:Ns) {
        for (j in 1:Nt) {
          theta[i, j, r] <- rnorm(1, beta0_true[j] + z[i, j, r], sqrt(tau2_true[j]))
          lambda[i, j, r] <- exp(theta[i, j, r])
          Y[i, j, r] <- rpois(1, n[i, j] * lambda[i, j, r])
        }
      }
    }
    
  } else if (model == "binom") {
    for (r in 1:n_reps) {
      for (i in 1:Ns) {
        for (j in 1:Nt) {
          theta[i, j, r] <- rnorm(1, beta0_true[j] + z[i, j, r], sqrt(tau2_true[j]))
          lambda[i, j, r] <- exp(theta[i, j, r]) / (1 + exp(theta[i, j, r]))
          Y[i, j, r] <- rbinom(1, n[i, j], lambda[i, j, r])
        }
      }
    }
  }

  # check
  if (FALSE) {
    matplot(t(apply(theta, 1:2, mean)), type = "l")
    beta0_true
    matplot(t(apply(theta - beta0_true - z, 1:2, sd)), type = "l")
    sqrt(tau2_true)
  }

  if (!is.null(file)) {
    save(Y, n, beta0_true, tau2_true, z, file = file)
  }
  # add Y and n to the global environment
  assign("Y", Y, envir = .GlobalEnv)
  assign("n", n, envir = .GlobalEnv)
  return(c(Y, n))
}
