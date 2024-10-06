# FIT MCAR MODEL



# setup ------------------------------------------------------------------------

library(here)
library(foreach)
library(doParallel)
load(here("application/data/hd-35-54.rdata"))
load(here("shapefiles/pa.rdata"))
label_sexrace
apply(Y, 2:3, sum)


y_data <- apply(Y, 1:2, sum)
n_data <- apply(n, 1:2, sum)
I <- dim(Y_data)[1]
K <- dim(Y_data)[2]
p <- 1
X <- array(1, dim = c(I, 1))
m <- num

source(here("fun/run_up_mcar_invw.R"))

# MODEL 1 (INV WISHART) --------------------------------------------------------


library(here)
library(MCMCpack)
library(mvtnorm)

set.seed(06272024)
Y_data <- Y[, , s]
n_data <- n[, , s]

# set suffix and path
suff <- paste0("pa_hd_mcar_35-64_", label_sexrace[s])
path <- here("application")
start <- paste0("start_", suff)

# get inits
if (TRUE) {
  get_inits_mcar_invw(
    Y = Y_data, n = n_data, model = "pois", start = start,
    beta0 = log((apply(Y_data, 2, sum) + 0.1) / apply(n_data, 2, sum)),
    sig2 = rep(1, K), rho = rep(0.1, K * (K - 1) / 2), tau2 = rep(0.1, K),
    path = path
  )
}

# run model
if (TRUE) {
  for (i in 1:1000) {
    run_up_mcar_invw(
      start = start, wfile = suff, path = path,
      T_inc = 50, thin = 10
    )
  } # 50,000 iterations

  for (i in 1:5) {
    run_up_mcar_invw(
      start = start, wfile = suff, path = path,
      T_inc = 10000, thin = 10
    )
  } # 100,000 iterations

  run_up_mcar_invw(
    start = start, wfile = suff, path = path,
    T_inc = 100000, thin = 10, get_DIC = TRUE
  )
}
