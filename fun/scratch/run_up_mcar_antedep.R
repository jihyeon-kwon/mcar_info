###############
# Run Up MCAR #
###############

# MCAR Model with Inverse Wishart Prior
# Covariance matrix is not structured

# Last Updated: Jun 25, 2024


# setup ------------------------------------------------------------------------
library(MCMCpack)
library(mvtnorm)
library(mcmc)


# get initial values -----------------------------------------------------------

get_inits_mcar_antedep <- function(
    Y = Y,
    n = n,
    model = c("pois", "binom"),
    start = "result",
    beta0 = numeric(Nt), # baseline
    z = array(0, dim = c(Ns, Nt)), # MCAR spatial random effect
    sig2 = numeric(Nt), # vector of spatial variances
    rho = numeric(Nt - 1), # vector of spatial correlations
    tau2 = numeric(Nt), # vector of non-spatial variances,
    theta = NULL, # vector of either log(prevalence) when pois or logit(incidence) when binom
    qt = array(1, dim = c(Ns, Nt)), # step sizes for theta
    t_accept = array(.4, dim = c(Ns, Nt)), # acceptance rate for theta
    qs = rep(5, Nt), # inverse stepsize for sig2
    s_accept = rep(.4, Nt), # acceptance rate for sig2
    qr = rep(5, Nt - 1), # inverse stepsize for rho
    r_accept = rep(.4, Nt - 1), # acceptance rate for rho,
    Tn = 1, # iteration
    path = "",
    prev = NA, # batch number of previous run
    numthres = 3, # number of common neighbors
    thres = Inf,
    info = NULL) { # informativeness

  if (length(rho) != (Nt - 1)) {
    stop(paste0("rho must be of length ", Nt-1))
  }

  if (is.na(prev)) { # do we have any values from previous runs?
    if (is.null(beta0) | is.null(sig2) | is.null(tau2)) {
      cat("You need to specify beta0, sig2, and tau2 \n")
      return()
    }

    if (!require(here)) {
      warning("here package is not available")
    } else {
      source(here("fun/info_functions.R"))
      source(here("fun/antedep_functions.R"))
    }

    Ns <- dim(Y)[1]
    Nt <- dim(Y)[2]

    # theta = log(lambda) or logit(pi)
    if (is.null(theta)) {
      if (model == "pois") {
        theta <- ifelse(
          Y != 0 & n != 0 & Y != n,
          log(Y / n),
          matrix(rep(beta0, each = Ns), nrow = Ns)
        )
        lam <- exp(theta)
      } else if (model == "binom") {
        theta <- ifelse(
          Y != 0 & n != 0 & Y != n,
          qlogis(Y / n),
          matrix(rep(beta0, each = Ns), nrow = Ns)
        )
        lam <- exp(theta) / (1 + exp(theta))
      }
    }

    # spatial random effect from MCAR
    if (is.null(z)) {
      z <- theta - matrix(rep(beta0, each = Ns), nrow = Ns)
    }

    # informativenss
    G <- antedep_cov_mat(sig2 = sig2, rho = rho)
    R <- diag(tau2)

    if (is.null(info)) {
      info <- info_value(G = G, R = R, numthres = numthres, dist = model)
    }

    # save it to mod
    mod <- list(
      Y = Y, n = n,
      beta0 = beta0, theta = theta, z = z, lam = lam,
      sig2 = sig2, rho = rho, tau2 = tau2,
      qt = qt, t_accept = t_accept,
      qs = qs, s_accept = s_accept,
      qr = qr, r_accept = r_accept,
      Tn = Tn, total = 1,
      thres = thres, numthres = numthres, info = info,
      seed = .Random.seed,
      time = 0
    )
  } else { # we take the previous run's results as the inits

    cat("loading results from previous run")

    load(
      file = paste0(path, "/output/", prev, ".rdata"),
      envir = sys.frame(sys.nframe())
    )
    mod$seed <- .Random.seed # keep the same see number so that it will work as one chain
  }

  save(mod,
    file = paste0(path, "/output/", start, ".rdata")
  )

  cat("here\n")
  cat(dim(mod$theta), "\n")
  cat("inits saved\n")
  rm(mod)
}


# run up mcar ------------------------------------------------------------------

run_up_mcar_antedep <- function(
    start = "",
    T_inc = 10, # number of iterations in this batch
    wfile = "current",
    path = "",
    a_s = 1, b_s = 1 / 7, # hyperparam for sig2
    a_t = 1, b_t = 1 / 100, # hyperparam for tau2
    a_r = 1.1, b_r = 1.1, # hyperparam for beta
    thin = 1,
    numthres = 3, # number of common neighbors
    get_DIC = FALSE) {
  if (!require(here)) {
    warning("here package is not available")
  } else {
    source(here("fun/antedep_functions.R"))
    source(here("fun/info_functions.R"))
  }

  load(
    file = paste0(path, "/output/", start, ".rdata"),
    envir = sys.frame(which = sys.nframe())
  )

  # read in the last values
  Tn <- mod$Tn + T_inc
  oldT <- mod$Tn # this should always be 1

  ## values that we save (copy from the previous batch)
  lam <- array(dim = c(Ns, Nt, Tn))
  lam[, , 1:mod$Tn] <- mod$lam # lam = exp(theta) for poisson
  theta <- mod$theta
  beta0 <- array(dim = c(Nt, Tn))
  beta0[, 1:mod$Tn] <- mod$beta0
  sig2 <- array(dim = c(Nt, Tn))
  sig2[, 1:mod$Tn] <- mod$sig2
  rho <- array(dim = c(Nt - 1, Tn))
  rho[, 1:mod$Tn] <- mod$rho
  tau2 <- array(dim = c(Nt, Tn))
  tau2[, 1:mod$Tn] <- mod$tau2
  info <- array(dim = c(Nt, Tn))
  info[, 1:mod$Tn] <- mod$info
  Y <- mod$Y
  n <- mod$n

  ## dimnames
  dimnames(lam) <- dimnames(Y)
  dimnames(beta0) <- dimnames(sig2) <- dimnames(tau2) <- dimnames(Y)[2]

  ## values that we don't save
  qt <- mod$qt
  qs <- mod$qs
  qr <- mod$qr
  z <- mod$z
  theta <- mod$theta
  G <- antedep_cov_mat(sig2 = sig2[, mod$Tn], rho = rho[, mod$Tn])
  G_inv <- antedep_cov_inv2(G)

  numthres <- mod$numthres
  thres <- mod$thres
  total <- mod$total # total number of iterations
  .Random.Seed <- mod$seed
  time <- mod$time
  if (get_DIC) {
    dtheta <- rep(NA, T_inc)
  } # for DIC

  # RUN THE MODEL --------------------------------------------------------------

  # tuning the proposal densities
  t_accept <- mod$t_accept
  t_accept <- ifelse(t_accept < 1 / 6, 1 / 6,
    ifelse(t_accept > 0.75, 0.75, t_accept)
  )
  qt <- qt * t_accept / 0.44
  t_accept <- array(0, dim = c(Ns, Nt))

  s_accept <- mod$s_accept
  s_accept <- ifelse(s_accept < 1 / 6, 1 / 6,
    ifelse(s_accept > 0.75, 0.75, s_accept)
  )
  qs <- qs * 0.44 / s_accept
  qs <- ifelse(qs > 1, qs, 1)
  qs <- ifelse(qs > 1000, 1000, qs)
  s_accept <- rep(0, Nt)

  r_accept <- mod$r_accept
  r_accept <- ifelse(r_accept < 1 / 6, 1 / 6,
    ifelse(r_accept > 0.75, 0.75, r_accept)
  )
  qr <- qr * 0.44 / r_accept
  qr <- ifelse(qr > 1, qr, 1)
  qr <- ifelse(qr > 1000, 1000, qr)
  r_accept <- rep(0, Nt - 1)

  time_start <- Sys.time()

  for (it in (oldT + 1):Tn) {
    # update beta0
    beta_var <- tau2[, it - 1] / Ns
    beta_mean <- apply(theta[, ] - z[, ], 2, sum) / Ns
    beta0[, it] <- rnorm(Nt, beta_mean, sd = sqrt(beta_var))

    # update z
    G_inv <- antedep_cov_inv1(sig2 = sig2[, it - 1], rho = rho[, it - 1])
    for (i in 1:Ns) {
      if (num[i] > 1) {
        muia <- apply(z[neigh[[i]], ], 2, mean)
      } else {
        muia <- z[neigh[[i]], ]
      }
      G_invi <- m[i] * G_inv

      z_prec <- diag(1 / tau2[, it - 1]) + G_invi
      z_var <- chol(z_prec)
      z_var <- chol2inv(z_var)
      z_mean <- z_var %*% (diag(1 / tau2[, it - 1]) %*%
        c(theta[i, ] - X[i, ] %*% t(beta0[, it])) +
        G_invi %*% muia)
      z[i, ] <- rmvnorm(1, z_mean, z_var)
    }
    z[, ] <- z[, ] - rep(apply(z[, ], 2, mean), each = Ns)

    # update theta (lam)
    ts <- rnorm(Ns * Nt, theta[, ], qt) # candidate
    ra <- Y[, ] * (ts - theta[, ])
    rb <- n * (exp(ts) - exp(theta[, ]))
    rc <- ((ts - X %*% t(beta0[, it]) - z[, ])^2 -
      (theta[, ] - X %*% t(beta0[, it]) - z[, ])^2)
    r <- exp(ra - rb - rc / (2 * rep(tau2[, it - 1], each = Ns)))
    accept <- (r > runif(Ns * Nt))
    t_accept <- t_accept + accept
    theta[, ] <- ifelse(accept, ts, theta[, ])

    # get lam = exp(theta)
    lam[, , it] <- exp(theta[, ])

    # get DIC
    if (get_DIC) {
      dtheta[it - oldT] <- -2 * sum(Y[, ] * theta[, ] - n * lam[, , it])
    }

    # update tau2
    tau2[, it] <- 1 / rgamma(
      Nt,
      Ns / 2 + a_t,
      apply((theta[, ] - X %*% t(beta0[, it]) - z[, ])^2, 2, sum) / 2 + b_t
    )

    # update sig2
    ## candidates
    sig2_star <- rgamma(Nt, qs, qs / sig2[, it - 1])
    rho_star <- rbeta(Nt - 1, qr * rho[it - 1] + 1, qr * (1 - rho[it - 1]) + 1)
    G_star <- antedep_cov_mat(sig2 = sig2_star, rho = rho_star)
    G_star_inv <- antedep_cov_inv2(G_star)

    ## assess ratio
    # sum <- 0
    # sum_star <- 0
    diff <- 0
    for (i in 1:Ns) {
      if (m[i] == 1) {
        zneigh <- z[neigh[[i]], ]
      } else {
        zneigh <- apply(z[neigh[[i]], ], 2, mean)
      }
      diff <- diff + t(z[i, ]) %*% (G_star_inv - G_inv) %*% (m[i] * (z[i, ] - zneigh))
    }

    R_mcar <- -(Ns - 1) / 2 * (log(det(G_star)) - log(det(G))) - 1 / 2 * diff
    R_ig <- sum((-a_s - 1) * (log(sig2_star) - log(sig2[, it - 1])) + (-b_s * (1 / sig2_star - 1 / sig2[, it - 1])))
    R_gam <- sum(dgamma(sig2[, it - 1], qs, qs / sig2[, it - 1], log = TRUE) - dgamma(sig2_star, qs, qs / sig2_star, log = TRUE))
    R_beta1 <- sum(dbeta(rho_star, a_r, b_r, log = TRUE) - dbeta(rho[, it - 1], a_r, b_r, log = TRUE))
    R_beta2 <- sum(dbeta(rho[, it - 1], qr * rho_star + 1, qr * (1 - rho_star) + 1, log = TRUE) -
                     dbeta(rho_star, qr * rho[, it - 1] + 1, qr * (1 - rho[, it - 1]) + 1, log = TRUE))
    r <- exp(R_mcar + R_ig + R_gam + R_beta1 + R_beta2)

    ## accept/reject
    accept <- (r > runif(1))
    if (accept) {
      sig2[, it] <- sig2_star
      rho[, it] <- rho_star
      G <- G_star
      G_inv <- G_star_inv
    } else {
      sig2[, it] <- sig2[, it - 1]
      rho[, it] <- rho[, it - 1]
      G <- antedep_cov_mat(sig2 = sig2[, it - 1], rho = rho[, it - 1])
      G_inv <- antedep_cov_inv2(G)
    }
    s_accept <- s_accept + rep(accept, Nt)
    r_accept <- r_accept + rep(accept, Nt-1)
    
    # compute the informativeness
    info[, it] <- info_value(
      G = G, R = diag(tau2[, it]),
      dist = "pois", numthres = numthres
    )

    # print the iteration
    if ((it + total) / 10 == floor((it + total) / 10)) {
      cat("Finished iteration:", (it + total), "info = ", (info[1, it]), "\n")
    }
  }

  time_end <- Sys.time()
  current_time <- as.numeric(difftime(time_end, time_start, units = "mins"))
  time <- time + current_time

  cat(
    "** Total Run Time:", time, "mins ",
    "Current Run Time:", current_time, "mins \n"
  )

  # save the results
  total <- total + T_inc
  t_accept <- t_accept / Tn
  s_accept <- s_accept / Tn
  r_accept <- r_accept / Tn

  # compute the DIC
  if (get_DIC) {
    Dbar <- mean(dtheta)
    lbar <- apply(lam, 1:2, mean)
    tbar <- log(lbar)
    Dhat <- -2 * sum(tbar * Y - n * lbar)
    pD <- Dbar - Dhat
    DIC <- Dbar + pD
  }

  # save the last values (initial values for the next batch)
  mod <- list(
    Y = Y, n = n,
    beta0 = beta0[, Tn], theta = theta, lam = lam[, , Tn],
    z = z, sig2 = sig2[, Tn], tau2 = tau2[, Tn], rho = rho[, Tn],
    info = info[, Tn], numthres = numthres, thres = thres,
    seed = .Random.seed,
    t_accept = t_accept, qt = qt,
    s_accept = s_accept, qs = qs,
    r_accept = r_accept, qr = qr,
    Tn = 1, total = total, time = time
  )

  save(mod,
    file = paste(path, "/output/", start, ".rdata", sep = ""),
    envir = sys.frame(which = sys.nframe())
  )
  rm(mod)

  # save samples from this batch (thinning)
  indices <- oldT + seq(thin, T_inc, thin)
  mod <- list(
    Y = Y, n = n,
    beta0 = beta0[, indices], theta = theta[, ], lam = lam[, , indices],
    z = z, sig2 = sig2[, indices], tau2 = tau2[, indices], rho = rho[, indices],
    info = info[, indices], numthres = numthres, thres = thres,
    seed = .Random.seed,
    t_accept = t_accept, qt = qt,
    s_accept = s_accept, qs = qs,
    r_accept = r_accept, qr = qr,
    Tn = T_inc, total = total, time = time
  )
  if (get_DIC) {
    mod$DIC <- DIC
  }

  calc_sum <- function(x) {
    c(quantile(x, probs = c(0.025, 0.5, 0.975)), mean(x), sd(x))
  }

  summary <- list(
    beta0 = apply(beta0[, indices], 1, calc_sum),
    lam = apply(lam[, , indices], 1:2, calc_sum),
    sig2 = apply(sig2[, indices], 1, calc_sum),
    rho = apply(rho[, indices], 1, calc_sum),
    tau2 = apply(tau2[, indices], 1, calc_sum),
    info = apply(info[, indices], 1, calc_sum)
  )

  rownames(summary$beta0) <- rownames(summary$sig2) <-
    rownames(summary$tau2) <- rownames(summary$rho) <-
    rownames(summary$info) <- dimnames(summary$lam)[[1]] <-
    c("2.5%", "50%", "97.5%", "mean", "sd")

  save(mod, summary,
    file = paste(path, "/output/", wfile, ".rdata", sep = ""),
    envir = sys.frame(which = sys.nframe())
  )
  return(summary)
}
