###############
# Run Up MCAR #
###############

# MCAR Model with Inverse Wishart Prior

# Last Updated: Aug 14


# setup ------------------------------------------------------------------------
library(MCMCpack)
library(mvtnorm)
library(mcmc)
library(MBESS)


# get initial values -----------------------------------------------------------

get_inits_mcar_invw <- function(
    Y = Y,
    n = n,
    model = c("pois", "binom"),
    start = "result",
    alpha0 = numeric(K), # baseline
    z = array(0, dim = c(I, K)), # MCAR spatial random effect
    sig2 = numeric(K), # vector of spatial variances
    rho = numeric(K * (K - 1) / 2), # vector of spatial correlations
    tau2 = numeric(K), # vector of non-spatial variances,
    theta = NULL, # vector of either log(prevalence) when pois or logit(incidence) when binom
    qt = array(1, dim = c(I, K)), # step sizes for theta
    t_accept = array(.4, dim = c(I, K)), # acceptance rate for theta
    Tn = 1, # iteration
    path = "",
    prev = NA, # batch number of previous run
    numthres = 3, # number of common neighbors
    thres = Inf,
    info = NULL) { # informativeness

  if (is.na(prev)) { # do we have any values from previous runs?
    if (is.null(alpha0) | is.null(sig2) | is.null(tau2)) {
      cat("You need to specify alpha0, sig2, and tau2 \n")
      return()
    }

    if (!require(here)) {
      warning("here package is not available")
    } else {
      source(here("fun/custom_functions.r"))
    }

    I <- dim(Y)[1]
    K <- dim(Y)[2]

    # theta = log(lambda) or logit(pi)
    if (is.null(theta)) {
      if (model == "pois") {
        theta <- ifelse(
          Y != 0 & n != 0 & Y != n,
          log(Y / n),
          matrix(rep(alpha0, each = I), nrow = I)
        )
        lam <- exp(theta)
      } else if (model == "binom") {
        theta <- ifelse(
          Y != 0 & n != 0 & Y != n,
          qlogis(Y / n),
          matrix(rep(alpha0, each = I), nrow = I)
        )
        lam <- exp(theta) / (1 + exp(theta))
      }
    }

    # spatial random effect from MCAR
    if (is.null(z)) {
      z <- theta - matrix(rep(alpha0, each = I), nrow = I)
    }

    # informativenss
    C <- diag(1, K)
    C[lower.tri(C)] <- rho
    C <- C + t(C) - diag(1, K)
    Sig <- MBESS::cor2cov(C, sqrt(sig2))
    Psi <- diag(tau2)

    if (is.null(info)) {
      info <- info_value(
        Sig = Sig, Psi = Psi, alpha0 = alpha0, numthres = numthres, dist = model
      )
    }

    # save it to mod
    mod <- list(
      Y = Y, n = n,
      alpha0 = alpha0, theta = theta, z = z, lam = lam,
      sig2 = sig2, rho = rho, tau2 = tau2,
      qt = qt, t_accept = t_accept,
      Tn = Tn, total = 1,
      thres = thres, numthres = numthres, info = info,
      model = model,
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


# run up mcar with inverse wishart prior ---------------------------------------

run_up_mcar_invw <- function(
    start = "",
    T_inc = 10, # number of iterations in this batch
    wfile = "current",
    path = "",
    thin = 1,
    a_t = 1, b_t = 1 / 100, 
    nu = K + 1, G0 = diag(1 / 7, K), # hyperparams
    numthres = 3, # number of common neighbors
    get_DIC = FALSE) {
  if (!require(here)) {
    warning("here package is not available")
  } else {
    source(here("fun/custom_functions.r"))
  }

  load(
    file = paste0(path, "/output/", start, ".rdata"),
    envir = sys.frame(which = sys.nframe())
  )

  # read in the last values
  Tn <- mod$Tn + T_inc
  oldT <- mod$Tn # this should always be 1

  ## values that we save (copy from the previous batch)
  lam <- array(dim = c(I, K, Tn))
  lam[, , 1:mod$Tn] <- mod$lam # lam = exp(theta) for pois
  theta <- mod$theta
  alpha0 <- array(dim = c(K, Tn))
  alpha0[, 1:mod$Tn] <- mod$alpha0
  sig2 <- array(dim = c(K, Tn))
  sig2[, 1:mod$Tn] <- mod$sig2
  rho <- array(dim = c(K * (K - 1) / 2, Tn))
  rho[, 1:mod$Tn] <- mod$rho
  tau2 <- array(dim = c(K, Tn))
  tau2[, 1:mod$Tn] <- mod$tau2
  info <- array(dim = c(K, Tn))
  info[, 1:mod$Tn] <- mod$info
  Y <- mod$Y
  n <- mod$n

  ## dimnames
  dimnames(lam) <- dimnames(Y)
  dimnames(alpha0) <- dimnames(sig2) <- dimnames(tau2) <- dimnames(Y)[2]

  ## values that we don't save
  qt <- mod$qt
  z <- mod$z
  theta <- mod$theta
  C <- diag(1, K)
  C[lower.tri(C)] <- rho[, mod$Tn]
  C <- C + t(C) - diag(1, K)
  Sig <- MBESS::cor2cov(C, sqrt(sig2[, mod$Tn]))
  Sig_inv <- mat_inv(Sig)

  numthres <- mod$numthres
  thres <- mod$thres
  total <- mod$total # total number of iterations
  model <- mod$model
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
  t_accept <- array(0, dim = c(I, K))

  qs <- qs
  s_accept <- 0.4 # not used in invw... but need it for unstr
  qr <- qr
  r_accept <- 0.4

  time_start <- Sys.time()

  for (it in (oldT + 1):Tn) {
    # update alpha0 -------------------------------------------------------------

    beta_var <- tau2[, it - 1] / I
    beta_mean <- apply(theta[, ] - z[, ], 2, mean)
    alpha0[, it] <- rnorm(K, beta_mean, sd = sqrt(beta_var))

    # eval the info
    if (model == "binom" && thres < Inf) { ## if restricted binom

      info_cand <- info_value(
        Sig = Sig, Psi = diag(tau2[, it - 1]), alpha0 = alpha0[, it],
        dist = model, numthres = numthres
      )

      # accept/reject
      if (prod(info_cand < thres) != 1) {
        alpha0[, it] <- alpha0[, it - 1] # then reject
      }
    }

    # update z -----------------------------------------------------------------

    for (i in 1:I) {
      if (num[i] > 1) {
        muia <- apply(z[neigh[[i]], ], 2, mean)
      } else {
        muia <- z[neigh[[i]], ]
      }

      Sig_invi <- m[i] * Sig_inv
      z_prec <- diag(1 / tau2[, it - 1]) + Sig_invi
      z_var <- mat_inv(z_prec)
      z_mean <- z_var %*% (diag(1 / tau2[, it - 1]) %*%
        c(theta[i, ] - X[i, ] %*% t(alpha0[, it])) +
        Sig_invi %*% muia)
      z[i, ] <- rmvnorm(1, z_mean, z_var)
    }

    z[, ] <- z[, ] - rep(apply(z[, ], 2, mean), each = I)



    # update theta -------------------------------------------------------------

    ts <- rnorm(I * K, theta[, ], qt) # candidate
    ra <- Y[, ] * (ts - theta[, ])
    if (model == "pois") {
      rb <- n * (exp(ts) - exp(theta[, ]))
    } else if (model == "binom") {
      rb <- n * (log(1 + exp(ts)) - log(1 + exp(theta[, ])))
    }
    rc <- ((ts - X %*% t(alpha0[, it]) - z[, ])^2 -
      (theta[, ] - X %*% t(alpha0[, it]) - z[, ])^2)
    r <- exp(ra - rb - rc / (2 * rep(tau2[, it - 1], each = I)))
    accept <- (r > runif(I * K))
    t_accept <- t_accept + accept
    theta[, ] <- ifelse(accept, ts, theta[, ])


    # get lambda (either prevalence for Poisson or rate for Binom) -------------

    if (model == "pois") {
      lam[, , it] <- exp(theta[, ])
    } else if (model == "binom") {
      lam[, , it] <- exp(theta[, ]) / (1 + exp(theta[, ]))
    }

    # get DIC ------------------------------------------------------------------

    if (get_DIC) {
      if (model == "pois") {
        dtheta[it - oldT] <- -2 * sum(Y[, ] * theta[, ] - n * lam[, , it])
      } else if (model == "binom") {
        dtheta[it - oldT] <- -2 * sum(Y[, ] * log(lam[, , it]) +
          (n - Y[, ]) * log(1 - lam[, , it]))
      }
    }


    # update tau2 (non-spatial variance) ---------------------------------------

    tau2[, it] <- 1 / rgamma(
      K,
      I / 2 + a_t,
      apply((theta[, ] - X %*% t(alpha0[, it]) - z[, ])^2, 2, sum) / 2 + b_t
    )

    Psi <- diag(tau2[, it])

    # eval info
    if (thres < Inf) {
      info_cand <- info_value(
        Sig =
          Sig, Psi = Psi, alpha0 = alpha0[, it], dist = model, numthres = numthres
      )
    } else {
      info_cand <- rep(0, K)
    }

    # accept/reject
    if (prod(info_cand < thres) != 1) {
      tau2[, it] <- tau2[, it - 1]
      Psi <- diag(tau2[, it])
    }

    # update Sig (spatial variance) ----------------------------------------------

    Ags <- G0

    for (i in 1:I) {
      if (m[i] == 1) {
        zneigh <- z[neigh[[i]], ]
      } else {
        zneigh <- apply(z[neigh[[i]], ], 2, mean)
      }
      Ags <- Ags + m[i] * (z[i, ] - zneigh) %*% t(z[i, ])
    }

    Sig_star <- riwish(I - 1 + nu, Ags)

    # eval info
    if (thres < Inf) {
      info_cand <- info_value(
        Sig = Sig_star, Psi = Psi, alpha0 = alpha0[, it], dist = model, numthres = numthres
      )
    } else {
      info_cand <- rep(0, K)
    }

    # accept/reject
    if (prod(info_cand < thres) == 1) {
      Sig <- Sig_star
      sig2[, it] <- diag(Sig[, ])
      C <- cov2cor(Sig[, ])
      rho[, it] <- C[lower.tri(C)]
      Sig_inv <- mat_inv(Sig[, ])
    } else {
      sig2[, it] <- sig2[, it - 1]
      rho[, it] <- rho[, it - 1]
    }



    # compute the informativeness ----------------------------------------------

    info[, it] <- info_value(
      Sig = Sig, Psi = Psi, alpha0 = alpha0[, it], dist = model, numthres = numthres
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

  # compute the DIC
  if (get_DIC) {
    Dbar <- mean(dtheta)
    lbar <- apply(lam, 1:2, mean)
    if (model == "pois") {
      tbar <- log(lbar)
      Dhat <- -2 * sum(tbar * Y - n * lbar)
    } else if (model == "binom") {
      tbar <- log(lbar / (1 - lbar))
      Dhat <- -2 * sum(Y * log(lbar) - (Y - n) * log(1 - lbar))
    }

    pD <- Dbar - Dhat
    DIC <- Dbar + pD
  }

  # save the last values (initial values for the next batch)
  mod <- list(
    Y = Y, n = n, model = model,
    alpha0 = alpha0[, Tn], theta = theta, lam = lam[, , Tn],
    z = z, sig2 = sig2[, Tn], tau2 = tau2[, Tn], rho = rho[, Tn],
    info = info[, Tn], numthres = numthres, thres = thres,
    seed = .Random.seed, t_accept = t_accept, qt = qt,
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
    Y = Y, n = n, model = model,
    alpha0 = alpha0[, indices], theta = theta[, ], lam = lam[, , indices],
    z = z, sig2 = sig2[, indices], tau2 = tau2[, indices], rho = rho[, indices],
    info = info[, indices], numthres = numthres, thres = thres,
    seed = .Random.seed, t_accept = t_accept, qt = qt,
    s_accept = s_accept, qs = qs, # not needed for invw but for unstr later
    r_accept = r_accept, qr = qr, # not needed for invw but for unstr later
    Tn = T_inc, total = total, time = time
  )
  if (get_DIC) {
    mod$DIC <- DIC
    mod$pD <- pD
    mod$Dbar <- Dbar
  }

  calc_sum <- function(x) {
    c(quantile(x, probs = c(0.025, 0.5, 0.975)), mean(x), sd(x))
  }

  if (dim(rho)[1] == 1) {
    summary <- list(
      alpha0 = apply(alpha0[, indices], 1, calc_sum),
      lam = apply(lam[, , indices], 1:2, calc_sum),
      sig2 = apply(sig2[, indices], 1, calc_sum),
      rho = calc_sum(rho[, indices]),
      tau2 = apply(tau2[, indices], 1, calc_sum),
      info = apply(info[, indices], 1, calc_sum)
    )
    
    rownames(summary$alpha0) <- rownames(summary$sig2) <-
      rownames(summary$tau2) <- names(summary$rho) <-
      rownames(summary$info) <- dimnames(summary$lam)[[1]] <-
      c("2.5%", "50%", "97.5%", "mean", "sd")
    
  } else {
    summary <- list(
      alpha0 = apply(alpha0[, indices], 1, calc_sum),
      lam = apply(lam[, , indices], 1:2, calc_sum),
      sig2 = apply(sig2[, indices], 1, calc_sum),
      rho = apply(rho[, indices], 1, calc_sum),
      tau2 = apply(tau2[, indices], 1, calc_sum),
      info = apply(info[, indices], 1, calc_sum)
    )
    
    rownames(summary$alpha0) <- rownames(summary$sig2) <-
      rownames(summary$tau2) <- rownames(summary$rho) <-
      rownames(summary$info) <- dimnames(summary$lam)[[1]] <-
      c("2.5%", "50%", "97.5%", "mean", "sd")
  }


  save(mod, summary,
    file = paste(path, "/output/", wfile, ".rdata", sep = ""),
    envir = sys.frame(which = sys.nframe())
  )
  return(summary)
}
