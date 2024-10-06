##############
# Run Up CAR #
##############

# Last Updated Aug 14

# setup ------------------------------------------------------------------------
library(MCMCpack)
library(mvtnorm)
library(mcmc)


# get initial values -----------------------------------------------------------

get_inits_car <- function(
    Y = Y,
    n = n,
    model = c("pois", "binom"),
    start = "result",
    beta0 = numeric(Nt), # baseline
    z = array(0, dim = c(Ns, Nt)), # CAR spatial random effect
    sig2 = numeric(Nt), # vector of spatial variances
    tau2 = numeric(Nt), # vector of non-spatial variances,
    theta = NULL, # vector of either log(prevalence) when pois or logit(incidence) when binom
    qt = array(1, dim = c(Ns, Nt)), # step sizes for theta
    t_accept = array(.4, dim = c(Ns, Nt)), # acceptance rate for theta
    Tn = 1, # iteration
    path = "",
    prev = NA, # batch number of previous run
    numthres = 3, # number of common neighbors
    thres = Inf,
    info = NULL) { # informativeness

  if (is.na(prev)) { # do we have any values from previous runs?
    if (is.null(beta0) | is.null(sig2) | is.null(tau2)) {
      cat("You need to specify beta0, sig2, and tau2 \n")
      return()
    }

    if (!require(here)) {
      warning("here package is not available")
    } else {
      source(here("fun/info_functions.R"))
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

    # spatial random effect from CAR
    if (is.null(z)) {
      z <- theta - matrix(rep(beta0, each = Ns), nrow = Ns)
    }

    # informativenss
    G <- diag(sig2)
    R <- diag(tau2)

    if (is.null(info)) {
      info <- info_value(
        G = G, R = R, beta0 = beta0[, it], numthres = numthres, dist = model
      )
    }

    # save it to mod
    mod <- list(
      Y = Y, n = n,
      model = model,
      beta0 = beta0, theta = theta, z = z, lam = lam,
      sig2 = sig2, tau2 = tau2,
      qt = qt, t_accept = t_accept,
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


# run up mcar with inverse wishart prior ---------------------------------------

run_up_car <- function(
    start = "",
    T_inc = 10, # number of iterations in this batch
    wfile = "current",
    path = "",
    thin = 1,
    a_t = 1, b_t = 1 / 100, a_s = 1, b_s = 1 / 7, # hyperparams
    numthres = 3, # number of common neighbors
    thres = Inf,
    get_DIC = FALSE) {
  if (!require(here)) {
    warning("here package is not available")
  } else {
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
  z <- mod$z
  theta <- mod$theta

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
  t_accept <- array(0, dim = c(Ns, Nt))

  time_start <- Sys.time()

  for (it in (oldT + 1):Tn) {
    # update beta0 -------------------------------------------------------------

    beta_var <- tau2[, it - 1] / Ns
    beta_mean <- apply(theta[, ] - z[, ], 2, sum) / Ns

    if (model == "binom" && thres < Inf) { ## if restricted (only applies to binomial model)

      tvar <- tau2[, it - 1] + (tau2[, it - 1] + sig2[, it - 1]) / numthres

      if (mod$total + it > 50) {
        plus <- ((1 - thres) + sqrt((thres - 1)^2 + 4 * (thres - 1 / tvar))) / 2
        plus <- ifelse(plus > 0, ifelse(plus < 1, plus, 1), 0)
        plus <- logit(plus)
        minu <- ((1 - thres) - sqrt((thres - 1)^2 + 4 * (thres - 1 / tvar))) / 2
        minu <- ifelse(minu > 0, ifelse(minu < 1, minu, 1), 0)
        minu <- logit(minu)
        u <- runif(Nt, pnorm(minu, beta.mean, sqrt(beta.var)), pnorm(plus, beta.mean, sqrt(beta.var)))
      } else {
        u <- runif(Nt, 0, 1)
      }

      beta0[, s, it] <- qnorm(u, beta.mean, sqrt(beta.var))
    } else { ## either standard model or Poisson restricted model

      beta0[, it] <- rnorm(Nt, beta_mean, sd = sqrt(beta_var))
    }

    # update z -----------------------------------------------------------------

    for (i in 1:Ns) {
      if (num[i] > 1) {
        muia <- apply(z[neigh[[i]], ], 2, mean)
      } else {
        muia <- z[neigh[[i]], ]
      }

      sigia <- sig2[, it - 1] / num[i]
      z_var <- 1 / (1 / tau2[, it - 1] + 1 / sigia)
      z_mean <- z_var * ((theta[i, ] - beta0[, it]) / tau2[, it - 1] + muia / sigia)
      z[i, ] <- rnorm(Nt, z_mean, sd = sqrt(z_var))
    }

    z[, ] <- z[, ] - rep(apply(z[, ], 2, mean), each = Ns)


    # update theta -------------------------------------------------------------

    ts <- rnorm(Ns * Nt, theta[, ], qt) # candidate
    ra <- Y[, ] * (ts - theta[, ])
    if (model == "pois") {
      rb <- n * (exp(ts) - exp(theta[, ]))
    } else if (model == "binom") {
      rb <- n * (log(1 + exp(ts)) - log(1 + exp(theta[, ])))
    }
    rc <- ((ts - X %*% t(beta0[, it]) - z[, ])^2 -
      (theta[, ] - X %*% t(beta0[, it]) - z[, ])^2)
    r <- exp(ra - rb - rc / (2 * rep(tau2[, it - 1], each = Ns)))
    accept <- (r > runif(Ns * Nt))
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

    ass <- Ns / 2 + a_t
    bss <- apply((theta[, ] - X %*% t(beta0[, it]) - z[, ])^2, 2, sum) / 2 + b_t

    if (thres < Inf && model == "pois") { # restricted Pois

      tau2_thres <- (log(1 / thres + 1) - sig2[, it - 1] / numthres) / (1 + 1 / numthres)
      tau2_thres <- ifelse(tau2_thres > 0, tau2_thres, 0)
      u <- runif(Nt, 0, pgamma(1 / tau2_thres, ass, bss))
      tau2[, it] <- 1 / qgamma(u, ass, bss)
    } else if (thres < Inf && model == "binom") { # restricted Binom

      mu <- exp(beta0[, it]) / (1 + exp(beta0[, it]))
      tau2_thres <- (1 / ((thres + mu) * (1 - mu)) - sig2[, it - 1] / numthres) / (1 + 1 / numthres)
      tau2_thres <- ifelse(tau2_thres > 0, tau_thres, 0)
      u <- runif(Nt, 0, pgamma(1 / tau2_thres, ass, bss))
      tau2[, it] <- 1 / qgamma(u, ass, bss)
    } else if (thres == Inf) { # standard

      tau2[, it] <- 1 / rgamma(Nt, ass, bss)
    }


    # update sig2 (spatial variance) -------------------------------------------

    bss <- rep(0, Nt)

    for (i in 1:Ns) {
      if (Nt == 1) {
        sum_z_j <- sum(z[neigh[[i]], ])
      } else {
        sum_z_j <- apply(z[neigh[[i]], ], 2, sum)
      }
      bss <- bss + z[i, ]^2 * num[i] - z[i, ] * sum_z_j
    }
    ass <- (Ns - 1) / 2 + a_s
    bss <- bss / 2 + b_s

    if (thres < Inf && model == "pois") { # restricted Pois

      sig2_thres <- (log(1 / thres + 1) - tau2[, it] * (1 + 1 / numthres)) * numthres
      sig2_thres <- ifelse(sig2_thres > 0, sig2_thres, 0)
      u <- runif(Nt, 0, pgamma(1 / sig2_thres, ass, bss))
      sig2[, it] <- 1 / qgamma(u, ass, bss)
    } else if (thres < Inf && model == "binom") { # restricted Binom

      mu <- exp(beta0[, it]) / (1 + exp(beta0[, it]))
      sig2_thres <- (1 / ((thres + mu) * (1 - mu)) - (1 + 1 / numthres) * tau2[, it]) * numthres
      sig2_thres <- ifelse(sig2_thres > 0, sig2_thres, 0)
      u <- runif(Nt, 0, pgamma(1 / sig2_thres, ass, bss))
      sig2[, it] <- 1 / qgamma(u, ass, bss)
    } else if (thres == Inf) { # standard

      sig2[, it] <- 1 / rgamma(Nt, ass, bss)
    }



    # compute the informativeness
    info[, it] <- info_value(
      G = diag(sig2[, it]), R = diag(tau2[, it]), beta0 = beta0[, it],
      dist = model, numthres = numthres
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
    beta0 = beta0[, Tn], theta = theta, lam = lam[, , Tn],
    z = z, sig2 = sig2[, Tn], tau2 = tau2[, Tn],
    info = info[, Tn], numthres = numthres, thres = thres,
    seed = .Random.seed, t_accept = t_accept, qt = qt,
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
    beta0 = beta0[, indices], theta = theta[, ], lam = lam[, , indices],
    z = z, sig2 = sig2[, indices], tau2 = tau2[, indices],
    info = info[, indices], numthres = numthres, thres = thres,
    seed = .Random.seed, t_accept = t_accept, qt = qt,
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

  summary <- list(
    beta0 = apply(beta0[, indices], 1, calc_sum),
    lam = apply(lam[, , indices], 1:2, calc_sum),
    sig2 = apply(sig2[, indices], 1, calc_sum),
    tau2 = apply(tau2[, indices], 1, calc_sum),
    info = apply(info[, indices], 1, calc_sum)
  )

  rownames(summary$beta0) <- rownames(summary$sig2) <-
    rownames(summary$tau2) <- rownames(summary$info) <-
    dimnames(summary$lam)[[1]] <- c("2.5%", "50%", "97.5%", "mean", "sd")

  save(mod, summary,
    file = paste(path, "/output/", wfile, ".rdata", sep = ""),
    envir = sys.frame(which = sys.nframe())
  )
  return(summary)
}
