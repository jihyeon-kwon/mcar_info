#########################################
# function to combine simulation output #
#########################################




# setup -----

library(here)



# function -----
# path <- here("simulation3/output/realpop_har_rho_0.1_dataset_")
# suff <- ""
# type <- "rdata"

combine_summary <- function(
    path = path, suff = c(),
    type = c("mcar", "car"),
    rho_loc = NULL, # specify what rho you want if it is not AR(1)
    save_to = NULL) {
  num_datasets <- length(suff)

  # create empty arrays to store results
  load(paste0(path, suff[1], ".rdata"))

  Ns <- dim(mod$lam)[[1]]
  Nt <- dim(mod$beta0)[[1]]

  if (type == "mcar") {
    rho <- array(NA, dim = c(5, num_datasets))
    dimnames(rho)[[1]] <- c("2.5%", "50%", "97.5%", "mean", "sd")
    dimnames(rho)[[2]] <- 1:num_datasets
  }

  sig2 <- array(NA, dim = c(5, Nt, num_datasets))
  tau2 <- array(NA, dim = c(5, Nt, num_datasets))
  beta0 <- array(NA, dim = c(5, Nt, num_datasets))
  info <- array(NA, dim = c(5, Nt, num_datasets))
  lam <- array(NA, dim = c(5, Ns, Nt, num_datasets))
  DIC <- pD <- Dbar <- rep(NA, num_datasets)

  dimnames(lam)[[2]] <- 1:Ns

  dimnames(sig2)[[1]] <- dimnames(tau2)[[1]] <- dimnames(beta0)[[1]] <-
    dimnames(info)[[1]] <- dimnames(lam)[[1]] <- c("2.5%", "50%", "97.5%", "mean", "sd")

  dimnames(sig2)[[2]] <- dimnames(tau2)[[2]] <- dimnames(beta0)[[2]] <-
    dimnames(info)[[2]] <- dimnames(lam)[[3]] <- 1:Nt

  dimnames(sig2)[[3]] <- dimnames(tau2)[[3]] <- dimnames(beta0)[[3]] <-
    dimnames(info)[[3]] <- dimnames(lam)[[4]] <- suff

  if (is.null(rho_loc)) {
    for (i in 1:num_datasets) {
      load(paste0(path, suff[i], ".rdata"))

      if (type == "mcar") {
        rho[, i] <- summary$rho[]
      }
      lam[, , , i] <- summary$lam[, , ]
      sig2[, , i] <- summary$sig2[, ]
      tau2[, , i] <- summary$tau2[, ]
      beta0[, , i] <- summary$beta0[, ]
      info[, , i] <- summary$info[, ]
      DIC[i] <- mod$DIC
      pD[i] <- mod$pD
      Dbar[i] <- mod$Dbar

      cat(paste0("Finished: ", suff[i]), "\n")
    }
  } else {
    for (i in 1:num_datasets) {
      load(paste0(path, suff[i], ".rdata"))

      if (type == "mcar") {
        rho[, i] <- summary$rho[, rho_loc]
      }
      lam[, , , i] <- summary$lam[, , ]
      sig2[, , i] <- summary$sig2[, ]
      tau2[, , i] <- summary$tau2[, ]
      beta0[, , i] <- summary$beta0[, ]
      info[, , i] <- summary$info[, ]
      DIC[i] <- mod$DIC
      pD[i] <- mod$pD
      Dbar[i] <- mod$Dbar

      cat(paste0("Finished: ", suff[i]), "\n")
    }
  }

  combined_summary <- list(
    lam = lam, sig2 = sig2, tau2 = tau2, beta0 = beta0, info = info,
    DIC = DIC, pD = pD, Dbar = Dbar
  )

  if (type == "mcar") {
    combined_summary$rho <- rho
  }
  if (is.null(save_to)) {
    return(combined_summary)
  } else {
    save(combined_summary, file = save_to)
  }
}


# combine_summary(
#   path = here("application/output/pa_teen_firearm_car_"),
#   suff = c("M-Black", "M-White", "F-Black", "F-White"),
#   type = "car"
# )

combine_samples <- function(
    path = path, suff = c(),
    type = c("mcar", "car"),
    rho_loc = NULL, # specify what rho you want if it is not AR(1)
    save_to = NULL) {
  num_datasets <- length(suff)

  # create empty arrays to store results
  load(paste0(path, suff[1], ".rdata"))

  Ns <- dim(mod$lam)[[1]]
  Nt <- dim(mod$beta0)[[1]]
  n_samples <- dim(mod$beta0)[[2]]

  if (type == "mcar") {
    rho <- array(NA, dim = c(n_samples, num_datasets))
    dimnames(rho)[[1]] <- 1:n_samples
    dimnames(rho)[[2]] <- 1:num_datasets
  }

  sig2 <- array(NA, dim = c(Nt, n_samples, num_datasets))
  tau2 <- array(NA, dim = c(Nt, n_samples, num_datasets))
  beta0 <- array(NA, dim = c(Nt, n_samples, num_datasets))
  info <- array(NA, dim = c(Nt, n_samples, num_datasets))
  lam <- array(NA, dim = c(Ns, Nt, n_samples, num_datasets))

  dimnames(lam)[[1]] <- 1:Ns

  dimnames(sig2)[[1]] <- dimnames(tau2)[[1]] <- dimnames(beta0)[[1]] <-
    dimnames(info)[[1]] <- dimnames(lam)[[2]] <- 1:Nt

  dimnames(sig2)[[2]] <- dimnames(tau2)[[2]] <- dimnames(beta0)[[2]] <-
    dimnames(info)[[2]] <- dimnames(lam)[[3]] <- 1:n_samples

  dimnames(sig2)[[3]] <- dimnames(tau2)[[3]] <- dimnames(beta0)[[3]] <-
    dimnames(info)[[3]] <- dimnames(lam)[[4]] <- suff

  if (is.null(rho_loc)) {
    for (i in 1:num_datasets) {
      load(paste0(path, suff[i], ".rdata"))

      if (type == "mcar") {
        rho[, i] <- mod$rho[]
      }
      lam[, , , i] <- mod$lam[, , ]
      sig2[, , i] <- mod$sig2[, ]
      tau2[, , i] <- mod$tau2[, ]
      beta0[, , i] <- mod$beta0[, ]
      info[, , i] <- mod$info[, ]

      cat(paste0("Finished: ", suff[i]), "\n")
    }
  } else {
    for (i in 1:num_datasets) {
      load(paste0(path, suff[i], ".rdata"))

      if (type == "mcar") {
        rho[, i] <- mod$rho[, rho_loc]
      }
      lam[, , , i] <- mod$lam[, , ]
      sig2[, , i] <- mod$sig2[, ]
      tau2[, , i] <- mod$tau2[, ]
      beta0[, , i] <- mod$beta0[, ]
      info[, , i] <- mod$info[, ]

      cat(paste0("Finished: ", suff[i]), "\n")
    }
  }

  combined_samples <- list(
    lam = lam, sig2 = sig2, tau2 = tau2, beta0 = beta0, info = info
  )

  if (type == "mcar") {
    combined_samples$rho <- rho
  }
  if (is.null(save_to)) {
    return(combined_samples)
  } else {
    save(combined_samples, file = save_to)
  }
}



combine_info <- function(
    path = path, suff = c(),
    save_to = NULL) {
  num_datasets <- length(suff)

  # create empty arrays to store results
  load(paste0(path, suff[1], ".rdata"))

  Ns <- dim(mod$lam)[[1]]
  Nt <- dim(mod$beta0)[[1]]
  n_samples <- dim(mod$beta0)[[2]]
  info <- array(NA, dim = c(Nt, n_samples, num_datasets))


  dimnames(info)[[1]] <- 1:Nt
  dimnames(info)[[2]] <- 1:n_samples
  dimnames(info)[[3]] <- suff

  for (i in 1:num_datasets) {
    load(paste0(path, suff[i], ".rdata"))
    info[, , i] <- mod$info[, ]
    cat(paste0("Finished: ", suff[i]), "\n")
  }

  combined_info <- list(
    info = info
  )
  if (is.null(save_to)) {
    return(combined_info)
  } else {
    save(combined_info, file = save_to)
  }
}
