# kinetics.R
# This file defines viral kinetics functions for various pathogens.
# Each function returns a fixed set of viral dynamics parameters.

kinetics_test <- function() {
  VL_min_infectious <- 0
  VL_floor <- 0
  VL_peak <- 10

  t_latent <- 0
  dt_proliferation <- 5
  dt_clearance <- 5
  t_peak <- t_latent + dt_proliferation
  t_clearance <- t_peak + dt_clearance
  t_Sx <- t_peak + runif(1, -5, 0)

  return(list(t_latent, t_peak, t_clearance, t_Sx, VL_floor, VL_peak, VL_min_infectious))
}

kinetics_flu <- function() {
  VL_min_infectious <- 4
  VL_floor <- 2.95
  VL_peak <- runif(1, 6, 8.5)

  t_latent <- runif(1, 0.5, 1.5)
  dt_proliferation <- runif(1, 1, 3)
  dt_clearance <- runif(1, 2, 3)
  t_peak <- t_latent + dt_proliferation
  t_clearance <- t_peak + dt_clearance
  t_Sx <- t_peak + runif(1, -2, 0)

  return(list(t_latent, t_peak, t_clearance, t_Sx, VL_floor, VL_peak, VL_min_infectious))
}

kinetics_rsv <- function() {
  VL_min_infectious <- 2.8
  VL_floor <- 2.8
  VL_peak <- runif(1, 4, 8)

  t_latent <- runif(1, 2, 4)
  dt_proliferation <- runif(1, 2, 4)
  dt_clearance <- runif(1, 3, 6)
  t_peak <- t_latent + dt_proliferation
  t_clearance <- t_peak + dt_clearance
  t_Sx <- t_peak + runif(1, -1, 1)

  return(list(t_latent, t_peak, t_clearance, t_Sx, VL_floor, VL_peak, VL_min_infectious))
}

kinetics_sarscov2_founder_naive <- function() {
  VL_min_infectious <- 5.5
  VL_floor <- 3
  VL_peak <- rlnorm(1, 1.99950529, 0.19915982)
  while (VL_peak < VL_floor) {
    VL_peak <- rlnorm(1, 1.99950529, 0.19915982)
  }

  t_latent <- runif(1, 2.5, 3.5)
  dt_proliferation <- rlnorm(1, 0.87328419, 0.78801765)
  while (dt_proliferation < 0.5 || dt_proliferation > 10) {
    dt_proliferation <- rlnorm(1, 0.87328419, 0.78801765)
  }
  dt_clearance <- rlnorm(1, 1.95305127, 0.61157974)
  while (dt_clearance < 0.5 || dt_clearance > 25) {
    dt_clearance <- rlnorm(1, 1.95305127, 0.61157974)
  }

  t_peak <- t_latent + dt_proliferation
  t_clearance <- t_peak + dt_clearance
  t_Sx <- t_peak + runif(1, 0, 3)

  return(list(t_latent, t_peak, t_clearance, t_Sx, VL_floor, VL_peak, VL_min_infectious))
}

kinetics_sarscov2_omicron_experienced <- function() {
  VL_min_infectious <- 5.5
  VL_floor <- 3
  VL_peak <- rlnorm(1, 1.87567526, 0.18128526)
  while (VL_peak < VL_floor) {
    VL_peak <- rlnorm(1, 1.87567526, 0.18128526)
  }

  t_latent <- runif(1, 2.5, 3.5)
  dt_proliferation <- rlnorm(1, 1.05318392, 0.68842803)
  while (dt_proliferation < 0.5 || dt_proliferation > 10) {
    dt_proliferation <- rlnorm(1, 1.05318392, 0.68842803)
  }
  dt_clearance <- rlnorm(1, 1.70427144, 0.49046478)
  while (dt_clearance < 0.5 || dt_clearance > 25) {
    dt_clearance <- rlnorm(1, 1.70427144, 0.49046478)
  }

  t_peak <- t_latent + dt_proliferation
  t_clearance <- t_peak + dt_clearance
  t_Sx <- t_peak + runif(1, -5, -1)

  return(list(t_latent, t_peak, t_clearance, t_Sx, VL_floor, VL_peak, VL_min_infectious))
}
