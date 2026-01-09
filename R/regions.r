# kinetics.R
# This file defines outbreak dynamics functions for cholera.
# Each function returns a fixed set of outbreak dynamics parameters.

outbreak_test <- function() {
  ir_min_infectious <- 0
  ir_floor <- 0
  ir_peak <- 10

  t_peak <- 5
  t_ <- t_latent + dt_proliferation
  t_clearance <- t_peak + dt_clearance
  t_Sx <- t_peak + runif(1, -5, 0)

  return(list(t_latent, t_peak, t_clearance, t_Sx, VL_floor, VL_peak, VL_min_infectious))
}

outbreak_ssa <- function() {

  ir_floor <- 0 # per 100,000 person
  ir_peak <- runif(1, 2, 8.5)
  d_peak <- runif(1, 2, 30) # delay to peak incidence
  d_end <- runif(1, 6, 100) # delay to the end of the outbreak
  t0 <- 0 # outbreak begin
  d_Vx <- runif(1, 2, 10) # delay to vaccination
  d_VE <- 3 # delay to vaccination effectiveness

  t_peak <- t0 + d_peak
  t_Vx <- t0 + d_Vx
  t_VE <- t0 + d_Vx + d_VE
  t_end <- t0 + d_end

  return(list(t0, t_peak, t_end, t_Vx, t_VE, ir_floor, ir_peak))
}
