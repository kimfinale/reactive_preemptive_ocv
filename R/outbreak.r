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


outbreak <- function(region_id) {
  ir_floor <- 0
  ir_peak <- runif(1, 6, 8.5)
  t0 <- 0
  d_peak <- runif(1, 3, 12)
  d_Vx <- runif(1, 3, 12)
  d_VE <- runif(1, 1, 4)  # delay to VE
  d_end <- runif(1, 6, 100)

  t_peak <- t0 + d_peak
  t_Vx <- t0 + d_Vx
  t_VE <- t_Vx + d_VE
  t_end <- t0 + d_end

  return(list(
    region_id = region_id,  # Link to region
    t0 = t0, t_peak = t_peak, t_end = t_end,
    t_Vx = t_Vx, t_VE = t_VE,
    ir_floor = ir_floor, ir_peak = ir_peak
  ))
}

