#' Compute instantaneous infection rate for tent-shaped outbreak
#'
#' Given a tent-shaped outbreak curve defined by start ($t_0$), peak ($t_p$), and end ($t_1$) times,
#' this function returns the infection rate at a specific time point $t$.
#'
#' @param t Time at which to compute the infection rate
#' @param t0 Time of outbreak start
#' @param tp Time of outbreak peak
#' @param t1 Time of outbreak end
#' @param r0 Baseline incidence rate (typically 0)
#' @param rp Peak incidence rate per capita
#'
#' @return Numeric value of the infection rate at time \code{t}
#' @examples
#' get_IR(t = 5, t0 = 0, tp = 10, t1 = 30, r0 = 0, rp = 0.02)
get_IR <- function(t, t0, tp, t1, r0, rp) {
  if (t <= tp) {
    r0 + (rp - r0) * (t - t0) / (tp - t0)
  } else {
    rp - (rp - r0) * (t - tp) / (t1 - tp)
  }
}

#' Compute area up to time t
#'
#' Integrates the incidence rate curve from \code{t0} to \code{t} assuming a triangular
#' (tent-shaped) outbreak with peak at \code{tp}. Uses piecewise linear incidence.
#'
#' @param t Time up to which area is calculated
#' @param t0 Time of outbreak start
#' @param tp Time of outbreak peak
#' @param t1 Time of outbreak end
#' @param r0 Baseline incidence rate
#' @param rp Peak incidence rate
#'
#' @return Area under incidence rate curve up to time t
#' @examples
#' get_area(t = 15, t0 = 0, tp = 10, t1 = 30, r0 = 0, rp = 0.02)
get_area <- function(t, t0, tp, t1, r0, rp) {
  if (t <= tp) {
    # Rising slope: triangle from t0 to t
    h <- r0 + (rp - r0) * (t - t0) / (tp - t0) - r0
    return(0.5 * (t - t0) * h)
  } else if (t <= t1) {
    # Rising triangle: t0 to tp
    A1 <- 0.5 * (tp - t0) * (rp - r0)

    # Falling trapezoid: tp to t
    r_t <- rp - (rp - r0) * (t - tp) / (t1 - tp) # incidence at time t
    A2 <- 0.5 * (t - tp) * (rp + r_t - 2 * r0) # subtract r0 for area above baseline

    return(A1 + A2)
  } else {
    # Total area under entire outbreak curve (from t0 to t1)
    A1 <- 0.5 * (tp - t0) * (rp - r0)
    A2 <- 0.5 * (t1 - tp) * (rp - r0)
    return(A1 + A2)
  }
}



#' Find the time at which cumulative incidence equals a given target number of cases
#'
#' @param C Target cumulative number of cases
#' @param N Population size
#' @param t0, tp, t1 Timepoints defining tent-shaped outbreak
#' @param r0, rp Baseline and peak incidence rate (per capita)
#'
#' @return Time t such that cumulative incidence equals C / N
get_time_for_CI <- function(C, N, t0, tp, t1, r0, rp) {
  area_target <- C / N # normalize with population size

  obj_fun <- function(t) {
    get_area(t, t0, tp, t1, r0, rp) - area_target
  }
  uniroot(obj_fun, lower = t0, upper = t1)$root
}

#' Generate a list of synthetic regions with population and outbreak probabilities
#'
#' @param n_region Number of regions to generate
#' @param pop_range Range of population sizes (min, max)
#' @param prop_U5_range Range of proportion under five in the population
#' @param outbreak_prob_range Range of outbreak probabilities
#' @return A data frame with one row per region and columns: id, pop_size, outbreak_prob
#' population and outbreak probability will be controlled in a more
#' sophisticated manner later
setup_region <- function(n_region = 1000,
                         pop_range = c(1000, 100000),
                         prop_U5_range = c(0.15, 0.18),
                         outbreak_prob_range = c(0.1, 0.9)) {
  tibble::tibble(
    id = 1:n_region,
    pop_size = sample(seq(pop_range[1], pop_range[2]), n_region, replace = TRUE),
    prop_U5 = runif(n_region, prop_U5_range[1], prop_U5_range[2]),
    outbreak_prob = runif(n_region, outbreak_prob_range[1], outbreak_prob_range[2])
  )
}


#' Generate outbreaks for each region
#'
#' This function simulates key outbreak parameters for each region:
#' start time \code{t0}, peak time \code{tp}, end time \code{t1}, baseline and peak
#' incidence rates \code{r0} and \code{rp}. It ensures that the total per-capita
#' incidence (area under the curve) remains less than 1.
#'
#' \deqn{ \text{Area} = \frac{1}{2}(t_1 - t_0)(r_p - r_0) \le 1 }
#'
#' Solving for the peak incidence rate:
#'
#' \deqn{ r_p \le r_0 + \frac{2}{t_1 - t_0} }
#'
#' @param regions Data frame with `id`, `pop_size`, `prop_U5`, `outbreak_prob`
#' @param peak_time_range Range to sample peak times (in days from t0)
#' @param duration_range Range to sample time from peak to end
#' @param peak_incidence_range Range of candidate values for peak incidence rate
#' @param max_attack_rate maximum attack rate allowed
#' @param r0 Baseline incidence rate (typically 0)
#'
#' @return Data frame with outbreak parameters and total cases
#' @export
generate_outbreaks <- function(regions,
                               start_time_range = c(0, 180),
                               peak_time_range = c(3 * 7, 20 * 7),
                               duration_range = c(5 * 7, 100 * 7),
                               peak_incidence_range = c(1e-6, 0.001), # peak incidence rate per capita
                               max_attack_rate = 0.1,
                               r0 = 0) {
  df <- purrr::map_dfr(1:nrow(regions), function(i) {
    region <- regions[i, ]
    if (runif(1) < region$outbreak_prob) {
      # t0 <- 0
      # tp <- runif(1, peak_time_range[1], peak_time_range[2])
      # t1 <- tp + runif(1, duration_range[1], duration_range[2])

      t0 <- runif(1, start_time_range[1], start_time_range[2])
      tp <- t0 + runif(1, peak_time_range[1], peak_time_range[2])
      t1 <- tp + runif(1, duration_range[1], duration_range[2])

      rp <- runif(1, peak_incidence_range[1], peak_incidence_range[2])

      # Area of triangle: A = (1/2)*(base)*(height)
      base <- t1 - t0
      rp_max <- r0 + 2 / base * max_attack_rate
      rp_sample <- runif(1, peak_incidence_range[1], peak_incidence_range[2])
      rp <- min(rp_sample, rp_max)
      height <- rp - r0
      area <- 0.5 * base * height # total cases per capita

      tibble::tibble(
        id = region$id,
        t0 = t0,
        tp = tp,
        t1 = t1,
        r0 = r0,
        rp = rp,
        attack_rate_per_capita = area,
        total_cases = area * region$pop_size # ignore dynamic change of population susceptible
      )
    } else {
      return(NULL)
    }
  })
  df <- left_join(regions, df, by = "id")
}


#' Generate stylized outbreak trajectories for each region
#'
#' This function simulates simplified outbreak dynamics at the region (grid-cell)
#' level. Outbreak occurrence is first determined by a Bernoulli trial with
#' region-specific outbreak probability. Conditional on an outbreak occurring,
#' the outbreak trajectory is represented by a triangular incidence curve
#' defined by a start time, peak time, and end time.
#'
#' Timing is parameterized using:
#' \itemize{
#'   \item a start time \eqn{t_0},
#'   \item a total outbreak duration \eqn{D = t_1 - t_0},
#'   \item a peak fraction \eqn{\phi \in (0,1)} such that
#'         \eqn{t_p = t_0 + \phi D}.
#' }
#'
#' This parameterization enforces the constraint \eqn{t_0 < t_p < t_1} by
#' construction and avoids rejection sampling.
#'
#' Incidence is assumed to increase linearly from a baseline rate \eqn{r_0}
#' at \eqn{t_0} to a peak incidence rate \eqn{r_p} at \eqn{t_p}, and then decrease
#' linearly back to \eqn{r_0} at \eqn{t_1}. The total per-capita attack rate
#' corresponds to the area under this triangular curve:
#'
#' \deqn{
#' A = \frac{1}{2} D (r_p - r_0).
#' }
#'
#' The peak incidence rate is truncated to ensure that the total attack rate
#' does not exceed \code{max_attack_rate}.
#'
#' @param regions Data frame containing at least:
#'   \code{id}, \code{pop_size}, and \code{outbreak_prob}.
#' @param start_time_range Range for outbreak start times.
#' @param total_duration_range Range for total outbreak duration.
#' @param peak_frac_shape Shape parameters (a, b) for Beta distribution of
#'   the peak fraction \eqn{\phi}.
#' @param peak_incidence_range Candidate range for peak incidence rate.
#' @param max_attack_rate Maximum per-capita attack rate.
#' @param r0 Baseline incidence rate (default = 0).
#'
#' @return Data frame with one row per region, including outbreak timing
#'   and severity parameters. Regions without outbreaks have NA values
#'   for outbreak-specific fields.
#'
#' @details
#' Outbreak shape is assumed independent of outbreak probability,
#' conditional on outbreak occurrence.
#'
#' @export
assign_outbreak_trajectory <- function(regions,
                                       start_time_range = c(0, 180),
                                       total_duration_range = c(10 * 7, 120 * 7),
                                       peak_frac_shape = c(2, 2),
                                       peak_incidence_range = c(1e-5, 5e-4),
                                       cfr_range = c(0.001, 0.1),
                                       max_attack_rate = 0.1,
                                       r0 = 0) {
  purrr::map_dfr(seq_len(nrow(regions)), function(i) {
    region <- regions[i, ]

    ## --- outbreak occurrence ---
    # if (runif(1) >= region$outbreak_prob)
    if (!region$Y) { # outbreak occurred
      return(
        tibble::tibble(
          id = region$id,
          t0 = NA_real_,
          tp = NA_real_,
          t1 = NA_real_,
          duration = NA_real_,
          peak_fraction = NA_real_,
          r0 = NA_real_,
          rp = NA_real_,
          attack_rate_per_capita = NA_real_,
          case_fatality_ratio = NA_real_
        )
      )
    }

    ## --- timing ---
    t0 <- runif(1, start_time_range[1], start_time_range[2])
    D <- runif(1, total_duration_range[1], total_duration_range[2])
    phi <- rbeta(1, peak_frac_shape[1], peak_frac_shape[2])

    tp <- t0 + phi * D
    t1 <- t0 + D

    ## --- peak incidence with attack-rate constraint ---
    rp_max <- r0 + 2 * max_attack_rate / D

    rp <- min(
      runif(1, peak_incidence_range[1], peak_incidence_range[2]),
      rp_max
    )

    ## --- total attack rate ---
    area <- 0.5 * D * (rp - r0)
    cfr <- runif(1, cfr_range[1], cfr_range[2])

    tibble::tibble(
      id = region$id,
      t0 = t0,
      tp = tp,
      t1 = t1,
      duration = D,
      peak_fraction = phi,
      r0 = r0,
      rp = rp,
      attack_rate_per_capita = area,
      case_fatality_ratio = cfr
    )
  }) |>
    dplyr::left_join(regions, by = "id")
}

#' Generate Exponential Outbreak Curve
#'
#' @param duration Total duration of the outbreak (days).
#' @param peak_time Time of the peak (days from start).
#' @param peak_incidence Peak incidence (cases per day).
#' @param t_vec Vector of time points to evaluate.
#'
#' @return Vector of incidence values.
#' @export
generate_exponential_curve <- function(duration, peak_time, peak_incidence, t_vec) {
  # Assume start at t=0 implies I(0) is small but non-zero for exponential fit logic,
  # but here we constrain the "width" by growth/decay rates.
  # For simplicity, we define start and end as where incidence drops to ~1% of peak or similar,
  # OR we just fit r so that it approximately matches the duration.

  # Let's solve r given we want the curve to be substantial within [0, duration].
  # We assume I(0) = I(duration) = epsilon.
  epsilon <- peak_incidence / 100

  # Growth phase: epsilon = I_peak * exp(r_g * (0 - t_peak))
  # ln(eps/I_peak) = -r_g * t_peak  =>  r_g = -ln(eps/I_peak) / t_peak
  if (peak_time > 0) {
    r_growth <- -log(epsilon / peak_incidence) / peak_time
  } else {
    r_growth <- 1 # steep
  }

  # Decay phase: epsilon = I_peak * exp(-r_d * (duration - t_peak))
  # ln(eps/I_peak) = -r_d * (duration - t_peak) => r_d = -ln(eps/I_peak) / (duration - t_peak)
  if (duration > peak_time) {
    r_decay <- -log(epsilon / peak_incidence) / (duration - peak_time)
  } else {
    r_decay <- 1 # steep
  }

  incidence <- numeric(length(t_vec))

  # Calculate for each t
  for (i in seq_along(t_vec)) {
    t <- t_vec[i]
    if (t < 0 || t > duration) {
      incidence[i] <- 0
    } else if (t <= peak_time) {
      incidence[i] <- peak_incidence * exp(r_growth * (t - peak_time))
    } else {
      incidence[i] <- peak_incidence * exp(-r_decay * (t - peak_time))
    }
  }

  return(incidence)
}
