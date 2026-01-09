#' Reactive vs. Preemptive Vaccination Simulation Functions
#'
#' This script contains functions to simulate outbreaks across subpopulations and evaluate
#' the relative impact of reactive and preemptive vaccination strategies. Outbreaks follow
#' tent-shaped incidence curves. Vaccination modifies the incidence based on coverage, efficacy,
#' and timing of response.
#'
#' Author: Jong-Hoon Kim
#' Date: 2025-06-27

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
    r_t <- rp - (rp - r0) * (t - tp) / (t1 - tp)  # incidence at time t
    A2 <- 0.5 * (t - tp) * (rp + r_t - 2 * r0)    # subtract r0 for area above baseline

    return(A1 + A2)
  } else {
    # Total area under entire outbreak curve (from t0 to t1)
    A1 <- 0.5 * (tp - t0) * (rp - r0)
    A2 <- 0.5 * (t1 - tp) * (rp - r0)
    return(A1 + A2)
  }
}



#' #' Find the time at which cumulative incidence equals a given target number of cases
#' #'
#' #' @param C Target cumulative number of cases
#' #' @param N Population size
#' #' @param t0, tp, t1 Timepoints defining tent-shaped outbreak
#' #' @param r0, rp Baseline and peak incidence rate (per capita)
#' #'
#' #' @return Time t such that cumulative incidence equals C / N
#' get_time_for_CI <- function(C, N, t0, tp, t1, r0, rp) {
#'   area_target <- C / N # normalize with population size
#'
#'   obj_fun <- function(t) {
#'     get_area(t, t0, tp, t1, r0, rp) - area_target
#'   }
#'   uniroot(obj_fun, lower = t0, upper = t1)$root
#' }
#'
#' #' Generate a list of synthetic regions with population and outbreak probabilities
#' #'
#' #' @param n_region Number of regions to generate
#' #' @param pop_range Range of population sizes (min, max)
#' #' @param prop_U5_range Range of proportion under five in the population
#' #' @param outbreak_prob_range Range of outbreak probabilities
#' #' @return A data frame with one row per region and columns: id, pop_size, outbreak_prob
#' #' population and outbreak probability will be controlled in a more
#' #' sophisticated manner later
#' setup_region <- function(n_region = 1000,
#'                          pop_range = c(1000, 100000),
#'                          prop_U5_range = c(0.15, 0.18),
#'                          outbreak_prob_range = c(0.1, 0.9)) {
#'   tibble::tibble(
#'     id = 1:n_region,
#'     pop_size = sample(seq(pop_range[1], pop_range[2]), n_region, replace = TRUE),
#'     prop_U5 = runif(n_region, prop_U5_range[1], prop_U5_range[2]),
#'     outbreak_prob = runif(n_region, outbreak_prob_range[1], outbreak_prob_range[2])
#'   )
#' }
#'
#'
#' #' Generate outbreaks for each region
#' #'
#' #' This function simulates key outbreak parameters for each region:
#' #' start time \code{t0}, peak time \code{tp}, end time \code{t1}, baseline and peak
#' #' incidence rates \code{r0} and \code{rp}. It ensures that the total per-capita
#' #' incidence (area under the curve) remains less than 1.
#' #'
#' #' \deqn{ \text{Area} = \frac{1}{2}(t_1 - t_0)(r_p - r_0) \le 1 }
#' #'
#' #' Solving for the peak incidence rate:
#' #'
#' #' \deqn{ r_p \le r_0 + \frac{2}{t_1 - t_0} }
#' #'
#' #' @param regions Data frame with `id`, `pop_size`, `prop_U5`, `outbreak_prob`
#' #' @param peak_time_range Range to sample peak times (in days from t0)
#' #' @param duration_range Range to sample time from peak to end
#' #' @param peak_incidence_range Range of candidate values for peak incidence rate
#' #' @param max_attack_rate maximum attack rate allowed
#' #' @param r0 Baseline incidence rate (typically 0)
#' #'
#' #' @return Data frame with outbreak parameters and total cases
#' #' @export
#' generate_outbreaks <- function(regions,
#'                                start_time_range = c(0, 180),
#'                                peak_time_range = c(3*7, 20*7),
#'                                duration_range = c(5*7, 100*7),
#'                                peak_incidence_range = c(1e-6, 0.001), # peak incidence rate per capita
#'                                max_attack_rate = 0.1,
#'                                r0 = 0) {
#' df  <- purrr::map_dfr(1:nrow(regions), function(i) {
#'     region <- regions[i, ]
#'     if (runif(1) < region$outbreak_prob) {
#'       # t0 <- 0
#'       # tp <- runif(1, peak_time_range[1], peak_time_range[2])
#'       # t1 <- tp + runif(1, duration_range[1], duration_range[2])
#'
#'       t0 <- runif(1, start_time_range[1], start_time_range[2])
#'       tp <- t0 + runif(1, peak_time_range[1], peak_time_range[2])
#'       t1 <- tp + runif(1, duration_range[1], duration_range[2])
#'
#'       rp <- runif(1, peak_incidence_range[1], peak_incidence_range[2])
#'
#'       # Area of triangle: A = (1/2)*(base)*(height)
#'       base <- t1 - t0
#'       rp_max <- r0 + 2 / base * max_attack_rate
#'       rp_sample <- runif(1, peak_incidence_range[1], peak_incidence_range[2])
#'       rp <- min(rp_sample, rp_max)
#'       height <- rp - r0
#'       area <- 0.5 * base * height  # total cases per capita
#'
#'       tibble::tibble(
#'         id = region$id,
#'         t0 = t0,
#'         tp = tp,
#'         t1 = t1,
#'         r0 = r0,
#'         rp = rp,
#'         attack_rate_per_capita = area,
#'         total_cases = area * region$pop_size # ignore dynamic change of population susceptible
#'       )
#'     } else {
#'       return(NULL)
#'     }
#'   })
#'   df <- left_join(regions, df, by="id")
#' }


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
#' @param regions A data frame containing at least the following columns:
#'   \itemize{
#'     \item \code{id}: unique region identifier.
#'     \item \code{pop_size}: population size of the region.
#'     \item \code{outbreak_prob}: probability that an outbreak occurs in the region.
#'   }
#'
#' @param start_time_range Numeric vector of length 2 giving the minimum and
#'   maximum possible outbreak start times (e.g., in days from the start of
#'   the planning horizon).
#'
#' @param total_duration_range Numeric vector of length 2 giving the range from
#'   which total outbreak durations are sampled.
#'
#' @param peak_frac_shape Numeric vector of length 2 giving the shape parameters
#'   \eqn{(a, b)} of the Beta distribution used to sample the peak fraction
#'   \eqn{\phi}. Values near \eqn{(2,2)} produce approximately symmetric outbreaks;
#'   skewed values correspond to early- or late-peaking outbreaks.
#'
#' @param peak_incidence_range Numeric vector of length 2 giving the candidate
#'   range for the peak incidence rate (per capita per unit time) prior to
#'   truncation by the attack-rate constraint.
#'
#' @param max_attack_rate Maximum allowable per-capita attack rate for any
#'   outbreak. Peak incidence rates are truncated to ensure this constraint
#'   is satisfied.
#'
#' @param r0 Baseline incidence rate outside the outbreak period. Defaults to 0.
#'
#' @return A data frame with one row per region, including both the original
#'   region-level information and simulated outbreak parameters. Regions in
#'   which no outbreak occurs have \code{NA} values for outbreak-specific fields.
#'
#' @details
#' This function is designed for use in simulation studies comparing
#' pre-emptive and reactive vaccination strategies. The explicit representation
#' of outbreak timing allows reactive vaccination effectiveness to depend on
#' response delays through the remaining area under the incidence curve.
#'
#' The outbreak shape is assumed to be independent of the outbreak probability,
#' conditional on an outbreak occurring.
#'
#' @examples
#' \dontrun{
#' regions <- data.frame(
#'   id = 1:10,
#'   pop_size = rep(1e5, 10),
#'   outbreak_prob = runif(10, 0.1, 0.3)
#' )
#'
#' outbreaks <- generate_outbreaks(regions)
#' }
#'
#' @export
generate_outbreaks <- function(regions,
                               start_time_range = c(0, 180),
                               total_duration_range = c(10 * 7, 120 * 7),
                               peak_frac_shape = c(2, 2),
                               peak_incidence_range = c(1e-6, 0.001),
                               max_attack_rate = 0.1,
                               r0 = 0) {

  purrr::map_dfr(seq_len(nrow(regions)), function(i) {

    region <- regions[i, ]

    ## Determine whether an outbreak occurs
    if (runif(1) >= region$outbreak_prob)
      return(NULL)

    ## Sample outbreak timing
    t0 <- runif(1, start_time_range[1], start_time_range[2])

    D <- runif(1,
               total_duration_range[1],
               total_duration_range[2])

    phi <- rbeta(1, peak_frac_shape[1], peak_frac_shape[2])

    tp <- t0 + phi * D
    t1 <- t0 + D

    ## Sample peak incidence rate subject to attack-rate constraint
    rp_max <- r0 + 2 * max_attack_rate / D

    rp <- min(
      runif(1, peak_incidence_range[1], peak_incidence_range[2]),
      rp_max
    )

    ## Total per-capita attack rate (area under incidence curve)
    area <- 0.5 * D * (rp - r0)

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
      total_cases = area * region$pop_size
    )
  }) |>
    dplyr::left_join(regions, by = "id")
}


#' Assign pre-emptive vaccination using imperfect prioritization based on rank correlation
#'
#' @param regions Data frame of regions (from setup_region or generate_outbreaks)
#' @param n_doses Total number of vaccine doses available
#' @param strategy Strategy for prioritizing regions:
#'        "by_risk", "by_expected_cases", "by_population", or "random"
#' @param coverage_eff Effective coverage per vaccinated person (default = 0.9)
#' @param target_rho Desired Spearman rank correlation between true risk and prioritization score
#'        (only used for "by_risk" and "by_expected_cases")
#'
#' The prioritization ranking is perturbed by performing random swaps on the true ranking.
#' The number of swaps is chosen such that the Spearman correlation is approximately `target_rho`.
#'
#' @return A data frame with added columns: `ranking`, `vaccinated`, `coverage`
#' @export
assign_preemptive_vaccination <- function(regions,
                                          n_doses,
                                          strategy = "by_risk",
                                          coverage_eff = 0.9,
                                          target_rho = 1.0) {

  # Helper: apply random swaps to achieve a target rank correlation
  perturb_to_target_correlation <- function(true_risk, target_rho, max_iter = 5000) {
    n <- length(true_risk)
    order_true <- order(-true_risk)
    best_score <- NULL
    best_rho <- -1

    for (n_swaps in 0:max_iter) {
      perturbed <- order_true
      for (i in seq_len(n_swaps)) {
        idx <- sample(1:n, 2)
        temp <- perturbed[idx[1]]
        perturbed[idx[1]] <- perturbed[idx[2]]
        perturbed[idx[2]] <- temp
      }
      score <- rep(NA_real_, n)
      score[perturbed] <- rev(seq_along(perturbed))
      rho <- suppressWarnings(cor(score, true_risk, method = "spearman"))

      if (is.finite(rho) && abs(rho - target_rho) < 0.01) {
        return(score)
      }
      if (is.finite(rho) && abs(rho - target_rho) < abs(best_rho - target_rho)) {
        best_score <- score
        best_rho <- rho
      }
    }
    return(best_score)
  }

  # Compute expected cases if needed
  regions <- regions %>%
    dplyr::mutate(expected_cases = pop_size * outbreak_prob)

  # Choose ranking strategy
  if (strategy == "by_risk") {
    ranking_var <- if (target_rho >= 1) regions$outbreak_prob else
      perturb_to_target_correlation(regions$outbreak_prob, target_rho)
  } else if (strategy == "by_expected_cases") {
    ranking_var <- if (target_rho >= 1) regions$expected_cases else
      perturb_to_target_correlation(regions$expected_cases, target_rho)
  } else if (strategy == "by_population") {
    ranking_var <- regions$pop_size
  } else if (strategy == "random") {
    ranking_var <- runif(nrow(regions))
  } else {
    stop("Invalid strategy")
  }

  # Sort and assign vaccination
  regions <- regions %>%
    dplyr::mutate(ranking = ranking_var) %>%
    dplyr::arrange(desc(ranking)) %>%
    dplyr::mutate(cumulative_doses = cumsum(pop_size),
                  vaccinated = cumulative_doses <= n_doses,
                  coverage = ifelse(vaccinated, coverage_eff, 0))

  return(regions)
}


#' Assign vaccination strategies to regions with timing and coverage info
#'
#' Computes vaccination start time and effective coverage for each region
#' based on strategy type and outbreak characteristics.
#'
#' @param df Data frame with outbreak features: must include columns `region_id`, `pop_size`, `t0`, `tp`, `t1`, `r0`, `rp`
#' @param strategy Character vector ("none", "preemptive", "reactive"); scalar or length n
#' @param coverage Numeric vector of coverage (0–1); scalar or length n
#' @param direct_vacc_efficacy_U5 Numeric vector of efficacy for 0-4 (0–1); scalar or length n
#' @param direct_vacc_efficacy_5P Numeric vector of efficacy for 5+ (0–1); scalar or length n
#' @param indirect_vacc_efficacy Numeric vector of indirect vaccine efficacy (0–1); scalar or length n
#' @param reactive_threshold Vector of total case count to trigger reactive vaccination; NA or 0 if not used
#' @param delay Vector of days between trigger and vaccine start; scalar or length n
#' @param duration Vector of vaccination program duration (for logging); scalar or length n
#'
#' @return Data frame with columns: region_id, strategy, threshold_case, vaccine_start, coverage, duration
#' @export
assign_vaccination_strategy <- function(df,
                                        strategy = "none",
                                        direct_vacc_efficacy_U5 = 0.35,
                                        direct_vacc_efficacy_5P = 0.7,
                                        indirect_vacc_efficacy = 0.6,
                                        coverage = 0.85,
                                        reactive_threshold = NA,
                                        delay_campaign = 6*7,
                                        delay_immunity = 7,
                                        campaign_duration = 10) {
  n <- nrow(df)

  # Recycle inputs as needed
  strategy <- rep_len(strategy, n)
  vaccine_coverage <- rep_len(coverage, n)
  reactive_threshold <- rep_len(reactive_threshold, n)
  delay_camp <- rep_len(delay_campaign, n)
  delay_imm <- rep_len(delay_immunity, n)
  direct_vacc_eff_U5 <- rep_len(direct_vacc_efficacy_U5, n)
  direct_vacc_eff_5P <- rep_len(direct_vacc_efficacy_5P, n)
  indirect_vacc_eff <- rep_len(indirect_vacc_efficacy, n)
  duration_camp <- rep_len(campaign_duration, n)

  # Initialize result columns
  region_id <- df$region_id
  vaccine_start <- rep(Inf, n)
  vaccine_eff_start <- rep(Inf, n)
  threshold_case <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    strat <- strategy[i]
    cov <- vaccine_coverage[i]
    thresh <- reactive_threshold[i]
    dur <- duration_camp[i]
    dcamp <- delay_camp[i]
    dimm <- delay_imm[i]
    t0 <- df$t0[i]; tp <- df$tp[i]; t1 <- df$t1[i]
    r0 <- df$r0[i]; rp <- df$rp[i]
    N <- df$pop_size[i]

    if (strat == "none") {
      vaccine_start[i] <- Inf
      vaccine_eff_start[i] <- Inf
      vaccine_coverage[i] <- NA
      threshold_case[i] <- NA

    } else if (strat == "preemptive") {
      vaccine_start[i] <- -Inf
      vaccine_eff_start[i] <- -Inf
      vaccine_coverage[i] <- cov
      threshold_case[i] <- NA
    } else if (strat == "reactive") {
      vaccine_coverage[i] <- cov
      threshold_case[i] <- thresh
      tryCatch({
        t_trigger <- get_time_for_CI(C = thresh, N = N, t0 = t0, tp = tp,
                                     t1 = t1, r0 = r0, rp = rp)
        vaccine_start[i] <- t_trigger + dcamp
        vaccine_eff_start[i] <- vaccine_start[i] + dur/2 + dimm
      }, error = function(e) {
        vaccine_start[i] <<- Inf
        vaccine_coverage[i] <<- 0
      })
    } else {
      stop(sprintf("Unknown strategy type for region %s: %s", region_id[i], strat))
    }
  }

  return(tibble::tibble(
    region_id = region_id,
    strategy = strategy,
    threshold_case = threshold_case,
    direct_vacc_efficacy_U5 = direct_vacc_eff_U5,
    direct_vacc_efficacy_5P = direct_vacc_eff_5P,
    indirect_vacc_efficacy = indirect_vacc_eff,
    vaccine_start = vaccine_start,
    vaccine_eff_start = vaccine_eff_start,
    vaccine_coverage = vaccine_coverage,
    duration_campaign = duration_camp
  ))
}

#' Simulate vaccine impact for multiple outbreaks with vaccination
#'
#' Computes total cases and deaths with and without vaccination for each outbreak.
#'
#' @param df Data frame with outbreak features and vaccination info:
#' must include `t0`, `tp`, `t1`, `r0`, `rp`, `vaccine_start`, `vaccinated`, `pop_size`
#' @param cfr Case fatality ratio (default = 0.01)
#'
#' @return Data frame with additional columns:
#' `cases_with_vaccine`, `cases_without_vaccine`, `cases_averted`,
#' `deaths_with_vaccine`, `deaths_without_vaccine`, `deaths_averted`
#' @export
simulate_vaccine_impact <- function(df, cfr = 0.01) {
  n <- nrow(df)

  # Preallocate vectors
  cases_with_vaccine <- numeric(n)
  cases_without_vaccine <- numeric(n)

  for (i in seq_len(n)) {
    with(df[i, ], {

      area_no_vx <- get_area(t1, t0, tp, t1, r0, rp)
      total_cases_no_vx <- area_no_vx * pop_size
      # Calculate fraction of cases to be averted via vaccination
      # Based on direct and indirect effects for under 5 and 5+ age groups
      frac_averted_U5 <- 1 -
          (vaccine_coverage * (1 - direct_vacc_efficacy_U5) * (1 - indirect_vacc_efficacy) +
                                  (1 - vaccine_coverage) * (1 - indirect_vacc_efficacy))
      frac_averted_5P <- 1 -
          (vaccine_coverage * (1 - direct_vacc_efficacy_5P) * (1 - indirect_vacc_efficacy)+
                                  (1 - vaccine_coverage)*(1-indirect_vacc_efficacy))
      if (strategy == "none"){
        total_cases_vx <- NA
      }

      else if (strategy == "preemptive"){
        total_cases_vx <- total_cases_no_vx * (prop_U5 * (1-frac_averted_U5) +
                                            (1-prop_U5) * (1-frac_averted_5P))
      }

      else if(strategy == "reactive") {
        if (vaccine_eff_start >= t1) {
          total_cases_vx <- total_cases_no_vx
        } else {
          area_pre_vx <- get_area(vaccine_eff_start, t0, tp, t1, r0, rp)
          area_post_vx <- area_no_vx - area_pre_vx
          area_post_vx_adj <- area_post_vx * (prop_U5 * (1-frac_averted_U5) +
                                            (1-prop_U5) * (1-frac_averted_5P))
          total_cases_vx <- (area_pre_vx + area_post_vx_adj) * pop_size
        }
      }
      cases_with_vaccine[i] <<- total_cases_vx
      cases_without_vaccine[i] <<- total_cases_no_vx
    })
  }

  # Add to data frame
  df$cases_with_vaccine <- cases_with_vaccine
  df$cases_without_vaccine <- cases_without_vaccine
  df$cases_averted <- df$cases_without_vaccine - df$cases_with_vaccine
  df$deaths_with_vaccine <- df$cases_with_vaccine * cfr
  df$deaths_without_vaccine <- df$cases_without_vaccine * cfr
  df$deaths_averted <- df$deaths_without_vaccine - df$deaths_with_vaccine

  return(df)
}


#' Generate Timestamp
#'
#' Creates a formatted timestamp with configurable components
#'
#' @param year Include year (default: TRUE)
#' @param month Include month (default: TRUE)
#' @param day Include day (default: TRUE)
#' @param hour Include hour (default: FALSE)
#' @param minute Include minute (default: FALSE)
#' @param second Include second (default: FALSE)
#' @return Formatted timestamp string
tstamp <- function(year=TRUE, month=TRUE, day=TRUE,
                   hour=FALSE, minute=FALSE, second=FALSE) {
  # Format date part
  date_format <- ""
  if (year) date_format <- paste0(date_format, "%Y")
  if (month) date_format <- paste0(date_format, "%m")
  if (day) date_format <- paste0(date_format, "%d")

  # Format time part
  time_format <- ""
  if (hour) time_format <- paste0(time_format, "%H")
  if (minute) time_format <- paste0(time_format, "%M")
  if (second) time_format <- paste0(time_format, "%S")

  # Create timestamp
  if (date_format == "" && time_format == "") {
    return("You'd better select parameters well.")
  }

  result <- if (date_format != "") format(Sys.time(), date_format) else ""

  if (time_format != "") {
    time_part <- format(Sys.time(), time_format)
    result <- if (result != "") paste0(result, "T", time_part) else time_part
  }

  return(result)
}
