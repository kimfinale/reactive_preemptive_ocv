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
                               peak_time_range = c(3*7, 20*7),
                               duration_range = c(5*7, 100*7),
                               peak_incidence_range = c(1e-6, 0.001), # peak incidence rate per capita
                               max_attack_rate = 0.1,
                               r0 = 0) {
df  <- purrr::map_dfr(1:nrow(regions), function(i) {
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
      area <- 0.5 * base * height  # total cases per capita

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
  df <- left_join(regions, df, by="id")
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


#' Compute optimal number of regions (k*) to preemptively vaccinate
#'
#' @param p Vector of outbreak probabilities for each region
#' @param C Vector of outbreak burdens for each region (if unvaccinated)
#' @param VE Vaccine efficacy (0 to 1)
#' @param Cov Vaccine coverage (0 to 1)
#' @param delta Reduction factor for reactive vaccination impact (0 to 1)
#'
#' @return A list with optimal k, benefit by k, and full benefit curve
optimal_k_preemptive <- function(p, C, VE, Cov, delta) {
  stopifnot(length(p) == length(C))
  M <- length(p)
  score <- p * C
  rank_idx <- order(score, decreasing = TRUE)
  ordered_score <- score[rank_idx]
  S_total <- sum(score)
  S_k <- cumsum(ordered_score)
  benefit_k <- VE * Cov * ((1 - delta) * S_k + delta * S_total)
  k_star <- which.max(benefit_k)

  return(list(
    k_star = k_star,
    max_benefit = benefit_k[k_star],
    benefit_curve = benefit_k,
    ranks = rank_idx
  ))
}


#' Generate a risk profile for regions or individuals
#'
#' This function samples a numeric risk score from a specified distribution,
#' commonly used to represent outbreak risk or infection probability across subpopulations.
#' Supported distributions include Beta, Uniform, and Truncated Normal (bounded in [0, 1]).
#'
#' The output vector is sorted in decreasing order to facilitate risk ranking.
#'
#' @param n Integer. Number of samples to generate.
#' @param dist Character. Distribution to sample from: `"beta"`, `"uniform"`, or `"normal"`.
#' @param params List. Distribution parameters:
#'   - For `"beta"`: `a`, `b` (shape parameters)
#'   - For `"uniform"`: `min`, `max` (bounds)
#'   - For `"normal"`: `mean`, `sd` (mean and standard deviation; truncated to [0,1])
#' @param seed Integer or NULL. Optional random seed for reproducibility.
#'
#' @return Numeric vector of length `n`, sorted in decreasing order.
#' @examples
#' generate_risk_profile(10, dist = "beta", params = list(a = 2, b = 5), seed = 42)
#' generate_risk_profile(10, dist = "uniform", params = list(min = 0, max = 1))
#' generate_risk_profile(10, dist = "normal", params = list(mean = 0.5, sd = 0.1))
#' @export
generate_risk_profile <- function(n, dist = "beta", params = list(), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Ensure truncnorm is available if needed
  if (dist == "normal" && !requireNamespace("truncnorm", quietly = TRUE)) {
    stop("Please install the 'truncnorm' package: install.packages('truncnorm')")
  }

  # Helper for default parameter values
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # Sample from the specified distribution
  p_raw <- switch(
    dist,
    "beta" = {
      a <- params$a %||% 2
      b <- params$b %||% 5
      rbeta(n, a, b)
    },
    "uniform" = {
      min_val <- params$min %||% 0
      max_val <- params$max %||% 1
      runif(n, min_val, max_val)
    },
    "normal" = {
      mu <- params$mean %||% 0.5
      sigma <- params$sd %||% 0.1
      truncnorm::rtruncnorm(n, a = 0, b = 1, mean = mu, sd = sigma)
    },
    stop("Unsupported distribution. Choose from 'beta', 'uniform', or 'normal'.")
  )

  # Return sorted risk profile (e.g., for prioritization)
  sort(p_raw, decreasing = TRUE)
}

#' Compute total expected cost under a mixed preemptive and reactive strategy
#'
#' This function calculates the total expected cost of a hybrid strategy where
#' a fraction `alpha` of the population is preemptively vaccinated (subject to
#' vaccine supply constraint `f`), and the remainder is managed reactively after
#' outbreaks occur. The risk profile should be sorted in decreasing order.
#'
#' Preemptively vaccinated regions incur a fixed cost `C` per region. Reactively
#' managed regions experience an expected cost proportional to their risk, weighted
#' by $(1 - r) \cdot D$, where `r` is the relative effectiveness of reactive
#' response compared to preemptive vaccination.
#'
#' @param alpha Numeric. Proportion of the population targeted for preemptive vaccination (between 0 and 1).
#' @param risk_profile Numeric vector. Outbreak risk scores for each region, sorted in decreasing order.
#' @param C Numeric. Cost of preemptive vaccination per region.
#' @param D Numeric. Cost of an unmitigated outbreak per region.
#' @param r Numeric. Relative effectiveness of reactive response compared to preemptive vaccination (0 = no effect, 1 = equally effective).
#' @param f Numeric. Vaccine supply availability as a fraction of the total population (e.g., `f = 1` = full coverage, `f = 0.1` = 10% coverage).
#'
#' @return Numeric. Total expected cost under the mixed intervention strategy.
#' @examples
#' risk <- generate_risk_profile(100)
#' compute_total_cost(alpha = 0.3, risk_profile = risk, C = 10, D = 1000, r = 0.5, f = 0.2)
#' @export
compute_total_cost <- function(alpha, risk_profile, C, D, r, f = 1) {
  n <- length(risk_profile)

  # Calculate number of regions that can be preemptively vaccinated,
  # constrained by vaccine supply fraction f
  k <- floor(alpha * f * n)

  # Preemptive cost: vaccinate k regions at cost C each
  cost_preemptive <- k * C

  # Reactive cost: sum of remaining outbreak risks × (1 - r) × D
  cost_reactive <- sum(risk_profile[(k + 1):n]) * (1 - r) * D

  # Total cost = preemptive + reactive
  cost_preemptive + cost_reactive
}

#' Optimize preemptive vaccination allocation (alpha) to minimize total cost
#'
#' This function identifies the optimal fraction `alpha` of the population to
#' preemptively vaccinate (subject to vaccine availability `f`) to minimize the
#' total expected cost of a mixed preemptive–reactive strategy.
#'
#' The optimization can be done via:
#' - `"optimize"`: continuous search using `stats::optimize()` (default).
#' - `"grid"`: discrete grid search over alpha in [0, 1] using specified step size.
#'
#' @param risk_profile Numeric vector. Risk scores for each region, sorted in decreasing order.
#' @param C Numeric. Cost of preemptive vaccination per region.
#' @param D Numeric. Cost of an unmitigated outbreak per region.
#' @param r Numeric. Relative effectiveness of reactive response compared to preemptive vaccination.
#' @param f Numeric. Vaccine supply availability as a fraction of the population (e.g., `f = 0.2` means 20% of population can be preemptively vaccinated).
#' @param method Character. Optimization method: `"optimize"` (default) or `"grid"`.
#' @param by Numeric. Grid step size for alpha (used only if `method = "grid"`). Default is 0.01.
#'
#' @return A list with:
#'   \item{alpha_opt}{Optimal alpha minimizing total cost}
#'   \item{cost_opt}{Minimum total cost achieved}
#' @examples
#' risk <- generate_risk_profile(100)
#' optimize_alpha(risk, C = 0.2, D = 1, r = 0.5, f = 0.2)
#' optimize_alpha(risk, C = 1, D = 1, r = 0.5, f = 0.2, method = "grid", by = 0.05)
#' @export
optimize_alpha <- function(risk_profile, C, D, r, f = 1, method = "optimize", by = 0.01) {
  # Define cost function to evaluate for a given alpha
  cost_fn <- function(alpha) {
    compute_total_cost(alpha, risk_profile, C, D, r, f)
  }

  if (method == "optimize") {
    # Use continuous optimization
    opt <- optimize(cost_fn, interval = c(0, 1), tol = 1e-4)
    return(list(alpha_opt = opt$minimum, cost_opt = opt$objective))
  }

  if (method == "grid") {
    # Perform discrete grid search over alpha
    alpha_grid <- seq(0, 1, by = by)
    costs <- sapply(alpha_grid, cost_fn)
    i_min <- which.min(costs)
    return(list(alpha_opt = alpha_grid[i_min], cost_opt = costs[i_min]))
  }

  stop("Invalid method. Use 'optimize' or 'grid'.")
}


#' Generate a noisy signal or rank preserving a target correlation
#'
#' This function perturbs a risk signal `p` to simulate imperfect prediction with
#' a specified correlation to the true signal. You can generate:
#'
#' - A continuous noisy version with **Pearson correlation** ≈ `target_rho`, or
#' - A permuted ranking with **Kendall's tau** ≈ `target_tau`.
#'
#' ## If `target_rho` is specified
#'
#' The function constructs:
#'
#' ```
#' z_hat <- rho * p_z + sqrt(1 - rho^2) * noise
#' ```
#'
#' where:
#' - `p_z` is a standardized version of the input signal (normal scores)
#' - `noise` is an independent standard normal vector
#' - `rho` is the target Pearson correlation
#'
#' ### Justification
#' The linear combination:
#'
#' $$
#' Z_{\text{hat}} = \rho Z_{\text{signal}} + \sqrt{1 - \rho^2} Z_{\text{noise}}
#' $$
#'
#' ensures that:
#'
#' - $\operatorname{Var}(Z_{\text{hat}}) = 1$
#' - $\operatorname{Cov}(Z_{\text{hat}}, Z_{\text{signal}}) = \rho$
#' - $\Rightarrow \operatorname{Cor}(Z_{\text{hat}}, Z_{\text{signal}}) = \rho$
#'
#' This is especially useful for modeling imperfect risk estimation with Gaussian-like noise.
#'
#' ## If `target_tau` is specified
#'
#' The function simulates noisy rankings by randomly swapping pairs of indices from a perfectly ordered
#' ranking of `p` with a probability proportional to $(1 - \tau)$.
#'
#' ### Justification
#' Kendall's tau is defined as:
#'
#' $$
#' \tau = \frac{\text{\# concordant pairs} - \text{\# discordant pairs}}{\binom{n}{2}}
#' $$
#'
#' Generate a Sequence with a Target Correlation
#'
#' Generates a new sequence from an input vector `p` that has a specified
#' Pearson's correlation (`rho`) or Kendall's rank correlation (`tau`).
#'
#' @details
#' This function has been revised to use a robust and efficient method based on a
#' Gaussian copula. This approach correctly generates a new vector that not only
#' meets the target correlation but also preserves the empirical distribution of
#' the original vector `p`.
#'
#' The core logic is as follows:
#' 1.  An underlying Pearson's correlation `rho` is determined. If `target_tau` is
#'     provided, it's converted using `rho = sin(target_tau * pi / 2)`. If
#'     `target_rho` is provided, it's used directly.
#' 2.  The input vector `p` is converted to normal scores.
#' 3.  A new, correlated normal score vector is generated.
#' 4.  This new vector is transformed back to the empirical distribution of `p`.
#'
#' This method replaces the previous, less reliable pairwise swapping algorithm
#' for `target_tau` and completes the logic for `target_rho`. The function now
#' also supports negative correlations.
#'
#' @param p Numeric vector. The original sequence of values.
#' @param target_rho Numeric in [-1, 1]. The target Pearson correlation between
#'   the output and `p`. Set this OR `target_tau`.
#' @param target_tau Numeric in [-1, 1]. The target Kendall's tau between the
#'   output and `p`. Set this OR `target_rho`.
#' @param seed Optional integer. A random seed for reproducibility.
#'
#' @return A numeric vector of the same length as `p`. The values in this
#'   vector are drawn from the same empirical distribution as `p` and will have
#'   a correlation with `p` that is approximately the target.
#'
#' @examples
#' set.seed(123) # for reproducibility
#' p <- 1:100
#'
#' # --- Example 1: Target Pearson's rho ---
#' noisy_rho <- generate_noisy_rank(p, target_rho = 0.8)
#' calculated_rho <- cor(noisy_rho, p, method = "pearson")
#' print(paste("Target Rho:", 0.8, "| Calculated Rho:", round(calculated_rho, 4)))
#'
#' plot(p, noisy_rho, main = paste("Pearson's r =", round(calculated_rho, 3)),
#'      xlab = "Original Sequence (p)", ylab = "Generated Sequence (y_rho)", pch=19)
#'
#' # --- Example 2: Target Kendall's tau ---
#' noisy_tau <- generate_noisy_rank(p, target_tau = 0.6)
#' calculated_tau <- cor(noisy_tau, p, method = "kendall")
#' print(paste("Target Tau:", 0.6, "| Calculated Tau:", round(calculated_tau, 4)))
#'
#' plot(p, noisy_tau, main = paste("Kendall's Tau =", round(calculated_tau, 3)),
#'      xlab = "Original Sequence (p)", ylab = "Generated Sequence (y_tau)", pch=19, col="blue")
#'
#' # --- Example 3: Negative Correlation ---
#' noisy_neg_tau <- generate_noisy_rank(p, target_tau = -0.5)
#' calculated_neg_tau <- cor(noisy_neg_tau, p, method = "kendall")
#' print(paste("Target Tau:", -0.5, "| Calculated Tau:", round(calculated_neg_tau, 4)))
#'
#' plot(p, noisy_neg_tau, main = paste("Kendall's Tau =", round(calculated_neg_tau, 3)),
#'      xlab = "Original Sequence (p)", ylab = "Generated Sequence (y_neg_tau)", pch=19, col="red")
#'
#' @export
generate_noisy_rank <- function(p, target_rho = NULL, target_tau = NULL, seed = NULL) {
  # --- 1. Input Validation and Setup ---
  if (!is.null(seed)) set.seed(seed)

  if (!is.null(target_rho) && !is.null(target_tau)) {
    stop("Specify only one of `target_rho` or `target_tau`, not both.")
  }
  if (is.null(target_rho) && is.null(target_tau)) {
    stop("You must specify either `target_rho` or `target_tau`.")
  }
  if (!is.numeric(p)) stop("'p' must be a numeric vector.")
  n <- length(p)
  if (n < 2) stop("'p' must have at least two elements.")

  # --- 2. Determine the underlying Pearson correlation `rho` needed ---
  if (!is.null(target_rho)) {
    if (abs(target_rho) > 1) stop("`target_rho` must be between -1 and 1.")
    rho <- target_rho
  } else { # target_tau must be non-NULL
    if (abs(target_tau) > 1) stop("`target_tau` must be between -1 and 1.")
    rho <- sin(target_tau * pi / 2)
  }

  # --- 3. Handle Edge Cases for Perfect Correlation ---
  if (rho == 1) {
    return(sort(p))
  }
  if (rho == -1) {
    return(rev(sort(p)))
  }

  # --- 4. Core Generation Logic (Gaussian Copula) ---

  # a. Transform p to normal scores
  p_ranks <- rank(p, ties.method = "random")
  p_uniform <- p_ranks / (n + 1)
  z_p <- qnorm(p_uniform)

  # b. Generate the second, correlated normal vector
  error_term <- rnorm(n)
  z_y <- rho * z_p + sqrt(1 - rho^2) * error_term

  # c. Transform the new normal vector back to the empirical distribution of p
  u_y <- pnorm(z_y)
  # Using quantile() with the original `p` maps the new ranks onto the original values
  y <- quantile(p, probs = u_y, type = 1, names = FALSE)

  return(y)
}
