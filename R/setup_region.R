#' Reactive vs. Preemptive Vaccination Simulation Functions
#'
#' This script contains functions to simulate outbreaks across subpopulations and evaluate
#' the relative impact of reactive and preemptive vaccination strategies. Outbreaks follow
#' tent-shaped incidence curves. Vaccination modifies the incidence based on coverage, efficacy,
#' and timing of response.
#'
#' Author: Jong-Hoon Kim
#' Date: 2025-06-27

#' Generate a list of synthetic regions with population and outbreak probabilities
#'
#' @param n_region Number of regions to generate
#' @param pop_range Range of population sizes (min, max)
#' @param prob_range Range of outbreak probabilities
#'
#' @return A data frame with one row per region and columns: id, pop_size, outbreak_prob
setup_region <- function(n_region = 100,
                         pop_range = c(1000, 10000),
                         prob_range = c(0.1, 0.9)) {
  tibble::tibble(
    id = 1:n_region,
    pop_size = sample(seq(pop_range[1], pop_range[2]), n_region, replace = TRUE),
    outbreak_prob = runif(n_region, prob_range[1], prob_range[2])
  )
}

#' Generate tent-shaped outbreaks for a list of regions
#'
#' @param regions A data frame from setup_region()
#' @param time_grid Vector of time points for outbreak curves
#'
#' @return A data frame with outbreak curves for each region that experiences an outbreak
simulate_outbreaks <- function(regions, time_grid = seq(0, 100, by = 1)) {
  purrr::map_dfr(1:nrow(regions), function(i) {
    region <- regions[i, ]
    has_outbreak <- runif(1) < region$outbreak_prob
    if (has_outbreak) {
      curve <- simulate_tent_curve(time_grid = time_grid)
      curve$region_id <- region$id
      curve$pop_size <- region$pop_size
      curve$outbreak_prob <- region$outbreak_prob
      return(curve)
    } else {
      return(NULL)
    }
  })
}

library(tibble)
library(dplyr)
library(purrr)

#' Simulate a tent-shaped outbreak incidence curve
#'
#' @param t0 Start time of outbreak
#' @param tp Peak time (optional, randomly drawn if NULL)
#' @param t1 End time (optional)
#' @param r0 Baseline incidence (default 0)
#' @param rp Peak incidence (optional)
#' @param duration_range Range for duration from peak to end
#' @param peak_time_range Range to draw peak timing
#' @param peak_incidence_range Range to draw peak incidence (per capita)
#' @param time_grid Time vector for simulation
#'
#' @return A tibble with columns: time, incidence, t0, tp, t1, r0, rp
simulate_tent_curve <- function(t0 = 0,
                                tp = NULL,
                                t1 = NULL,
                                r0 = 0,
                                rp = NULL,
                                duration_range = c(15, 60),
                                peak_time_range = c(3, 20),
                                peak_incidence_range = c(0.005, 0.05),
                                time_grid = seq(0, 100, by = 1)) {

  if (is.null(tp)) tp <- runif(1, peak_time_range[1], peak_time_range[2])
  if (is.null(t1)) t1 <- tp + runif(1, duration_range[1], duration_range[2])
  if (is.null(rp)) rp <- runif(1, peak_incidence_range[1], peak_incidence_range[2])

  incidence <- numeric(length(time_grid))
  for (i in seq_along(time_grid)) {
    t <- time_grid[i]
    if (t <= t0 || t >= t1) {
      incidence[i] <- r0
    } else if (t > t0 && t <= tp) {
      incidence[i] <- r0 + (rp - r0) * (t - t0) / (tp - t0)
    } else if (t > tp && t < t1) {
      incidence[i] <- rp - (rp - r0) * (t - tp) / (t1 - tp)
    }
  }
  tibble(time = time_grid, incidence = incidence, t0 = t0, tp = tp, t1 = t1, r0 = r0, rp = rp)
}

#' Assign a vaccination strategy to an outbreak
#'
#' @param outbreak_df A data frame with outbreak incidence
#' @param strategy One of 'none', 'preemptive', or 'reactive'
#' @param coverage Vaccine coverage (0-1)
#' @param efficacy Vaccine efficacy (0-1)
#' @param reactive_threshold Cumulative incidence threshold to trigger reactive vaccination
#' @param delay Delay between trigger and vaccine start
#'
#' @return Modified outbreak_df with `vaccinated` and `vaccine_start` columns
assign_vaccination_strategy <- function(outbreak_df,
                                        strategy = c("none", "preemptive", "reactive"),
                                        coverage = 0.7,
                                        efficacy = 0.8,
                                        reactive_threshold = 0.01,
                                        delay = 7) {
  strategy <- match.arg(strategy)
  outbreak_df$cum_incidence <- cumsum(outbreak_df$incidence)

  if (strategy == "none") {
    outbreak_df$vaccinated <- 0
    outbreak_df$vaccine_start <- Inf
  }

  if (strategy == "preemptive") {
    outbreak_df$vaccinated <- coverage * efficacy
    outbreak_df$vaccine_start <- min(outbreak_df$time) - delay
  }

  if (strategy == "reactive") {
    trigger_time <- outbreak_df$time[which(outbreak_df$cum_incidence >= reactive_threshold)[1]]
    if (is.na(trigger_time)) {
      outbreak_df$vaccine_start <- Inf
      outbreak_df$vaccinated <- 0
    } else {
      vaccine_start <- trigger_time + delay
      outbreak_df$vaccine_start <- vaccine_start
      outbreak_df$vaccinated <- ifelse(outbreak_df$time >= vaccine_start, coverage * efficacy, 0)
    }
  }
  return(outbreak_df)
}

#' Compute disease outcomes including cases and deaths with and without vaccine
#'
#' @param outbreak_df A data frame including `incidence` and `vaccinated`
#' @param cfr Case-fatality ratio
#'
#' @return Data frame with additional columns for cases and deaths averted
simulate_disease_impact <- function(outbreak_df, cfr = 0.01) {
  with_vaccine <- outbreak_df$incidence * (1 - outbreak_df$vaccinated)
  without_vaccine <- outbreak_df$incidence

  outbreak_df$cases_with_vaccine <- with_vaccine
  outbreak_df$cases_without_vaccine <- without_vaccine
  outbreak_df$cases_averted <- without_vaccine - with_vaccine

  outbreak_df$deaths_with_vaccine <- cfr * outbreak_df$cases_with_vaccine
  outbreak_df$deaths_without_vaccine <- cfr * outbreak_df$cases_without_vaccine
  outbreak_df$deaths_averted <- outbreak_df$deaths_without_vaccine - outbreak_df$deaths_with_vaccine

  return(outbreak_df)
}

#' Run a full simulation across multiple regions and strategies
#'
#' @param n_region Number of regions to simulate
#' @param strategies Vector of strategies (e.g., c("none", "preemptive", "reactive"))
#' @param time_grid Time vector for simulation
#' @param ... Additional arguments passed to strategy assignment
#'
#' @return A combined data frame of all results
run_scenario <- function(n_region = 100,
                         strategies = c("none", "preemptive", "reactive"),
                         time_grid = seq(0, 100, by = 1),
                         ...) {
  results <- purrr::map_dfr(1:n_region, function(i) {
    outbreak <- simulate_tent_curve(time_grid = time_grid)
    purrr::map_dfr(strategies, function(strategy) {
      vaccinated <- assign_vaccination_strategy(outbreak, strategy = strategy, ...)
      outcome <- simulate_disease_impact(vaccinated)
      outcome$region_id <- i
      outcome$strategy <- strategy
      outcome
    })
  })
  return(results)
}
