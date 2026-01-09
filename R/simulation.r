# simulation.R
# This file defines all functions used to run simulations in R, inspired by the Python version.

library(digest)      # for hashing
library(readr)       # for reading/writing CSVs
library(dplyr)       # for data manipulation
library(tibble)      # for tidy outputs
library(fs)          # for path management

source("montecarlo.R")
source("kinetics.R")

conduct_simulation <- function(inputs, N_draws) {
  list2env(setNames(inputs, c("wait", "supply", "Q", "tat", "c", "f", "T", "L",
                               "wait_exit", "Q_exit", "f_exit", "L_exit",
                               "kinetics", "testing", "isolation")), envir = environment())

  opts_testing <- list(
    test_regular = test_regular,
    test_post_exposure = test_post_exposure,
    test_post_symptoms = test_post_symptoms
  )
  testing_function <- opts_testing[[testing]]

  opts_kinetics <- list(
    kinetics_test = kinetics_test,
    kinetics_flu = kinetics_flu,
    kinetics_sarscov2_founder_naive = kinetics_sarscov2_founder_naive,
    kinetics_sarscov2_omicron_experienced = kinetics_sarscov2_omicron_experienced,
    kinetics_rsv = kinetics_rsv
  )
  kinetics_function <- opts_kinetics[[kinetics]]

  results <- vector("list", N_draws)
  for (i in seq_len(N_draws)) {
    results[[i]] <- get_sample(wait, supply, Q, tat, c, f, T, L,
                               wait_exit, Q_exit, f_exit, L_exit,
                               kinetics_function, testing_function, isolation)
  }
  return(results)
}

save_simulation <- function(simulation, outputs, output_columns, output_path) {
  df <- as_tibble(do.call(rbind, outputs))
  names(df) <- output_columns
  zip_path <- path(output_path, paste0(simulation, ".zip"))
  write_csv(df, zip_path)
}

summarize_simulation <- function(df) {
  R_0 <- mean(df$I0, na.rm = TRUE)
  R_testing_perfect_isolation <- mean(df$Itest, na.rm = TRUE)
  R_post_exit <- mean(df$Iexit, na.rm = TRUE)
  R_testing <- R_testing_perfect_isolation + R_post_exit
  ascertainment <- compute_ascertainment(df$tDx)
  n_tests_dx <- mean(df$n_tests, na.rm = TRUE)
  n_tests_exit <- mean(df$n_tests_exit, na.rm = TRUE)
  T_isolation <- compute_isolation_time(df$tDx, df$tExit)
  TE <- 1 - R_testing / R_0

  return(c(TE, ascertainment, n_tests_dx, R_0, R_testing, R_post_exit, n_tests_exit, T_isolation))
}

compute_ascertainment <- function(tDx) {
  N <- length(tDx)
  misses <- sum(is.infinite(tDx))
  return(1 - misses / N)
}

compute_isolation_time <- function(tDx, tExit) {
  T_isolation <- tExit - tDx
  T_nonzero <- T_isolation[T_isolation > 0]
  if (length(T_nonzero) == 0) return(0)
  return(mean(T_nonzero, na.rm = TRUE))
}

local_hash <- function(inputs) {
  input_string <- paste(unlist(inputs), collapse = ",")
  hashed <- digest(input_string, algo = "sha256")
  return(hashed)
}

append_simulation_to_log <- function(log_path, inputs, simulation, outcomes) {
  if (simulation != local_hash(inputs)) {
    stop("Hash failure")
  }
  line <- c(simulation, unlist(inputs), outcomes)
  write(paste(line, collapse = ","), file = log_path, append = TRUE)
}

initialize_log <- function(log_path, input_columns, outcome_columns) {
  if (!file_exists(log_path)) {
    header <- c("simulation", input_columns, outcome_columns)
    write(paste(header, collapse = ","), file = log_path)
  }
}

is_simulation_redundant <- function(simulation, output_folder) {
  fname <- path(output_folder, paste0(simulation, ".zip"))
  return(file_exists(fname))
}
