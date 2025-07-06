# Load libraries
library(yaml)
library(tibble)
library(purrr)
library(readr)
library(fs)
library(stringr)

# Custom source files
source("R/simulation.R")  # this should include functions like conduct_simulation(), initialize_log(), etc.

# Set paths
output_path <- path_real(path("..", "output"))
input_path  <- path_real(path("..", "input"))

# Argument parser
args <- commandArgs(trailingOnly = TRUE)
experiment_file <- args[1]

# YAML → list
yaml_path <- path(input_path, experiment_file)
parsed_yaml <- tryCatch(
  yaml::read_yaml(yaml_path),
  error = function(e) stop("YAML parsing error: ", e)
)

# Convert list to variables in global env
for (name in names(parsed_yaml)) {
  value <- parsed_yaml[[name]]
  if (is.list(value) || is.vector(value)) {
    assign(name, value, envir = .GlobalEnv)
  } else if (name %in% c("N_draws", "log_name")) {
    assign(name, value, envir = .GlobalEnv)
  } else {
    assign(name, list(value), envir = .GlobalEnv)
  }
}

# Fixed input and output columns
input_columns <- c("wait", "supply", "Q", "tat", "c", "f", "T", "L",
                   "wait_exit", "Q_exit", "f_exit", "L_exit",
                   "kinetics", "testing", "isolation")

log_columns <- c("tDx", "tExit", "n_tests", "n_tests_exit", "I0", "Itest",
                 "Iexit", "A", "P", "B", "tSx", "m", "M", "first", "last", "infectious_threshold")

outcome_columns <- c("TE", "ascertainment", "n_tests_dx", "R_no_testing",
                     "R_testing", "R_post_exit", "n_tests_exit", "T_isolation")

# Initialize log file
log_path <- path(output_path, paste0(log_name, ".csv"))
initialize_log(log_path, input_columns, outcome_columns)

# Generate all parameter combinations
param_ranges <- list(wait, supply, Q, tat, c, f, T, L,
                     wait_exit, Q_exit, f_exit, L_exit,
                     kinetics, testing, isolation)

param_grid <- expand.grid(param_ranges, stringsAsFactors = FALSE)

# Iterate through combinations
walk(1:nrow(param_grid), function(i) {
  args <- as.list(param_grid[i, ])
  simulation_name <- local_hash(args)

  if (is_simulation_redundant(simulation_name, output_path)) {
    cat("🚀", simulation_name, "\n")
  } else {
    outputs <- conduct_simulation(args, N_draws)
    save_simulation(simulation_name, outputs, log_columns, output_path)
    cat("🧑‍💻", simulation_name, "\n")
  }

  df <- read_csv(path(output_path, paste0(simulation_name, ".zip")))
  outcomes <- summarize_simulation(df)
  cat(sprintf("\tTE=%.4f\n", outcomes[[1]]))
  append_simulation_to_log(log_path, args, simulation_name, outcomes)
})
