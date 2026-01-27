source("R/load_all.R")

log_file <- "debug_test.log"
sink(log_file)

print("Starting Test...")

N_test <- 1000
coverage_test <- 0.5
dur_test <- 5

print(paste("Config: N", N_test, "Cov", coverage_test, "Dur", dur_test))

# Manually calc expected rho
rho_expected <- -log(1 - coverage_test) / dur_test
print(paste("Expected Rho:", rho_expected))

out <- try(run_sir_vaccination(
    time_step = 0.01, duration = 10, N = N_test,
    I0 = 0, R0 = 0, gamma = 0, # No infection
    coverage = coverage_test, vacc_start = 0, vacc_duration = dur_test, ve = 1
), silent = FALSE)

if (inherits(out, "try-error")) {
    print("Error in run_sir_vaccination")
    print(out)
} else {
    final_S <- tail(out$S, 1)

    print(paste("Final S:", final_S))
    print(paste("Target S:", N_test * (1 - coverage_test)))

    if (abs(final_S - 500) < 1) {
        print("PASS")
    } else {
        print("FAIL")
    }
}

sink()
