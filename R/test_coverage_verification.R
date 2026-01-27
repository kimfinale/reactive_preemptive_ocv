source("R/load_all.R")

# Verification:
# If N=1000, coverage=0.5, vacc_duration=5, and beta=0 (no infection),
# S should decrease from 1000 to 500.

print("Testing Coverage Logic...")

N_test <- 1000
coverage_test <- 0.5
dur_test <- 5

tryCatch(
    {
        out <- run_sir_vaccination(
            time_step = 0.01, duration = 10, N = N_test,
            I0 = 0, R0 = 0, gamma = 0, # No infection
            coverage = coverage_test, vacc_start = 0, vacc_duration = dur_test, ve = 1
        )

        print("Structure of out:")
        str(out)

        final_S <- tail(out$S, 1)
        final_S <- tail(out$S, 1)

        print(paste("Final S:", final_S))
        print(paste("Target S:", N_test * (1 - coverage_test)))

        # Expected S is approx N * (1 - coverage) = 500
        if (abs(final_S - 500) < 1) {
            print("Coverage verification passed: ~500 S remaining")
        } else {
            print("Coverage verification failed")
            quit(status = 1)
        }
    },
    error = function(e) {
        print(paste("Error:", e$message))
        quit(status = 1)
    }
)
