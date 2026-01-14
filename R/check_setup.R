# check_setup.R

# Source all R files
r_files <- list.files(".", pattern = "\\.R$", full.names = TRUE)
# filter out check_setup.R to avoid infinite recursion if I source current dir
r_files <- r_files[!grepl("check_setup.R", r_files)]

for (f in r_files) {
    tryCatch(
        {
            source(f)
            cat("Sourced:", f, "\n")
        },
        error = function(e) {
            cat("Error sourcing", f, ":", e$message, "\n")
        }
    )
}

# Test Analytical Functions
cat("\nTesting Analytical Functions...\n")
tryCatch(
    {
        c_pre <- cost_pre_one(p = 0.5, R = 10, r = 0.6, nu = 0.8)
        cat("Cost Pre One:", c_pre, "\n")
    },
    error = function(e) {
        cat("Error in Analytical Function:", e$message, "\n")
    }
)

# Test Simulation Functions
cat("\nTesting Simulation Functions...\n")
tryCatch(
    {
        library(dplyr)
        library(tibble)

        # Setup dummy regions
        regions <- tibble(
            id = 1:5,
            pop_size = 1000,
            outbreak_prob = 0.5,
            Y = c(TRUE, FALSE, TRUE, FALSE, TRUE) # Manually adding Y as expected by assign_outbreak_trajectory
        )

        # Assign trajectory
        outbreaks <- assign_outbreak_trajectory(regions)
        print(head(outbreaks))

        # Assign reactive vaccination
        rv_outbreaks <- assign_reactive_vaccination(
            outbreaks,
            t_now = 50,
            max_targets = 2
        )
        print(head(rv_outbreaks))
    },
    error = function(e) {
        cat("Error in Simulation Function:", e$message, "\n")
    }
)
