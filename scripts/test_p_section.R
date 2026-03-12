library(tidyverse)
library(MASS)

ts_outbreak <- readr::read_rds("data/outbreak_files_Oct2024/time_series_outbreak_extraction.rds")

obs_years <- 2023 - 2010 + 1  # 14 years

ts23 <- ts_outbreak %>%
  mutate(admin_level = str_count(location, "::")) %>%
  filter(admin_level %in% c(2, 3))

outbreak_events <- ts23 %>%
  group_by(location, outbreak_number) %>%
  summarise(
    year_start     = year(min(TL)),
    year_end       = year(max(TR)),
    duration_weeks = as.numeric(difftime(max(TR), min(TL), units = "weeks")),
    .groups = "drop"
  )

location_p <- ts23 %>%
  group_by(location) %>%
  summarise(n_outbreaks = n_distinct(outbreak_number), .groups = "drop") %>%
  mutate(
    lambda_hat = n_outbreaks / obs_years,
    p_annual   = 1 - exp(-lambda_hat)
  )

beta_fit <- fitdistr(location_p$p_annual, "beta",
                     start = list(shape1 = 1, shape2 = 5))
sh1 <- beta_fit$estimate["shape1"]
sh2 <- beta_fit$estimate["shape2"]

cat(sprintf("Beta fit: shape1 = %.3f, shape2 = %.3f\n", sh1, sh2))
cat(sprintf("Implied mean p = %.3f, median p = %.3f\n",
            sh1 / (sh1 + sh2),
            qbeta(0.5, sh1, sh2)))
cat(sprintf("Locations: %d | Outbreaks: %d | Median outbreaks/location: %.1f\n",
            nrow(location_p),
            sum(location_p$n_outbreaks),
            median(location_p$n_outbreaks)))
cat(sprintf("Median outbreak duration: %.0f weeks\n", median(outbreak_events$duration_weeks)))
cat(sprintf("IQR: %.0f-%.0f weeks\n",
            quantile(outbreak_events$duration_weeks, 0.25),
            quantile(outbreak_events$duration_weeks, 0.75)))
cat("Done.\n")
