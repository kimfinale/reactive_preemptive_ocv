library(tidyverse)
library(MASS)

ts_outbreak <- readr::read_rds("data/outbreak_files_Oct2024/time_series_outbreak_extraction.rds")
obs_years <- 14

ts23 <- ts_outbreak %>%
  mutate(admin_level = str_count(location, "::")) %>%
  filter(admin_level %in% c(2, 3))

location_p <- ts23 %>%
  group_by(location) %>%
  summarise(n_outbreaks = n_distinct(outbreak_number), .groups = "drop") %>%
  mutate(lambda_hat = n_outbreaks / obs_years,
         p_annual   = 1 - exp(-lambda_hat))

outbreak_events <- ts23 %>%
  group_by(location, outbreak_number) %>%
  summarise(duration_weeks = as.numeric(difftime(max(TR), min(TL), units = "weeks")),
            .groups = "drop")

x <- location_p$p_annual
n <- length(x)

fits <- list(
  Beta        = fitdistr(x, "beta",        start = list(shape1 = 1, shape2 = 5)),
  `Log-normal` = fitdistr(x, "lognormal"),
  Gamma       = fitdistr(x, "gamma",       lower = 1e-6),
  Weibull     = fitdistr(x, "weibull"),
  Exponential = fitdistr(x, "exponential")
)

fit_stats <- imap_dfr(fits, function(f, nm) {
  k  <- length(f$estimate)
  ll <- as.numeric(logLik(f))
  tibble(Distribution = nm, k = k, LogLik = ll,
         AIC = -2*ll + 2*k, BIC = -2*ll + k*log(n),
         Parameters = paste(sprintf("%s=%.3f", names(f$estimate), f$estimate), collapse=", "))
}) %>% arrange(AIC)

print(fit_stats)
cat("\nBest fit:", fit_stats$Distribution[1], "\n")
