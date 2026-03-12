library(tidyverse)

ts_outbreak <- readr::read_rds("data/outbreak_files_Oct2024/time_series_outbreak_extraction.rds")
obs_years <- 14

ts23 <- ts_outbreak %>%
  mutate(admin_level = str_count(location, "::")) %>%
  filter(admin_level %in% c(2, 3))

outbreak_events <- ts23 %>%
  group_by(location, outbreak_number) %>%
  summarise(
    year_start = year(min(TL)), year_end = year(max(TR)),
    .groups = "drop"
  )

all_years <- 2010:2023
all_locs  <- unique(ts23$location)

# location × year: count outbreak events active in each calendar year
loc_year <- outbreak_events %>%
  mutate(years = map2(year_start, pmin(year_end, 2023), seq)) %>%
  unnest(years) %>%
  rename(year = years) %>%
  filter(year %in% all_years) %>%
  count(location, year, name = "n_active")

# Full grid (786 locs × 14 years), zeros for unobserved
loc_year_full <- expand_grid(location = all_locs, year = all_years) %>%
  left_join(loc_year, by = c("location", "year")) %>%
  replace_na(list(n_active = 0))

cat("loc_year_full:", nrow(loc_year_full), "rows\n")
cat("Non-zero entries:", sum(loc_year_full$n_active > 0), "\n")

# Cross-validation function
rho_cv <- function(T_val, k_val) {
  valid_splits <- all_years[all_years >= (min(all_years) + T_val - 1) &
                             (all_years + k_val) <= max(all_years)]
  if (length(valid_splits) == 0) return(NULL)

  rho_vals <- map_dbl(valid_splits, function(s) {
    train_yrs <- (s - T_val + 1):s
    test_yrs  <- (s + 1):(s + k_val)

    train_score <- loc_year_full %>%
      filter(year %in% train_yrs) %>%
      group_by(location) %>%
      summarise(score = sum(n_active), .groups = "drop")

    test_out <- loc_year_full %>%
      filter(year %in% test_yrs) %>%
      group_by(location) %>%
      summarise(obs = as.integer(sum(n_active) > 0), .groups = "drop")

    combined <- inner_join(train_score, test_out, by = "location")
    cor(combined$score, combined$obs, method = "spearman")
  })

  tibble(T = T_val, k = k_val,
         rho      = mean(rho_vals, na.rm = TRUE),
         rho_lo   = quantile(rho_vals, 0.25, na.rm = TRUE),
         rho_hi   = quantile(rho_vals, 0.75, na.rm = TRUE),
         n_splits = length(rho_vals))
}

cv_grid <- expand_grid(T_val = 1:13, k_val = c(1, 2, 3, 5)) %>%
  filter(T_val + k_val <= 14)

cv_results <- pmap_dfr(cv_grid, rho_cv)

print(cv_results, n = 40)

# Quick plot
p <- cv_results %>%
  mutate(k = factor(paste0("k = ", k, " yr"), levels = paste0("k = ", c(1,2,3,5), " yr"))) %>%
  ggplot(aes(x = T, y = rho, colour = k, fill = k)) +
  geom_ribbon(aes(ymin = rho_lo, ymax = rho_hi), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:13) +
  labs(x = "Training window (years)", y = expression(rho ~ "(Spearman)"),
       colour = "Forecast", fill = "Forecast") +
  theme_minimal()

ggsave("scripts/test_rho_output.png", p, width = 7, height = 4, dpi = 100)
cat("Plot saved.\n")
