library(tidyverse)
ts <- readr::read_rds("data/outbreak_files_Oct2024/time_series_outbreak_extraction.rds")

# Admin level via :: count
ts <- ts %>% mutate(admin_level = stringr::str_count(location, "::"))

# Filter admin 2 & 3 (2 or 3 :: separators)
ts23 <- ts %>% filter(admin_level %in% c(2, 3))
cat("Rows at admin 2&3:", nrow(ts23), "\n")
cat("Unique locations:", n_distinct(ts23$location), "\n")

# Examine outbreak_number and year columns
cat("\nSample of outbreak_number:\n")
ts23 %>% select(location, TL, TR, year, outbreak_number) %>% head(15) %>% print()

# How many distinct outbreaks per location?
cat("\nOutbreaks per location summary:\n")
ts23 %>%
  group_by(location) %>%
  summarise(n_outbreaks = n_distinct(outbreak_number), .groups="drop") %>%
  summary() %>% print()

# Year range
cat("\nYear range:", range(ts23$year), "\n")

# Distinct outbreak events: define by location + outbreak_number
outbreak_events <- ts23 %>%
  group_by(location, outbreak_number) %>%
  summarise(
    year_start = year(min(TL)),
    year_end   = year(max(TR)),
    duration_weeks = as.numeric(difftime(max(TR), min(TL), units="weeks")),
    .groups = "drop"
  )
cat("\nDuration of outbreaks (weeks):\n")
summary(outbreak_events$duration_weeks) %>% print()
cat("\nOutbreaks spanning multiple years:", sum(outbreak_events$year_start != outbreak_events$year_end), "\n")
