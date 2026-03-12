library(tidyverse)
library(sf)

ts <- readr::read_rds("data/outbreak_files_Oct2024/time_series_outbreak_extraction.rds")

# How many location_period_ids per location?
loc_per_id <- ts %>%
  distinct(location, location_period_id) %>%
  count(location, name = "n_period_ids")

cat("Locations with >1 period ID:\n")
loc_per_id %>% filter(n_period_ids > 1) %>% print()
cat("Total unique locations:", nrow(loc_per_id), "\n")
cat("Total unique location_period_ids:", n_distinct(ts$location_period_id), "\n")

# Admin level distribution of shapefile-matched locations
ts_loc_map <- ts %>%
  distinct(location_period_id, location) %>%
  mutate(admin_level = str_count(location, "::"))

cat("\nAdmin level distribution in ts:\n")
print(table(ts_loc_map$admin_level))
