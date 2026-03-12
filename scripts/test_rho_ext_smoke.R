library(tidyverse); library(sf)
ts <- readr::read_rds("data/outbreak_files_Oct2024/time_series_outbreak_extraction.rds")
ts23 <- ts %>% mutate(admin_level = str_count(location, "::")) %>% filter(admin_level %in% c(2,3))
ob <- ts23 %>% group_by(location, outbreak_number) %>%
  summarise(year_start = year(min(TL)), year_end = year(max(TR)), .groups = "drop")
all_years <- 2010:2023; all_locs <- sort(unique(ts23$location))
lye <- ob %>% filter(year_end %in% all_years) %>%
  count(location, year_end, name = "n_ended") %>% rename(year = year_end) %>%
  right_join(expand_grid(location = all_locs, year = all_years), by = c("location","year")) %>%
  replace_na(list(n_ended = 0L))
cat("loc_year_end_full rows:", nrow(lye), "non-zero:", sum(lye$n_ended > 0), "\n")

shp <- readr::read_rds("data/outbreak_files_Oct2024/outbreak_shapefiles.rds")
loc_id_map <- ts23 %>% distinct(location_period_id, location) %>%
  group_by(location) %>% slice_max(as.integer(location_period_id), n=1, with_ties=FALSE) %>% ungroup()
shp_a23 <- shp %>% left_join(loc_id_map, by=c("lctn_pr"="location_period_id")) %>%
  filter(!is.na(location), str_count(location, "::") %in% c(2,3))
coords <- shp_a23 %>% st_make_valid() %>% st_point_on_surface() %>% st_coordinates()
cat("coords rows:", nrow(coords), "  OK\n")
