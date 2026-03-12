library(tidyverse)
library(sf)

shp <- readr::read_rds("data/outbreak_files_Oct2024/outbreak_shapefiles.rds")
ts  <- readr::read_rds("data/outbreak_files_Oct2024/time_series_outbreak_extraction.rds")

# lctn_pr is character or numeric?
cat("lctn_pr class:", class(shp$lctn_pr), "\n")
cat("location_period_id class:", class(ts$location_period_id), "\n")
cat("Sample lctn_pr:", head(shp$lctn_pr, 5), "\n")
cat("Sample location_period_id:", head(ts$location_period_id, 5), "\n")

# Check overlap
shp_ids <- unique(shp$lctn_pr)
ts_ids  <- unique(as.character(ts$location_period_id))
cat("Shapefile IDs:", length(shp_ids), "\n")
cat("Time series IDs:", length(ts_ids), "\n")
cat("Overlap:", sum(shp_ids %in% ts_ids), "\n")

# Get location name for a few shapefile IDs
ts_loc_map <- ts %>%
  distinct(location_period_id, location) %>%
  mutate(location_period_id = as.character(location_period_id))
joined <- shp %>%
  sf::st_drop_geometry() %>%
  left_join(ts_loc_map, by = c("lctn_pr" = "location_period_id"))
cat("\nSample join:\n")
print(head(joined[!is.na(joined$location),], 10))
cat("NAs in location after join:", sum(is.na(joined$location)), "/", nrow(joined), "\n")
