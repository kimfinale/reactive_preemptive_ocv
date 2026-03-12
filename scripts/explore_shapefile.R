library(tidyverse)
shp <- readr::read_rds("data/outbreak_files_Oct2024/outbreak_shapefiles.rds")
cat("Class:", class(shp), "\n")
cat("Names:", paste(names(shp), collapse=", "), "\n")
cat("Rows:", nrow(shp), "\n")
# Look at first few rows (non-geometry columns)
if (inherits(shp, "sf")) {
  library(sf)
  cat("CRS:", st_crs(shp)$input, "\n")
  print(head(sf::st_drop_geometry(shp), 10))
} else {
  print(head(shp, 5))
}
