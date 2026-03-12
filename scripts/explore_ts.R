ts <- readr::read_rds("data/outbreak_files_Oct2024/time_series_outbreak_extraction.rds")
cat("TL class:", class(ts$TL), "\n")
cat("TR class:", class(ts$TR), "\n")
cat("Sample TL:\n"); print(head(ts$TL, 5))
cat("\nSample locations:\n"); print(head(unique(ts$location), 20))
cat("\nSpatial scale table:\n"); print(table(ts$spatial_scale))
# Count :: separators in location
n_colons <- stringr::str_count(ts$location, "::")
cat("\nAdmin level (:: count) distribution:\n"); print(table(n_colons))
# Show admin2/3 examples
cat("\nAdmin 2 examples:\n"); print(head(ts$location[n_colons == 2], 10))
cat("\nAdmin 3 examples:\n"); print(head(ts$location[n_colons == 3], 10))
