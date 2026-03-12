library(tidyverse)
library(sf)
library(MASS)

ts_outbreak <- readr::read_rds("data/outbreak_files_Oct2024/time_series_outbreak_extraction.rds")
obs_years <- 14
k_vals <- c(1, 2, 3, 5)

ts23 <- ts_outbreak %>%
  mutate(admin_level = str_count(location, "::")) %>%
  filter(admin_level %in% c(2, 3))

location_p <- ts23 %>%
  group_by(location) %>%
  summarise(n_outbreaks = n_distinct(outbreak_number), .groups = "drop") %>%
  mutate(lambda_hat = n_outbreaks / obs_years,
         p_annual   = 1 - exp(-lambda_hat))

shp_all <- readr::read_rds("data/outbreak_files_Oct2024/outbreak_shapefiles.rds")

loc_id_map <- ts23 %>%
  distinct(location_period_id, location) %>%
  group_by(location) %>%
  slice_max(as.integer(location_period_id), n = 1, with_ties = FALSE) %>%
  ungroup()

shp_named <- shp_all %>%
  left_join(loc_id_map, by = c("lctn_pr" = "location_period_id"))

shp_admin23 <- shp_named %>%
  filter(!is.na(location), str_count(location, "::") %in% c(2, 3))

cat("shp_admin23 rows:", nrow(shp_admin23), "\n")
cat("Matched locations:", n_distinct(shp_admin23$location), "\n")

shp_maps <- map_dfr(k_vals, function(k_val) {
  shp_admin23 %>%
    left_join(
      location_p %>% transmute(location, p_k = 1 - exp(-k_val * lambda_hat)),
      by = "location"
    ) %>%
    mutate(k_label = factor(
      paste0("k = ", k_val, ifelse(k_val == 1, " year", " years")),
      levels = paste0("k = ", k_vals, ifelse(k_vals == 1, " year", " years"))
    ))
})

cat("shp_maps rows:", nrow(shp_maps), "\n")
cat("p_k range k=1:", range(shp_maps$p_k[shp_maps$k_label == "k = 1 year"], na.rm=TRUE), "\n")
cat("p_k range k=5:", range(shp_maps$p_k[shp_maps$k_label == "k = 5 years"], na.rm=TRUE), "\n")

# Test plot
p <- ggplot() +
  geom_sf(data = shp_all, fill = "grey92", colour = "grey70", linewidth = 0.1) +
  geom_sf(data = shp_maps, aes(fill = p_k), colour = NA) +
  facet_wrap(~k_label, nrow = 2) +
  scale_fill_viridis_c(option = "plasma", na.value = "grey85",
                       labels = scales::label_percent(accuracy = 1),
                       name = "P(>=1 outbreak)") +
  coord_sf(xlim = c(-20, 55), ylim = c(-35, 38), expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), legend.position = "bottom")

ggsave("scripts/test_map_output.png", p, width = 10, height = 8, dpi = 100)
cat("Map saved to scripts/test_map_output.png\n")
