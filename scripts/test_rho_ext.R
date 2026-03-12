library(tidyverse)
library(sf)

# ── data ──────────────────────────────────────────────────────────────────────
ts_outbreak <- readr::read_rds("data/outbreak_files_Oct2024/time_series_outbreak_extraction.rds")
shp_all     <- readr::read_rds("data/outbreak_files_Oct2024/outbreak_shapefiles.rds")

ts23 <- ts_outbreak %>%
  mutate(admin_level = str_count(location, "::")) %>%
  filter(admin_level %in% c(2, 3))

outbreak_events <- ts23 %>%
  group_by(location, outbreak_number) %>%
  summarise(year_start = year(min(TL)), year_end = year(max(TR)), .groups = "drop")

all_years <- 2010:2023
all_locs  <- sort(unique(ts23$location))

loc_year <- outbreak_events %>%
  mutate(yrs = map2(year_start, pmin(year_end, 2023L), seq)) %>%
  unnest(yrs) %>% rename(year = yrs) %>%
  filter(year %in% all_years) %>%
  count(location, year, name = "n_active")

loc_year_full <- expand_grid(location = all_locs, year = all_years) %>%
  left_join(loc_year, by = c("location", "year")) %>%
  replace_na(list(n_active = 0L))

# (1) Completed-only matrix: count endings, not active years
loc_year_end_full <- outbreak_events %>%
  filter(year_end %in% all_years) %>%
  count(location, year_end, name = "n_ended") %>%
  rename(year = year_end) %>%
  {right_join(expand_grid(location = all_locs, year = all_years), .,
              by = c("location", "year"))} %>%
  replace_na(list(n_ended = 0L))

cat("Ongoing outbreaks (year_end > 2023):",
    sum(outbreak_events$year_end > 2023), "\n")
cat("Total outbreaks:", nrow(outbreak_events), "\n")
cat("Outbreaks ending in 2023:", sum(outbreak_events$year_end == 2023), "\n")

# (2) Spatial neighbor list from centroids
loc_id_map <- ts23 %>%
  distinct(location_period_id, location) %>%
  group_by(location) %>%
  slice_max(as.integer(location_period_id), n = 1, with_ties = FALSE) %>%
  ungroup()

shp_admin23 <- shp_all %>%
  left_join(loc_id_map, by = c("lctn_pr" = "location_period_id")) %>%
  filter(!is.na(location), str_count(location, "::") %in% c(2, 3)) %>%
  arrange(match(location, all_locs))

coords_mat <- shp_admin23 %>%
  st_make_valid() %>%
  st_point_on_surface() %>%
  st_coordinates()
rownames(coords_mat) <- shp_admin23$location

dist_mat <- as.matrix(dist(coords_mat))
rownames(dist_mat) <- colnames(dist_mat) <- shp_admin23$location
diag(dist_mat) <- Inf

K_NBR_MAX <- 10L
nbr_list <- setNames(
  lapply(shp_admin23$location,
         function(loc) names(sort(dist_mat[loc, ]))[seq_len(K_NBR_MAX)]),
  shp_admin23$location
)
cat("Neighbor list built for", length(nbr_list), "locations\n")

# ── cross-validation function ─────────────────────────────────────────────────
rho_cv_ext <- function(T_val, k_val, L_val = 0L,
                       excl_ongoing = FALSE, n_nbr = 0L) {
  valid_splits <- all_years[
    all_years >= (min(all_years) + T_val - 1L) &
    (all_years + L_val + k_val) <= max(all_years)
  ]
  if (length(valid_splits) == 0L) return(NULL)

  rho_vals <- map_dbl(valid_splits, function(s) {
    train_yrs <- seq(s - T_val + 1L, s)
    test_yrs  <- seq(s + L_val + 1L, s + L_val + k_val)

    # Training score
    if (excl_ongoing) {
      # Count outbreak *completions* inside the training window (year_end in train_yrs)
      # → outbreaks still ongoing at s are automatically excluded
      train_score <- loc_year_end_full %>%
        filter(year %in% train_yrs) %>%
        group_by(location) %>%
        summarise(score = sum(n_ended), .groups = "drop")
    } else {
      train_score <- loc_year_full %>%
        filter(year %in% train_yrs) %>%
        group_by(location) %>%
        summarise(score = sum(n_active), .groups = "drop")
    }

    # Spatial smoothing: average score with n_nbr nearest neighbours
    if (n_nbr > 0L) {
      sv <- setNames(train_score$score, train_score$location)
      train_score <- train_score %>%
        mutate(score = map_dbl(location, function(loc) {
          nbrs <- nbr_list[[loc]][seq_len(n_nbr)]
          nbrs <- nbrs[nbrs %in% names(sv)]
          mean(c(sv[[loc]], sv[nbrs]))
        }))
    }

    # Test outcome: any outbreak active in test window?
    test_out <- loc_year_full %>%
      filter(year %in% test_yrs) %>%
      group_by(location) %>%
      summarise(obs = as.integer(sum(n_active) > 0L), .groups = "drop")

    combined <- inner_join(train_score, test_out, by = "location")
    cor(combined$score, combined$obs, method = "spearman")
  })

  tibble(T = T_val, k = k_val, L = L_val,
         excl_ongoing = excl_ongoing, n_nbr = n_nbr,
         rho      = mean(rho_vals, na.rm = TRUE),
         rho_lo   = quantile(rho_vals, 0.25, na.rm = TRUE),
         rho_hi   = quantile(rho_vals, 0.75, na.rm = TRUE),
         n_splits = length(rho_vals))
}

# ── run grids ─────────────────────────────────────────────────────────────────
cat("Running ongoing-exclusion grid...\n")
cv_excl <- expand_grid(T_val = 1:13, k_val = c(1,2,3,5),
                        L_val = 0L, excl_ongoing = c(FALSE, TRUE), n_nbr = 0L) %>%
  filter(T_val + k_val <= 14L) %>%
  pmap_dfr(rho_cv_ext)

cat("Running lag grid...\n")
cv_lag <- expand_grid(T_val = 1:10, k_val = c(1,2,3,5),
                       L_val = 0:2, excl_ongoing = TRUE, n_nbr = 0L) %>%
  filter(T_val + L_val + k_val <= 14L) %>%
  pmap_dfr(rho_cv_ext)

cat("Running spatial grid...\n")
cv_nbr <- expand_grid(T_val = 1:10, k_val = c(1,2,3,5),
                       L_val = 0L, excl_ongoing = TRUE, n_nbr = c(0L,3L,5L,10L)) %>%
  filter(T_val + k_val <= 14L) %>%
  pmap_dfr(rho_cv_ext)

cat("Done. Rows:", nrow(cv_excl), nrow(cv_lag), nrow(cv_nbr), "\n")

# ── quick check ───────────────────────────────────────────────────────────────
cat("\nTop rho rows (excl grid):\n")
cv_excl %>% arrange(desc(rho)) %>% head(6) %>% print()
cat("\nTop rho rows (lag grid):\n")
cv_lag  %>% arrange(desc(rho)) %>% head(6) %>% print()
cat("\nTop rho rows (spatial grid):\n")
cv_nbr  %>% arrange(desc(rho)) %>% head(6) %>% print()

# ── plots ──────────────────────────────────────────────────────────────────────
k_levels <- paste0("k = ", c(1,2,3,5), ifelse(c(1,2,3,5)==1," yr"," yrs"))
k_lab <- function(x) factor(paste0("k = ", x, ifelse(x==1," yr"," yrs")), levels=k_levels)

p1 <- cv_excl %>%
  mutate(k_label    = k_lab(k),
         condition  = ifelse(excl_ongoing, "Exclude ongoing", "Include ongoing")) %>%
  ggplot(aes(x=T, y=rho, colour=condition, fill=condition)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
  geom_ribbon(aes(ymin=rho_lo, ymax=rho_hi), alpha=0.15, colour=NA) +
  geom_line(linewidth=0.9) + geom_point(size=2) +
  facet_wrap(~k_label, nrow=2) +
  scale_colour_manual(values=c("steelblue","firebrick")) +
  scale_fill_manual(values=c("steelblue","firebrick")) +
  labs(title="Effect of excluding ongoing outbreaks",
       x="Training window T (years)", y=expression(rho), colour=NULL, fill=NULL) +
  theme_minimal(base_size=11) + theme(legend.position="bottom")

p2 <- cv_lag %>%
  mutate(k_label   = k_lab(k),
         lag_label = factor(paste0("Lag L = ", L, " yr"), levels=paste0("Lag L = ",0:2," yr"))) %>%
  ggplot(aes(x=T, y=rho, colour=lag_label, fill=lag_label)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
  geom_ribbon(aes(ymin=rho_lo, ymax=rho_hi), alpha=0.15, colour=NA) +
  geom_line(linewidth=0.9) + geom_point(size=2) +
  facet_wrap(~k_label, nrow=2) +
  scale_colour_brewer(palette="Dark2") +
  scale_fill_brewer(palette="Dark2") +
  labs(title="Effect of prediction lag (excl. ongoing)",
       x="Training window T (years)", y=expression(rho), colour=NULL, fill=NULL) +
  theme_minimal(base_size=11) + theme(legend.position="bottom")

p3 <- cv_nbr %>%
  mutate(k_label   = k_lab(k),
         nbr_label = factor(paste0("n_nbr = ", n_nbr),
                            levels=paste0("n_nbr = ", c(0,3,5,10)))) %>%
  ggplot(aes(x=T, y=rho, colour=nbr_label, fill=nbr_label)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
  geom_ribbon(aes(ymin=rho_lo, ymax=rho_hi), alpha=0.15, colour=NA) +
  geom_line(linewidth=0.9) + geom_point(size=2) +
  facet_wrap(~k_label, nrow=2) +
  scale_colour_brewer(palette="Set2") +
  scale_fill_brewer(palette="Set2") +
  labs(title="Effect of spatial neighbourhood smoothing (excl. ongoing, L = 0)",
       x="Training window T (years)", y=expression(rho), colour=NULL, fill=NULL) +
  theme_minimal(base_size=11) + theme(legend.position="bottom")

library(patchwork)
ggsave("scripts/rho_ext_excl.png",    p1, width=8, height=6, dpi=100)
ggsave("scripts/rho_ext_lag.png",     p2, width=8, height=6, dpi=100)
ggsave("scripts/rho_ext_spatial.png", p3, width=8, height=6, dpi=100)
cat("Plots saved.\n")
