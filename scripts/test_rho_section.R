library(tidyverse)

ts_outbreak <- readr::read_rds("data/outbreak_files_Oct2024/time_series_outbreak_extraction.rds")

ts23 <- ts_outbreak %>%
  mutate(admin_level = str_count(location, "::")) %>%
  filter(admin_level %in% c(2, 3))

outbreak_events <- ts23 %>%
  group_by(location, outbreak_number) %>%
  summarise(year_start = year(min(TL)), year_end = year(max(TR)), .groups = "drop")

all_years <- 2010:2023
all_locs  <- unique(ts23$location)

loc_year <- outbreak_events %>%
  mutate(years = map2(year_start, pmin(year_end, 2023L), seq)) %>%
  unnest(years) %>%
  rename(year = years) %>%
  filter(year %in% all_years) %>%
  count(location, year, name = "n_active")

loc_year_full <- expand_grid(location = all_locs, year = all_years) %>%
  left_join(loc_year, by = c("location", "year")) %>%
  replace_na(list(n_active = 0L))

rho_cv <- function(T_val, k_val) {
  valid_splits <- all_years[all_years >= (min(all_years) + T_val - 1L) &
                              (all_years + k_val) <= max(all_years)]
  if (length(valid_splits) == 0L) return(NULL)
  rho_vals <- map_dbl(valid_splits, function(s) {
    train_score <- loc_year_full %>%
      filter(year %in% seq(s - T_val + 1L, s)) %>%
      group_by(location) %>% summarise(score = sum(n_active), .groups = "drop")
    test_out <- loc_year_full %>%
      filter(year %in% seq(s + 1L, s + k_val)) %>%
      group_by(location) %>% summarise(obs = as.integer(sum(n_active) > 0L), .groups = "drop")
    combined <- inner_join(train_score, test_out, by = "location")
    cor(combined$score, combined$obs, method = "spearman")
  })
  tibble(T=T_val, k=k_val, rho=mean(rho_vals,na.rm=TRUE),
         rho_lo=quantile(rho_vals,0.25,na.rm=TRUE),
         rho_hi=quantile(rho_vals,0.75,na.rm=TRUE), n_splits=length(rho_vals))
}

cv_grid    <- expand_grid(T_val = 1:13, k_val = c(1, 2, 3, 5)) %>% filter(T_val + k_val <= 14L)
cv_results <- pmap_dfr(cv_grid, rho_cv)

k_levels <- paste0("k = ", c(1,2,3,5), ifelse(c(1,2,3,5)==1," year"," years"))
p <- cv_results %>%
  mutate(k_label = factor(paste0("k = ", k, ifelse(k==1," year"," years")), levels=k_levels)) %>%
  ggplot(aes(x=T, y=rho, colour=k_label, fill=k_label)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
  geom_ribbon(aes(ymin=rho_lo, ymax=rho_hi), alpha=0.18, colour=NA) +
  geom_line(linewidth=1) + geom_point(size=2.5) +
  scale_colour_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") +
  labs(x="Training window (years)", y="rho (Spearman)") +
  theme_minimal()
ggsave("scripts/test_rho_section.png", p, width=7, height=4, dpi=100)
cat("Done. Max rho:", round(max(cv_results$rho),3), "at T=", cv_results$T[which.max(cv_results$rho)],
    "k=", cv_results$k[which.max(cv_results$rho)], "\n")
