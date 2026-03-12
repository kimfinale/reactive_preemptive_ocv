# scripts/fake_outbreak_timeseries.R
# ─────────────────────────────────────────────────────────────────────────────
# Simulate weekly incidence for 50 regions over 10 years.
# Plots:
#   Top panel    – P(≥1 outbreak in next k years) for k = 1, 2, 5
#   Bottom panel – incidence heatmap (one row per region)
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})

set.seed(42)

N_REG <- 50L
N_YR  <- 10L
WPY   <- 52L                  # weeks per year
N_WK  <- N_YR * WPY           # 520 total weeks
WEEKS <- seq_len(N_WK)

# ── 1.  Inter-outbreak interval per region ────────────────────────────────────
# Log-uniform on [0.8, 5.5] yr gives a right-skewed distribution of outbreak
# counts (many regions with 2-4 outbreaks; a few with >10 over 10 years).
# Sort descending so region 1 = most frequent (top of heatmap).
ival_yr <- sort(
  exp(runif(N_REG, log(0.8), log(5.5))),
  decreasing = TRUE
)

# ── 2.  Generate one region's time series ────────────────────────────────────
gen_ts <- function(reg, ival) {
  ival_wk <- ival * WPY
  ts      <- rpois(N_WK, 1)            # sparse endemic background

  # Poisson-process outbreaks: exponential inter-arrival times
  t0 <- rexp(1, 1 / ival_wk)
  while (t0 < N_WK) {
    i0  <- max(1L, min(N_WK, round(t0)))
    dur <- max(5L, min(24L, round(rlnorm(1, log(11), 0.40))))  # weeks
    pk  <- max(2L, round(dur * rbeta(1, 2, 4)))   # peak early in outbreak
    amp <- rlnorm(1, log(75), 0.55)                # peak case count

    # Normal-shaped epidemic curve, normalised to amplitude
    crv <- dnorm(seq_len(dur), mean = pk, sd = dur / 5)
    crv <- crv / max(crv) * amp

    i1  <- min(i0 + dur - 1L, N_WK)
    len <- i1 - i0 + 1L
    ts[i0:i1] <- ts[i0:i1] + rpois(len, pmax(0.5, crv[seq_len(len)]))

    t0 <- t0 + rexp(1, 1 / ival_wk)
  }
  tibble(region = reg, week = WEEKS, incidence = ts, ival_yr = ival)
}

ts_df <- map2_dfr(seq_len(N_REG), ival_yr, gen_ts)

# ── 3.  Define outbreak: ≥4 consecutive weeks above mean + 2·SD ──────────────
# Step 1: flag weeks that exceed the region-specific threshold
ts_df <- ts_df |>
  group_by(region) |>
  mutate(
    reg_mean     = mean(incidence),
    reg_sd       = sd(incidence),
    above_thresh = incidence > reg_mean + 2 * reg_sd
  ) |>
  ungroup()

# Step 2: within each region, mark a run of above-threshold weeks as an
# outbreak only if it contains ≥ 4 consecutive weeks (rle-based).
mark_outbreaks <- function(above) {
  r   <- rle(above)
  # qualify only TRUE runs of length ≥ 4
  out <- r$values & (r$lengths >= 4L)
  rep(out, r$lengths)
}

ts_df <- ts_df |>
  group_by(region) |>
  mutate(in_outbreak = mark_outbreaks(above_thresh)) |>
  ungroup()

# Count distinct outbreak episodes per region
n_ep <- ts_df |>
  group_by(region) |>
  summarise(
    n_episodes = sum(diff(c(FALSE, in_outbreak)) == 1L),
    ival_yr    = first(ival_yr),
    .groups    = "drop"
  )

cat("Outbreak episodes per region (summary):\n")
print(summary(n_ep$n_episodes))
cat("\nDistribution across 50 regions:\n")
print(table(n_ep$n_episodes))
cat("\n")

# Rank regions by actual observed episodes (descending: rank 1 = most frequent)
# Ties broken by original region index (stable sort).
n_ep <- n_ep |>
  arrange(desc(n_episodes), region) |>
  mutate(r_ord = row_number())

ts_df <- ts_df |>
  left_join(n_ep |> select(region, r_ord, n_episodes), by = "region")

# ── 4.  P(outbreak in next k years) – circular wrapping ──────────────────────
# outbreak matrix: rows = regions, cols = weeks
ob_mat <- matrix(
  ts_df |> arrange(region, week) |> pull(in_outbreak),
  nrow = N_REG, ncol = N_WK, byrow = TRUE
)

# At time t, look ahead k_wk weeks (wrapping back to start when past end).
# A region "has an outbreak" in the window if ANY week is in_outbreak.
p_next_k <- function(k_yr) {
  k_wk <- round(k_yr * WPY)
  vapply(WEEKS, \(t) {
    idx <- ((t - 1L + seq_len(k_wk)) %% N_WK) + 1L
    mean(rowSums(ob_mat[, idx, drop = FALSE]) > 0L)
  }, numeric(1L))
}

K_VALS  <- c(1, 2, 5)
prob_df <- map_dfr(
  K_VALS,
  \(k) tibble(week = WEEKS, k_label = paste0("k = ", k, " yr"), prob = p_next_k(k))
)

# ── 5.  Per-region normalised incidence for heatmap ──────────────────────────
ts_df <- ts_df |>
  group_by(region) |>
  mutate(inc_norm = {
    li <- log1p(incidence)
    (li - min(li)) / (max(li) - min(li) + 1e-9)
  }) |>
  ungroup()

# ── 6.  Build plot ────────────────────────────────────────────────────────────
xbr  <- seq(1L, N_WK, by = WPY)
xlab <- paste0("Y", seq_len(N_YR))

k_cols <- c(
  "k = 1 yr" = "#4393c3",
  "k = 2 yr" = "#d6604d",
  "k = 5 yr" = "#762a83"
)

## 6a.  Top panel – probability lines
p_top <- ggplot(prob_df, aes(week, prob, colour = k_label)) +
  geom_line(linewidth = 0.85) +
  scale_x_continuous(breaks = xbr, labels = xlab, expand = c(0, 0)) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = label_percent(1),
    name   = expression(hat(P)(T[next] <= k ~ "yr"))
  ) +
  scale_colour_manual(values = k_cols, name = NULL) +
  labs(
    title    = "Estimated probability that a randomly chosen region has its next outbreak within k years",
    subtitle = "Outbreak: \u22654 consecutive weeks above mean + 2\u00b7SD  |  Circular window at series end"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.title.x     = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    plot.title       = element_text(size = 10, face = "bold"),
    plot.subtitle    = element_text(size = 8, colour = "grey50")
  )

## 6b.  Bottom panel – incidence heatmap
# r_ord = 1 → most outbreaks (top); r_ord = 50 → fewest (bottom)
p_bot <- ggplot(ts_df, aes(week, r_ord, fill = inc_norm)) +
  geom_raster() +
  scale_fill_viridis_c(
    option = "inferno",
    name   = "Incidence\n(norm. per region)",
    labels = c("low", "", "high"),
    breaks = c(0, 0.5, 1)
  ) +
  scale_x_continuous(breaks = xbr, labels = xlab, expand = c(0, 0)) +
  scale_y_reverse(
    breaks = c(1, 25, 50),
    labels = c("most\nfrequent", "25", "rarest"),
    name   = "Region (sorted by # outbreaks \u2193)"
  ) +
  labs(x = "Year") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid   = element_blank(),
    axis.ticks.y = element_line(colour = "grey70", linewidth = 0.3)
  )

## 6c.  Combine and save
fig <- (p_top / p_bot) + plot_layout(heights = c(1.5, 5))

out_path <- "outbreak_timeseries.png"
ggsave(out_path, fig, width = 14, height = 10, dpi = 150)
message("Saved: ", normalizePath(out_path))

# Save data objects for use by vaccination_simulation.R
saveRDS(ts_df, "data/ts_df_outbreaks.rds")
saveRDS(n_ep,  "data/n_ep_outbreaks.rds")
message("Data saved to data/ts_df_outbreaks.rds and data/n_ep_outbreaks.rds")
