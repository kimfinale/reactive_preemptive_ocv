# scripts/vaccination_simulation.R
# ─────────────────────────────────────────────────────────────────────────────
# Vaccination allocation simulations for 50 regions, grounded in the simulated
# outbreak time series.
#
# Prerequisite: run scripts/fake_outbreak_timeseries.R first to generate
#   data/ts_df_outbreaks.rds  and  data/n_ep_outbreaks.rds
#
# Outputs: vaccination_simulation.png
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
  library(MASS)
})
source("R/analytic_utils.R")

N_REG <- 50L
N_YR  <- 10L

# ── 1. Load time series data ──────────────────────────────────────────────────
if (!file.exists("data/ts_df_outbreaks.rds")) {
  stop("Run scripts/fake_outbreak_timeseries.R first.")
}
ts_df <- readRDS("data/ts_df_outbreaks.rds")
n_ep  <- readRDS("data/n_ep_outbreaks.rds")

# ── 2. Extract outbreak episode characteristics ───────────────────────────────
# For each qualifying episode: duration (weeks), peak incidence, time to peak

extract_episodes <- function(region_df) {
  r_le   <- rle(region_df$in_outbreak)
  ends   <- cumsum(r_le$lengths)
  starts <- c(1L, ends[-length(ends)] + 1L)
  ep_idx <- which(r_le$values)
  if (length(ep_idx) == 0L) return(NULL)

  map_dfr(ep_idx, function(i) {
    wks   <- starts[i]:ends[i]
    edata <- region_df[wks, ]
    pk    <- which.max(edata$incidence)
    tibble(
      duration_wk     = length(wks),
      peak_inc        = edata$incidence[pk],
      time_to_peak_wk = pk,
      peak_frac       = pk / length(wks)
    )
  })
}

episode_df <- ts_df |>
  arrange(region, week) |>
  group_split(region) |>
  map_dfr(extract_episodes)

cat("─── Episode characteristics ─────────────────────────────────────────\n")
cat(sprintf("  n episodes total : %d\n", nrow(episode_df)))
cat(sprintf("  Duration (wk)    : median %.1f, IQR [%.0f, %.0f]\n",
  median(episode_df$duration_wk),
  quantile(episode_df$duration_wk, 0.25),
  quantile(episode_df$duration_wk, 0.75)))
cat(sprintf("  Peak timing      : mean frac = %.2f (%.0f%% of way through)\n",
  mean(episode_df$peak_frac), 100 * mean(episode_df$peak_frac)))
cat(sprintf("  Peak incidence   : median %.0f cases/wk\n\n",
  median(episode_df$peak_inc)))

# ── 2b. Inter-epidemic intervals (IEI) ───────────────────────────────────────
# IEI = time between starts of consecutive episodes within the same region

extract_iei <- function(region_df) {
  r_le      <- rle(region_df$in_outbreak)
  ends      <- cumsum(r_le$lengths)
  starts    <- c(1L, ends[-length(ends)] + 1L)
  ep_idx    <- which(r_le$values)
  if (length(ep_idx) < 2L) return(NULL)
  ep_starts <- starts[ep_idx]
  tibble(iei_yr = diff(ep_starts) / 52)
}

iei_df  <- ts_df |>
  arrange(region, week) |>
  group_split(region) |>
  map_dfr(extract_iei)

IEI_med <- median(iei_df$iei_yr, na.rm = TRUE)
IEI_mn  <- mean(iei_df$iei_yr,   na.rm = TRUE)

cat(sprintf("─── Inter-epidemic interval (IEI) ───────────────────────────────────\n"))
cat(sprintf("  IEI  median = %.2f yr  |  mean = %.2f yr\n", IEI_med, IEI_mn))
cat(sprintf("  IEI  10th–90th pct: [%.2f, %.2f] yr\n\n",
  quantile(iei_df$iei_yr, 0.10), quantile(iei_df$iei_yr, 0.90)))

# ── 3. Per-region annual outbreak probability ─────────────────────────────────
# Annual Poisson rate lambda_i = n_episodes / 10 years
# Annual probability p_i = 1 - exp(-lambda_i)
region_p <- n_ep |>
  mutate(
    lambda_yr = n_episodes / N_YR,
    p_annual  = 1 - exp(-lambda_yr)
  )

p_bar <- mean(region_p$p_annual)
cat(sprintf("─── Per-region outbreak probability ─────────────────────────────────\n"))
cat(sprintf("  Mean p̄ = %.3f   SD = %.3f\n", p_bar, sd(region_p$p_annual)))

# ── 4. Fit parametric distribution to p_i ────────────────────────────────────
# (a) Constrained Beta(1, θ) consistent with existing framework
theta_mme <- (1 - p_bar) / p_bar          # method-of-moments

# (b) Unconstrained Beta(a, b) via MLE (for comparison)
p_fit_vals <- pmax(pmin(region_p$p_annual, 0.999), 0.001)
beta_mle <- tryCatch(
  MASS::fitdistr(p_fit_vals, "beta",
                 start = list(shape1 = 1.5, shape2 = 6)),
  error = function(e) NULL
)

cat(sprintf("  Beta(1, θ) fit   : θ = %.2f (mean = %.3f)\n",
  theta_mme, 1 / (1 + theta_mme)))
if (!is.null(beta_mle)) {
  a_mle <- beta_mle$estimate["shape1"]
  b_mle <- beta_mle$estimate["shape2"]
  cat(sprintf("  Beta(a, b) MLE   : a = %.2f, b = %.2f (mean = %.3f)\n\n",
    a_mle, b_mle, a_mle / (a_mle + b_mle)))
}

# ── 5. Calibrate endogenous simulation parameters from episode stats ──────────
# Duration: 10th–90th percentile range in days (×7 to convert weeks → days)
dur_range_days <- as.numeric(quantile(episode_df$duration_wk, c(0.10, 0.90))) * 7L

# Peak fraction: fit Beta(a, b) to observed time_to_peak / duration
pk_fit <- tryCatch(
  MASS::fitdistr(
    pmax(pmin(episode_df$peak_frac, 0.99), 0.01),
    "beta",
    start = list(shape1 = 2, shape2 = 4)
  ),
  error = function(e) {
    message("Peak frac Beta fit failed; using Beta(2, 4) defaults.")
    list(estimate = c(shape1 = 2, shape2 = 4))
  }
)
pk_shape <- as.numeric(pk_fit$estimate)   # c(shape1, shape2)

# Response delay: 2–6 weeks (reactive vaccination campaigns)
delay_range_days <- c(14, 42)

# Growth/decay rates: not directly extractable from weekly data; use defaults
rg_range <- c(0.03, 0.10)   # /day
rd_range <- c(0.01, 0.05)   # /day

cat(sprintf("─── Calibrated simulation parameters ────────────────────────────────\n"))
cat(sprintf("  Outbreak duration : [%.0f, %.0f] days\n",
  dur_range_days[1], dur_range_days[2]))
cat(sprintf("  Peak frac Beta    : (%.2f, %.2f)  [mean = %.2f]\n",
  pk_shape[1], pk_shape[2], pk_shape[1] / sum(pk_shape)))
cat(sprintf("  Response delay    : [%d, %d] days\n\n",
  delay_range_days[1], delay_range_days[2]))

# ── 6. Vaccination simulations ────────────────────────────────────────────────
# For each (f, R): compute expected cost across alpha grid, find optimal alpha*
# Uses endogenous r: r_i = OVE × frac_remaining(t_delay)
# VE parameters: VEd=0.7, VEi=0.3, coverage f_cov=0.7 (cholera benchmarks)
# Targeting: rho=1 (perfect), rho=0.7 (realistic), rho=0 (random/no targeting)

F_VALS     <- c(0.2, 0.3, 0.5)
R_VALS     <- c(2, 5, 10)
RHO_VALS   <- c(1.0, 0.7, 0.0)
ALPHA_GRID <- seq(0, 1, by = 0.1)
N_SIM      <- 150L

cat("Running simulations ...\n")
sim_res <- map_dfr(RHO_VALS, function(rho_val) {
  map_dfr(F_VALS, function(f_val) {
    map_dfr(R_VALS, function(R_val) {
      cat(sprintf("  rho=%.1f  f=%.1f  R=%2.0f\n", rho_val, f_val, R_val))
      sim <- sim_cost_heterorisk_endogenous(
        alpha           = ALPHA_GRID,
        n               = N_REG,
        p               = p_bar,
        f               = f_val,
        R               = R_val,
        VEd             = 0.7,
        VEi             = 0.3,
        f_cov           = 0.7,
        nu              = 1,
        rho             = rho_val,
        duration_range  = dur_range_days,
        peak_frac_shape = pk_shape,
        rg_range        = rg_range,
        rd_range        = rd_range,
        delay_range     = delay_range_days,
        n_sim           = N_SIM,
        seed            = 42
      )
      sim |> mutate(
        f     = f_val,  R     = R_val,  rho   = rho_val,
        f_lbl = paste0("f = ", f_val),
        R_lbl = paste0("R = ", R_val),
        rho_lbl = dplyr::case_when(
          rho_val == 1.0 ~ "\u03c1 = 1 (perfect)",
          rho_val == 0.7 ~ "\u03c1 = 0.7 (realistic)",
          rho_val == 0.0 ~ "\u03c1 = 0 (random)"
        )
      )
    })
  })
})

# ── 7. Summarise ──────────────────────────────────────────────────────────────
cost_smry <- sim_res |>
  group_by(f, R, rho, f_lbl, R_lbl, rho_lbl, alpha) |>
  summarise(
    cost_mean = mean(cost),
    cost_lo   = quantile(cost, 0.10),
    cost_hi   = quantile(cost, 0.90),
    mean_r    = mean(mean_r, na.rm = TRUE),
    .groups   = "drop"
  )

policy_tbl <- cost_smry |>
  group_by(f, R, rho, f_lbl, R_lbl, rho_lbl) |>
  summarise(
    alpha_star      = alpha[which.min(cost_mean)],
    cost_star       = min(cost_mean),
    cost_reactive   = cost_mean[alpha == 0],
    cost_preemptive = cost_mean[alpha == 1],
    save_vs_react   = (cost_reactive - cost_star) / cost_reactive * 100,
    save_vs_pre     = (cost_preemptive - cost_star) / cost_preemptive * 100,
    eff_r           = mean_r[alpha == min(alpha[!is.na(mean_r)])],
    .groups         = "drop"
  )

cat("\n─── Policy summary (rho=1, perfect targeting) ───────────────────────\n")
policy_tbl |>
  filter(rho == 1) |>
  dplyr::select(f, R, alpha_star, save_vs_react, save_vs_pre, eff_r) |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  print()

# ── 8. Effective r from simulation ────────────────────────────────────────────
eff_r_val <- cost_smry |>
  filter(rho == 1, alpha > 0) |>
  pull(mean_r) |>
  mean(na.rm = TRUE)
cat(sprintf("\n  Effective reactive r from trajectories: %.3f\n\n", eff_r_val))

# ── 9. p_crit thresholds (analytical) ────────────────────────────────────────
p_crit_tbl <- expand_grid(f = F_VALS, R = R_VALS) |>
  mutate(
    r_eff    = eff_r_val,
    p_crit   = map2_dbl(R, f, \(Rv, fv)
      p_star_multi(R = Rv, r = eff_r_val, f = fv, nu = 1)),
    pct_regions_above = map_dbl(p_crit, \(pc)
      mean(region_p$p_annual > pc, na.rm = TRUE) * 100)
  )

cat("─── p_crit thresholds vs observed regional probabilities ────────────\n")
print(p_crit_tbl |> mutate(across(where(is.numeric), \(x) round(x, 3))))

# ── 9b. Carry-over p_crit analysis ───────────────────────────────────────────
# From simulation.qmd Section: Carry-Over Immunity and the Multi-Period Decision
# Modified p_crit satisfies: m*nu*R*p^2 - [1+R(nu-r)+m*nu*R]*p + 1 = 0
# Smaller root (physically meaningful) is:
#   p_crit^(m) = {[1+R(nu-r)+mnu*R] - sqrt(discriminant)} / (2*m*nu*R)

p_crit_carryover <- function(m, R, r, nu = 1) {
  if (m <= 0) return(1 / (1 + R * (nu - r)))
  a <- m * nu * R
  b <- -(1 + R * (nu - r) + m * nu * R)
  (-b - sqrt(b^2 - 4 * a)) / (2 * a)
}

# Carry-over multiplier (undiscounted): m = floor(T_immune / IEI)
T_IMM_VALS <- c(2, 4, 6)   # three immunity-window scenarios (years)

carryover_tbl <- expand_grid(T_immune = T_IMM_VALS, R = R_VALS) |>
  mutate(
    m           = floor(T_immune / IEI_med),
    p_crit0     = map_dbl(R, \(Rv)
                    p_crit_carryover(0, Rv, r = eff_r_val)),
    p_crit_m    = map2_dbl(m, R, \(mv, Rv)
                    p_crit_carryover(mv, Rv, r = eff_r_val)),
    pct_above_0 = map_dbl(p_crit0,  \(pc)
                    mean(region_p$p_annual > pc) * 100),
    pct_above_m = map_dbl(p_crit_m, \(pc)
                    mean(region_p$p_annual > pc) * 100)
  )

cat(sprintf("\n─── Carry-over p_crit  (IEI median = %.2f yr; r = %.3f) ─────────────\n",
  IEI_med, eff_r_val))
carryover_tbl |>
  dplyr::select(T_immune, R, m, p_crit0, p_crit_m, pct_above_0, pct_above_m) |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  print(n = Inf)
cat("\n")

# ── 10. Plots ─────────────────────────────────────────────────────────────────
f_cols  <- c("f = 0.2" = "#d6604d", "f = 0.3" = "#4393c3", "f = 0.5" = "#4d9221")
rho_lty <- c("\u03c1 = 1 (perfect)"    = "solid",
             "\u03c1 = 0.7 (realistic)" = "dashed",
             "\u03c1 = 0 (random)"      = "dotted")

## Panel A: Per-region p distribution + Beta fits
x_seq <- seq(0, 0.55, length.out = 300)
fit_df <- bind_rows(
  tibble(p = x_seq, d = dbeta(x_seq, 1, theta_mme), fit = "Beta(1, \u03b8) MME"),
  if (!is.null(beta_mle)) {
    tibble(p = x_seq,
           d = dbeta(x_seq, beta_mle$estimate["shape1"], beta_mle$estimate["shape2"]),
           fit = sprintf("Beta(%.1f, %.1f) MLE", beta_mle$estimate["shape1"],
                         beta_mle$estimate["shape2"]))
  }
)
p_A <- ggplot(region_p, aes(x = p_annual)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.05,
                 fill = "#4393c3", alpha = 0.65, colour = "white") +
  geom_line(data = fit_df, aes(p, d, colour = fit, linetype = fit), linewidth = 0.9) +
  geom_vline(xintercept = p_bar, linetype = "dotted", colour = "grey40") +
  annotate("text", x = p_bar + 0.02, y = 4.5,
           label = sprintf("p\u0305 = %.2f", p_bar), size = 3, hjust = 0) +
  scale_x_continuous(labels = label_percent(1), expand = c(0, 0),
                     name = "Annual outbreak probability (p\u1d62)") +
  scale_colour_manual(values = c("#b2182b", "#762a83"), name = NULL) +
  scale_linetype_manual(values = c("dashed", "solid"), name = NULL) +
  labs(y = "Density",
       title = "A. Per-region outbreak risk distribution") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = c(0.65, 0.75),
        legend.background = element_rect(fill = "white", colour = NA))

## Panel B: Episode duration + peak timing
p_B1 <- ggplot(episode_df, aes(x = duration_wk)) +
  geom_histogram(binwidth = 2, fill = "#4393c3", alpha = 0.7, colour = "white") +
  scale_x_continuous(name = "Duration (weeks)") +
  labs(y = "Count", title = "B. Outbreak episode characteristics",
       subtitle = "Duration distribution") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey50"))

p_B2 <- ggplot(episode_df, aes(x = peak_frac)) +
  geom_histogram(binwidth = 0.1, fill = "#d6604d", alpha = 0.7, colour = "white") +
  stat_function(fun = dbeta, args = list(shape1 = pk_shape[1], shape2 = pk_shape[2]),
                colour = "#762a83", linewidth = 1) +
  scale_x_continuous(labels = label_percent(1), name = "Peak timing (fraction of episode)") +
  labs(y = "Count", subtitle = "Peak timing (curve = fitted Beta)") +
  theme_minimal(base_size = 10) +
  theme(plot.subtitle = element_text(size = 8, colour = "grey50"))

p_B <- p_B1 / p_B2 + plot_layout(heights = c(1, 1))

## Panel C: Cost vs alpha (R=5, both rho values)
p_C <- cost_smry |>
  filter(R == 5) |>
  ggplot(aes(alpha, cost_mean,
             colour = f_lbl, linetype = rho_lbl, group = interaction(f_lbl, rho_lbl))) +
  geom_ribbon(aes(ymin = cost_lo, ymax = cost_hi, fill = f_lbl),
              alpha = 0.10, colour = NA,
              data = \(d) filter(d, rho == 1)) +
  geom_line(linewidth = 0.85) +
  geom_point(
    data = cost_smry |> filter(R == 5) |>
      group_by(f_lbl, rho_lbl) |> slice_min(cost_mean, n = 1),
    shape = 18, size = 4, show.legend = FALSE
  ) +
  scale_x_continuous(labels = label_percent(1),
                     name = "Pre-emptive fraction (\u03b1)") +
  scale_colour_manual(values = f_cols, name = "Capacity") +
  scale_fill_manual(values = f_cols, name = "Capacity") +
  scale_linetype_manual(values = rho_lty, name = "Targeting") +
  labs(y = "Expected cost / C\u1d5b",
       title = "C. Cost vs. allocation strategy (R = 5)",
       subtitle = "Diamonds = \u03b1*; shading = 10\u201390th percentile (perfect targeting)") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey50"),
        legend.position = "right")

## Panel D: Optimal alpha* by R and f (perfect targeting only)
p_D <- policy_tbl |>
  filter(rho == 1) |>
  ggplot(aes(R, alpha_star, colour = f_lbl, group = f_lbl)) +
  geom_hline(yintercept = p_bar, linetype = "dotted", colour = "grey60") +
  annotate("text", x = 10.2, y = p_bar + 0.02,
           label = sprintf("p\u0305 = %.2f", p_bar), size = 3, hjust = 0) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = R_VALS,
                     name = "Cost ratio (R = C\u1d35 / C\u1d5b)") +
  scale_y_continuous(limits = c(0, 1), labels = label_percent(1),
                     name = "Optimal pre-emptive fraction (\u03b1*)") +
  scale_colour_manual(values = f_cols, name = "Capacity") +
  labs(title = "D. Optimal \u03b1* by capacity and cost ratio",
       subtitle = sprintf("p\u0305 = %.2f; effective r \u2248 %.2f; perfect targeting",
                          p_bar, eff_r_val)) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey50"),
        legend.position = "right")

## Panel E: Carry-over p_crit vs T_immune (step function due to floor())
T_range   <- seq(0, 8, by = 0.05)
R_cols_E  <- c("R = 2" = "#4393c3", "R = 5" = "#d6604d", "R = 10" = "#762a83")

panel_e_df <- expand_grid(T_immune = T_range, R = R_VALS) |>
  mutate(
    m        = floor(T_immune / IEI_med),
    p_crit_m = map2_dbl(m, R, \(mv, Rv)
                 p_crit_carryover(mv, Rv, r = eff_r_val)),
    R_lbl    = paste0("R = ", R)
  )

p_E <- ggplot(panel_e_df, aes(T_immune, p_crit_m, colour = R_lbl)) +
  geom_step(linewidth = 0.85) +
  geom_hline(yintercept = p_bar, linetype = "dotted", colour = "grey40") +
  geom_vline(xintercept = IEI_med, linetype = "dashed", colour = "grey60",
             linewidth = 0.6) +
  annotate("text", x = IEI_med + 0.12, y = 0.94, hjust = 0, size = 2.8,
           colour = "grey40",
           label = sprintf("IEI = %.1f yr", IEI_med)) +
  annotate("text", x = 7.8, y = p_bar + 0.03, hjust = 1, size = 2.8,
           colour = "grey40",
           label = sprintf("p\u0305 = %.2f", p_bar)) +
  scale_x_continuous(name = expression(T[immune] ~ "(years)"),
                     limits = c(0, 8)) +
  scale_y_continuous(limits = c(0, 1), labels = label_percent(1),
                     name = expression(italic(p)[crit]^{(m)})) +
  scale_colour_manual(values = R_cols_E, name = NULL) +
  labs(
    title    = "E. Carry-over immunity lowers critical outbreak probability",
    subtitle = sprintf(
      "m = floor(T_immune / IEI); IEI (median) = %.1f yr; r \u2248 %.2f",
      IEI_med, eff_r_val)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(size = 8, colour = "grey50"),
    legend.position = "right"
  )

## Assemble
fig <- (p_A | p_B) / (p_C | p_D | p_E) +
  plot_annotation(
    title = "Vaccination allocation policy for 50 regions",
    subtitle = paste0(
      "Simulated outbreak time series (10 yr) \u2192 empirical p\u1d62 \u2192 ",
      "Beta(1, \u03b8)-calibrated endogenous simulation (n = 150 reps)"
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  )

out <- "vaccination_simulation.png"
ggsave(out, fig, width = 20, height = 11, dpi = 150)
message("Saved: ", normalizePath(out))

# ── 11. Concise policy table ──────────────────────────────────────────────────
cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat("POLICY RECOMMENDATION TABLE\n")
cat("Mean outbreak probability p̄ =", round(p_bar, 3), "\n")
cat("Effective reactive effectiveness r ≈", round(eff_r_val, 2), "\n")
cat("───────────────────────────────────────────────────────────────────────\n")
policy_tbl |>
  filter(rho == 1) |>
  transmute(
    Capacity   = f_lbl,
    `R=C_I/C_V`  = R,
    `α* (pre-emp %)` = paste0(round(alpha_star * 100), "%"),
    `Save vs reactive` = paste0(round(save_vs_react, 1), "%"),
    `Save vs pre-emp`  = paste0(round(save_vs_pre,  1), "%")
  ) |>
  print(n = Inf)
cat("═══════════════════════════════════════════════════════════════════════\n")
