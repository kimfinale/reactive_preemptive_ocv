# scripts/vaccination_sim_chrono.R
# ─────────────────────────────────────────────────────────────────────────────
# Chronological (round-based) vaccination simulation with carry-over immunity.
#
# Each "round" = one inter-epidemic interval (IEI_med ≈ 2.08 yr).
# m = floor(T_immune / IEI_med): number of future rounds protected by one dose.
#
# This directly operationalises the carry-over framework in simulation.qmd
# §Carry-Over Immunity and the Multi-Period Decision.
#
# Prerequisite: run scripts/fake_outbreak_timeseries.R first.
# Outputs: vaccination_sim_chrono.png
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
  library(MASS)
})
source("R/analytic_utils.R")

N_REG    <- 50L
N_YR     <- 10L
N_ROUNDS <- 5L     # epidemic rounds (≈ T_SIM / IEI_med; 5 × 2.08 yr ≈ 10.4 yr)

# ── 1. Load data ───────────────────────────────────────────────────────────────
if (!file.exists("data/ts_df_outbreaks.rds")) {
  stop("Run scripts/fake_outbreak_timeseries.R first.")
}
ts_df <- readRDS("data/ts_df_outbreaks.rds")
n_ep  <- readRDS("data/n_ep_outbreaks.rds")

# ── 2. IEI estimate ────────────────────────────────────────────────────────────
extract_iei <- function(region_df) {
  r_le   <- rle(region_df$in_outbreak)
  ends   <- cumsum(r_le$lengths)
  starts <- c(1L, ends[-length(ends)] + 1L)
  ep_idx <- which(r_le$values)
  if (length(ep_idx) < 2L) return(NULL)
  tibble(iei_yr = diff(starts[ep_idx]) / 52)
}

iei_df  <- ts_df |> arrange(region, week) |> group_split(region) |> map_dfr(extract_iei)
IEI_MED <- median(iei_df$iei_yr)

cat(sprintf("IEI median = %.2f yr  →  %d rounds ≈ %.1f yr\n\n",
  IEI_MED, N_ROUNDS, N_ROUNDS * IEI_MED))

# ── 3. Per-region outbreak probability: per-round Beta fit ─────────────────────
# Convert annual p_i to per-round p_i using Poisson consistency:
#   p_round = 1 - (1 - p_annual)^IEI_MED
region_p <- n_ep |>
  mutate(
    lambda_yr = n_episodes / N_YR,
    p_annual  = 1 - exp(-lambda_yr),
    p_round   = 1 - (1 - pmax(pmin(p_annual, 0.999), 0.001))^IEI_MED
  )

p_bar_round <- mean(region_p$p_round)
cat(sprintf("Annual p̄ = %.3f  →  per-round p̄ = %.3f\n", mean(region_p$p_annual), p_bar_round))

beta_fit <- tryCatch(
  MASS::fitdistr(region_p$p_round, "beta", start = list(shape1 = 1.5, shape2 = 6)),
  error = function(e) {
    message("Beta MLE failed; using MME.")
    list(estimate = c(shape1 = 1, shape2 = (1 - p_bar_round) / p_bar_round))
  }
)
P_SH <- beta_fit$estimate   # shape1, shape2 for per-round probabilities

cat(sprintf("Beta fit (per-round): a = %.2f, b = %.2f  (mean = %.3f)\n\n",
  P_SH[1], P_SH[2], P_SH[1] / sum(P_SH)))

# ── 4. Fixed parameters ────────────────────────────────────────────────────────
R_EFF <- 0.27   # endogenous reactive effectiveness (from vaccination_simulation.R)
NU    <- 1.0    # pre-emptive vaccine effectiveness

# ── 5. Analytical p_crit with carry-over ──────────────────────────────────────
# Derived in simulation.qmd §Carry-Over Immunity (smaller root of quadratic)
p_crit_carryover <- function(m, R, r, nu = 1) {
  if (m <= 0) return(1 / (1 + R * (nu - r)))
  a <- m * nu * R
  b <- -(1 + R * (nu - r) + m * nu * R)
  (-b - sqrt(b^2 - 4 * a)) / (2 * a)
}

# ── 6. Simulation functions ────────────────────────────────────────────────────

#' One simulation replicate: round-based chronological model.
#'
#' Each round k:
#'   (1) Pre-emptive doses allocated to top susceptible regions by targeting score.
#'   (2) Outbreaks (Bernoulli) processed; reactive doses allocated first-come-first-served.
#'
#' Same outbreak realisations reused across all alpha_grid values (variance reduction).
#'
#' @return tibble(alpha, cost) where cost = total_cost / (n * n_rounds)
sim_chrono_rep <- function(alpha_grid, f, R, nu, r_eff, m, n_rounds, n,
                            p_sh1, p_sh2, rho) {
  # Per-region per-round outbreak probability (fixed across rounds)
  p_i <- pmax(pmin(rbeta(n, p_sh1, p_sh2), 0.999), 0.001)

  # Targeting scores: Gaussian copula correlated with p_i at Spearman rho
  z_p     <- qnorm(rank(p_i, ties.method = "average") / (n + 1))
  z_sc    <- rho * z_p + sqrt(1 - rho^2) * rnorm(n)
  risk_ord <- order(z_sc, decreasing = TRUE)  # regions sorted by targeting score

  # Pre-generate outbreak events: n × n_rounds Bernoulli matrix
  ob_mat  <- matrix(runif(n * n_rounds) < rep(p_i, n_rounds), nrow = n)

  # Compute cost for each alpha (reuse ob_mat across all alphas)
  map_dfr(alpha_grid, function(alpha_val) {
    n_pre  <- round(alpha_val * f * n)
    n_reac <- round((1 - alpha_val) * f * n)

    immune_till <- rep(0L, n)   # immune through end of round immune_till[i]
    total_cost  <- 0

    for (k in seq_len(n_rounds)) {
      # ── Pre-emptive allocation ────────────────────────────────────────────
      # Select top susceptible (not yet immune) regions by targeting score
      if (n_pre > 0L) {
        suc_mask     <- immune_till < k               # logical: susceptible?
        suc_in_order <- risk_ord[suc_mask[risk_ord]]  # susceptible, sorted by risk
        n_vac        <- min(n_pre, length(suc_in_order))
        if (n_vac > 0L) {
          idx <- suc_in_order[seq_len(n_vac)]
          immune_till[idx] <- k + m   # immune through round k+m (carry-over)
          total_cost <- total_cost + n_vac
        }
      }

      # ── Outbreak events ───────────────────────────────────────────────────
      ob_idx   <- which(ob_mat[, k])
      reac_rem <- n_reac

      for (i in ob_idx) {
        if (immune_till[i] >= k) {
          # Immune: pay residual infection cost (fraction 1-nu)
          total_cost <- total_cost + R * (1 - nu)
        } else if (reac_rem > 0L) {
          # Susceptible, reactive capacity available
          reac_rem   <- reac_rem - 1L
          total_cost <- total_cost + 1 + R * (1 - r_eff)
          immune_till[i] <- k + m
        } else {
          # No vaccine: full infection cost
          total_cost <- total_cost + R
        }
      }
    }

    tibble(alpha = alpha_val, cost = total_cost / (n * n_rounds))
  })
}

#' Run N_SIM replicates for one (m, f, R, rho) scenario.
sim_chrono_batch <- function(m, f, R, rho, nu, r_eff, n_rounds, n,
                              alpha_grid, p_sh1, p_sh2, n_sim, seed) {
  set.seed(seed)
  map_dfr(seq_len(n_sim), \(s)
    sim_chrono_rep(alpha_grid, f, R, nu, r_eff, m, n_rounds, n,
                   p_sh1, p_sh2, rho)
  )
}

# ── 7. Run simulations ─────────────────────────────────────────────────────────
T_IMM_VALS <- c(2, 4, 6)     # immunity-window scenarios (years) → m = 0, 1, 2
F_VALS     <- c(0.3, 0.5)
R_VALS     <- c(2, 5, 10)
RHO_VALS   <- c(1.0, 0.7)
ALPHA_GRID <- seq(0, 1, by = 0.1)
N_SIM      <- 300L

cat("Running chronological simulations ...\n")
sim_res <- map_dfr(T_IMM_VALS, function(T_imm) {
  m_val <- floor(T_imm / IEI_MED)
  cat(sprintf("  T_immune = %d yr  (m = %d rounds carry-over)\n", T_imm, m_val))
  map_dfr(F_VALS, function(f_val) {
    map_dfr(R_VALS, function(R_val) {
      map_dfr(RHO_VALS, function(rho_val) {
        cat(sprintf("    f=%.1f  R=%2.0f  rho=%.1f ...\n", f_val, R_val, rho_val))
        sim_chrono_batch(
          m        = m_val,
          f        = f_val,
          R        = R_val,
          rho      = rho_val,
          nu       = NU,
          r_eff    = R_EFF,
          n_rounds = N_ROUNDS,
          n        = N_REG,
          alpha_grid = ALPHA_GRID,
          p_sh1    = P_SH[1],
          p_sh2    = P_SH[2],
          n_sim    = N_SIM,
          seed     = 42L
        ) |> mutate(
          T_immune = T_imm,
          m        = m_val,
          f        = f_val,
          R        = R_val,
          rho      = rho_val
        )
      })
    })
  })
})

# ── 8. Summarise ───────────────────────────────────────────────────────────────
cost_smry <- sim_res |>
  group_by(T_immune, m, f, R, rho, alpha) |>
  summarise(
    cost_mean = mean(cost),
    cost_lo   = quantile(cost, 0.10),
    cost_hi   = quantile(cost, 0.90),
    .groups   = "drop"
  ) |>
  mutate(
    T_lbl = paste0(T_immune, " yr  (m=", m, ")"),
    f_lbl = paste0("f = ", f),
    R_lbl = paste0("R = ", R)
  )

policy_tbl <- cost_smry |>
  group_by(T_immune, m, f, R, rho, T_lbl, f_lbl, R_lbl) |>
  summarise(
    alpha_star      = alpha[which.min(cost_mean)],
    cost_star       = min(cost_mean),
    cost_reactive   = cost_mean[alpha == 0],
    cost_preemptive = cost_mean[alpha == 1],
    save_vs_react   = (cost_reactive - cost_star) / cost_reactive * 100,
    .groups = "drop"
  )

cat("\n─── Policy summary (rho=1, f=0.5) ──────────────────────────────────────\n")
policy_tbl |>
  filter(rho == 1, f == 0.5) |>
  dplyr::select(T_immune, m, R, alpha_star, cost_reactive, cost_star, save_vs_react) |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  print(n = Inf)

# ── 9. Analytical comparison table ────────────────────────────────────────────
analytic_df <- expand_grid(T_immune = T_IMM_VALS, R = R_VALS) |>
  mutate(
    m         = floor(T_immune / IEI_MED),
    p_crit_m  = map2_dbl(m, R, \(mv, Rv) p_crit_carryover(mv, Rv, r = R_EFF)),
    pct_above = map_dbl(p_crit_m, \(pc) mean(region_p$p_round > pc) * 100)
  )

cat("\n─── Analytical p_crit^(m) vs simulation α* (rho=1, f=0.5) ──────────────\n")
policy_tbl |>
  filter(rho == 1, f == 0.5) |>
  left_join(analytic_df |> dplyr::select(T_immune, R, m, p_crit_m, pct_above),
            by = c("T_immune", "R", "m")) |>
  dplyr::select(T_immune, m, R, p_crit_m, pct_above, alpha_star, save_vs_react) |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  print(n = Inf)

# ── 10. Plots ──────────────────────────────────────────────────────────────────
T_COLS <- c(
  "2 yr  (m=0)" = "#4393c3",
  "4 yr  (m=1)" = "#d6604d",
  "6 yr  (m=2)" = "#762a83"
)

## Panel A: Cost vs α — R=5, f=0.5, rho=1
pA_dat     <- cost_smry |> filter(R == 5, f == 0.5, rho == 1)
pA_opt     <- policy_tbl |> filter(R == 5, f == 0.5, rho == 1)

pA <- ggplot(pA_dat, aes(alpha, cost_mean, colour = T_lbl, fill = T_lbl)) +
  geom_ribbon(aes(ymin = cost_lo, ymax = cost_hi), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1.0) +
  geom_point(
    data = pA_opt,
    aes(x = alpha_star, y = cost_star, colour = T_lbl),
    shape = 18, size = 4.5, show.legend = FALSE
  ) +
  scale_x_continuous(labels = label_percent(1),
                     name = "Pre-emptive fraction (\u03b1)") +
  scale_y_continuous(name = "Expected cost (C\u1d5b \u00b7 region\u207b\u00b9 \u00b7 round\u207b\u00b9)") +
  scale_colour_manual(values = T_COLS, name = expression(T[immune])) +
  scale_fill_manual(  values = T_COLS, name = expression(T[immune])) +
  labs(
    title    = "A. Carry-over immunity shifts optimal pre-emptive allocation",
    subtitle = "R = 5, f = 0.5, \u03c1 = 1; diamonds = \u03b1*; shading = 10\u201390th pctile"
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey50"),
        legend.position = "right")

## Panel B: Cost vs α — by R, f=0.5, rho=1 (small multiples)
pB_dat <- cost_smry |> filter(f == 0.5, rho == 1)
pB_opt <- policy_tbl |> filter(f == 0.5, rho == 1)

pB <- ggplot(pB_dat, aes(alpha, cost_mean, colour = T_lbl)) +
  geom_line(linewidth = 0.8) +
  geom_point(
    data = pB_opt,
    aes(x = alpha_star, y = cost_star, colour = T_lbl),
    shape = 18, size = 3.5, show.legend = FALSE
  ) +
  facet_wrap(~ R_lbl, scales = "free_y", nrow = 1) +
  scale_x_continuous(labels = label_percent(1),
                     name = "Pre-emptive fraction (\u03b1)") +
  scale_y_continuous(name = "Expected cost") +
  scale_colour_manual(values = T_COLS, name = expression(T[immune])) +
  labs(
    title    = "B. Cost vs. \u03b1 across cost ratios  (f = 0.5, \u03c1 = 1)",
    subtitle = "Diamonds = \u03b1*"
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey50"),
        strip.text    = element_text(face = "bold"),
        legend.position = "right")

## Panel C: Optimal α* vs R — by T_immune, linetype by f, rho=1
f_lty <- c("f = 0.3" = "dashed", "f = 0.5" = "solid")

pC <- policy_tbl |>
  filter(rho == 1) |>
  ggplot(aes(R, alpha_star, colour = T_lbl, linetype = f_lbl)) +
  geom_hline(yintercept = p_bar_round, linetype = "dotted", colour = "grey40") +
  annotate("text", x = 10.2, y = p_bar_round + 0.01, hjust = 0, size = 2.8,
           colour = "grey40",
           label = sprintf("p\u0305\u1d63\u2099\u1d48 = %.2f", p_bar_round)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = R_VALS,
                     name = "Cost ratio (R = C\u1d35/C\u1d5b)") +
  scale_y_continuous(limits = c(0, 1), labels = label_percent(1),
                     name = "Optimal pre-emptive fraction (\u03b1*)") +
  scale_colour_manual(values = T_COLS, name = expression(T[immune])) +
  scale_linetype_manual(values = f_lty, name = "Capacity") +
  labs(
    title    = "C. Optimal \u03b1* rises with R and with T_immune",
    subtitle = "\u03c1 = 1 (perfect targeting)"
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey50"),
        legend.position = "right")

## Panel D: % cost saving vs T_immune — by R, f=0.5, rho=1
R_COLS <- c("R = 2" = "#4393c3", "R = 5" = "#d6604d", "R = 10" = "#762a83")

pD <- policy_tbl |>
  filter(rho == 1, f == 0.5) |>
  mutate(m_lbl = paste0("m=", m)) |>
  ggplot(aes(as.factor(T_immune), save_vs_react, colour = R_lbl, group = R_lbl)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_x_discrete(
    name   = expression(T[immune] ~ "(years)"),
    labels = paste0(T_IMM_VALS, " yr\n(m=", floor(T_IMM_VALS / IEI_MED), ")")
  ) +
  scale_y_continuous(labels = label_percent(1), limits = c(0, NA),
                     name = "Cost saving vs. pure reactive (%)") +
  scale_colour_manual(values = R_COLS, name = "R") +
  labs(
    title    = "D. Value of carry-over immunity: savings over pure reactive",
    subtitle = "f = 0.5, \u03c1 = 1"
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey50"),
        legend.position = "right")

## Assemble
fig <- (pA | pC) / (pB | pD) +
  plot_annotation(
    title = "Chronological simulation: carry-over immunity and pre-emptive vaccination policy",
    subtitle = sprintf(
      "Round-based model | IEI = %.1f yr | %d rounds (\u2248%.0f yr) | n = %d regions | N\u2090\u1d35\u2098 = %d",
      IEI_MED, N_ROUNDS, N_ROUNDS * IEI_MED, N_REG, N_SIM
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  )

out <- "vaccination_sim_chrono.png"
ggsave(out, fig, width = 17, height = 11, dpi = 150)
message("Saved: ", normalizePath(out))
