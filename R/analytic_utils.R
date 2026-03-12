#' Pre-emptive cost (normalized by C_V)
#' @export
cost_pre_one <- function(p, R, r, nu = 1) {
    # p and r unused in the original (nu=1) case but kept for interface.
    # With nu < 1: Cost = 1 + p * (1 - nu) * R
    1 + p * (1 - nu) * R
}

#' Reactive cost (normalized by C_V)
#' @export
cost_react_one <- function(p, R, r, nu = 1) {
    p * (1 + (1 - r) * R)
}

#' Threshold probability for single population
#' @export
p_star_one <- function(R, r, nu = 1) {
    # Solve: 1 + p(1-nu)R = p(1 + (1-r)R)
    # 1 + pR - p*nu*R = p + pR - prR
    # 1 = p (1 + R - rR - R + nu*R) = p(1 + R(nu - r))
    1 / (1 + R * (nu - r))
}

#' Multi-population pre-emptive cost
#' @export
cost_pre_multi <- function(p, R, r, f, nu = 1) {
    # Pre-emptive group (f): Cost = f * (1 + p(1-nu)R)
    # Unvaccinated group (1-f): Cost = (1-f) * pR
    # Total = f + fpR - fp*nu*R + pR - fpR
    #       = f + pR - fp*nu*R
    #       = f + pR(1 - f*nu)
    f + p * R * (1 - f * nu)
}

#' Multi-population reactive cost
#' @export
cost_react_multi <- function(p, R, r, f, nu = 1) {
    ifelse(
        f < p,
        # reactive-limited
        f + (p - f * r) * R,
        # reactive-rich
        p * (1 + (1 - r) * R)
    )
}

#' Cost difference
#' @export
cost_diff_multi <- function(p, R, r, f, nu = 1) {
    c_pre <- cost_pre_multi(p, R, r, f, nu)
    c_react <- cost_react_multi(p, R, r, f)
    c_pre - c_react
}

#' Threshold probability for multi-population
#' @export
p_star_multi <- function(R, r, f, nu = 1, tol = 1e-10) {
    # Candidate 1: reactive-limited regime (p > f)
    # c_pre = c_react_lim
    # f + pR(1 - f*nu) = f + (p - f*r)R
    # pR - pf*nu*R = pR - frR
    # -pf*nu*R = -frR
    # p * nu = r  => p* = r / nu

    if (nu == 0) {
        p_star_RL <- NA_real_ # Avoid division by zero
    } else {
        p_star_RL <- r / nu
    }

    valid_RL <- !is.na(p_star_RL) && (p_star_RL > f) && (p_star_RL >= 0) && (p_star_RL <= 1)
    if (valid_RL) {
        diff_val <- cost_diff_multi(p_star_RL, R = R, r = r, f = f, nu = nu)
        if (abs(diff_val) < 1e-6) {
            return(p_star_RL)
        }
    }

    # Candidate 2: reactive-rich regime (p <= f)
    # c_pre = c_react_rich
    # f + pR(1 - f*nu) = p(1 + (1-r)R)
    # f + pR - pf*nu*R = p + pR - prR
    # f = p (1 + R - rR - R + f*nu*R)
    # f = p (1 + R(f*nu - r))
    # p* = f / (1 + R(f*nu - r))

    denom <- 1 + R * (f * nu - r)
    if (abs(denom) < tol) {
        p_star_RR <- NA_real_
    } else {
        p_star_RR <- f / denom
    }

    valid_RR <- !is.na(p_star_RR) &&
        (p_star_RR >= 0) && (p_star_RR <= f) && (p_star_RR <= 1)

    if (valid_RR) {
        diff_val <- cost_diff_multi(p_star_RR, R = R, r = r, f = f, nu = nu)
        if (abs(diff_val) < 1e-6) {
            return(p_star_RR)
        }
    }

    # Fallback: numeric search
    p_grid <- seq(1e-6, 1 - 1e-6, length.out = 2001)
    vals <- cost_diff_multi(p_grid, R = R, r = r, f = f, nu = nu)

    # Check for sign change
    idx <- which(vals[-1] * vals[-length(vals)] <= 0)
    if (length(idx) == 0) {
        return(NA_real_)
    }

    i <- idx[1]
    lower <- p_grid[i]
    upper <- p_grid[i + 1]

    uniroot(cost_diff_multi,
        lower = lower, upper = upper,
        R = R, r = r, f = f, nu = nu, tol = tol
    )$root
}

#' Simulate Lambda and convert to P
#' @export
sim_lambda_P <- function(M, theta, seed = 1) {
    set.seed(seed)
    lambda <- rexp(M, rate = theta) # Lambda ~ Exp(theta)
    P <- 1 - exp(-lambda) # P = 1 - exp(-Lambda)
    list(lambda = lambda, P = P)
}

#' Given given sigma, compute Spearman rank correlation between Lambda and S
#' @export
spearman_rho_given_sigma <- function(lambda, sigma, seed = 1) {
    set.seed(seed)
    S <- lambda + rnorm(length(lambda), mean = 0, sd = sigma)
    suppressWarnings(cor(lambda, S, method = "spearman"))
}

#' Calibrate sigma so that Spearman cor(Lambda, S) ~= rho_target
#' @export
calibrate_sigma_for_rho <- function(
    lambda, rho_target,
    sigma_hi = 50, seed = 1) {
    # monotone: rho decreases as sigma increases
    f_obj <- function(log_sigma) {
        sigma <- exp(log_sigma)
        rho_hat <- spearman_rho_given_sigma(lambda, sigma, seed = seed)
        rho_hat - rho_target
    }

    # bracket in log-space

    lo <- log(1e-6)
    hi <- log(sigma_hi)

    # ensure bracket contains root

    f_lo <- f_obj(lo)
    f_hi <- f_obj(hi)
    if (f_lo < 0) {
        return(0)
    } # already below target even at tiny sigma
    if (f_hi > 0) {
        return(exp(hi))
    } # still above target even at huge sigma

    uniroot(f_obj, lower = lo, upper = hi)$root |> exp()
}

#' Given q, compute score cutoff s_cut so that Pr(S >= s_cut) = q
#' @export
score_cutoff <- function(S, q) {
    # top q => (1-q) quantile
    as.numeric(stats::quantile(S, probs = 1 - q, names = FALSE, type = 7))
}

#' Estimate p_pre and p_rem by Monte Carlo under calibrated sigma
#' @export
estimate_pre_rem_means <- function(lambda, P, sigma, q, seed = 1) {
    set.seed(seed)
    S <- lambda + rnorm(length(lambda), mean = 0, sd = sigma)

    if (q <= 0) {
        return(list(p_pre = NA_real_, p_rem = mean(P)))
    }
    if (q >= 1) {
        return(list(p_pre = mean(P), p_rem = NA_real_))
    }

    s_cut <- score_cutoff(S, q = q)
    sel <- S >= s_cut

    p_pre <- mean(P[sel])
    p_rem <- mean(P[!sel])

    list(p_pre = p_pre, p_rem = p_rem, s_cut = s_cut)
}

#' Mixed-strategy cost under Option 2
#' @export
cost_mix_heterorisk <- function(alpha, f, r, R, theta, rho, nu = 1,
                                lambda, P, seed_sigma = 1, seed_score = 2,
                                sigma_cache = NULL) {
    alpha <- max(0, min(1, alpha))
    q <- alpha * f

    # calibrate sigma if not provided

    sigma <- if (is.null(sigma_cache)) {
        calibrate_sigma_for_rho(lambda,
            rho_target = rho,
            seed = seed_sigma
        )
    } else {
        sigma_cache
    }

    # estimate remaining mean risk p' = p_rem(q, rho)
    # AND estimate pre-emptive mean risk p_pre for the cost calculation

    pre_rem <- estimate_pre_rem_means(lambda, P, sigma = sigma, q = q, seed = seed_score)
    p_prime <- pre_rem$p_rem
    p_pre_mean <- pre_rem$p_pre

    # effective reactive capacity among remaining group

    frac_rem <- 1 - q
    f_react <- (1 - alpha) * f
    f_prime <- f_react / frac_rem

    # pre-emptive cost (normalized by C_V)
    # With nu < 1, the pre-emptive group (size q) incurs vaccination cost (q)
    # PLUS expected outbreak cost: sum(P_i * (1-nu) * R) for i in pre-emptive
    # Expected outbreak cost = q * p_pre_mean * (1 - nu) * R

    # Handle q=0 case where p_pre_mean might be NA or undefined
    if (q < 1e-9) {
        c_pre <- 0
    } else {
        c_pre <- q + q * p_pre_mean * (1 - nu) * R
    }

    # remaining-group per-pop cost (same two regimes as equal-risk subproblem)
    # Note: Reactive group still has effectiveness r, so 'r' parameter is used here.

    if (f_prime < p_prime) {
        c_rem_per <- f_prime + (p_prime - f_prime * r) * R
    } else {
        c_rem_per <- p_prime * (1 + (1 - r) * R)
    }

    # If frac_rem is 0, then just c_pre
    if (frac_rem < 1e-9) {
        return(c_pre)
    }

    c_pre + frac_rem * c_rem_per
}

#' Grid-search for alpha*
#' @export
opt_alpha_heterorisk <- function(
    f, r, R, theta, rho, nu = 1,
    lambda, P,
    grid_len = 501,
    seed_sigma = 1, seed_score = 2) {
    # calibrate sigma once per scenario for speed/consistency

    sigma <- calibrate_sigma_for_rho(lambda, rho_target = rho, seed = seed_sigma)

    alpha_grid <- seq(0, 1, length.out = grid_len)
    cost_grid <- sapply(alpha_grid, function(a) {
        cost_mix_heterorisk(a,
            f = f, r = r, R = R, theta = theta, rho = rho, nu = nu,
            lambda = lambda, P = P, seed_sigma = seed_sigma,
            seed_score = seed_score, sigma_cache = sigma
        )
    })

    idx <- which.min(cost_grid)

    list(
        alpha_star = alpha_grid[idx],
        cost_star = cost_grid[idx],
        alpha_grid = alpha_grid,
        cost_grid = cost_grid,
        sigma = sigma
    )
}

#' Mixed strategy cost c_mix^{(n)}(alpha) = C_mix^{(n)} / C_V
#' @export
cost_mix_multi <- function(alpha, p, R, r, f, nu = 1) {
    alpha <- pmax(pmin(alpha, 1), 0)

    f_pre <- alpha * f
    f_react <- (1 - alpha) * f

    # Reactive-rich among non-pre-emptive populations
    cond_rich <- f_react >= p * (1 - f_pre)

    c_rich <- f_pre * (1 + p * (1 - nu) * R) + (1 - f_pre) * p * (1 + (1 - r) * R)
    c_lim <- f + R * (p - f_pre * p * nu - f_react * r)

    ifelse(cond_rich, c_rich, c_lim)
}

#' Numerical optimizer over a grid (useful for plotting / validation)
#' @export
opt_alpha_equalrisk <- function(p, R, r, f, nu = 1, grid_len = 1001) {
    alpha_grid <- seq(0, 1, length.out = grid_len)
    costs <- cost_mix_multi(alpha = alpha_grid, p = p, R = R, r = r, f = f, nu = nu)

    idx_min <- which.min(costs)

    list(
        alpha_star = alpha_grid[idx_min],
        cost_star  = costs[idx_min],
        alpha_grid = alpha_grid,
        cost_grid  = costs
    )
}

#' Closed-form alpha^* for equal-risk multi-population case
#' @export
alpha_star_equalrisk <- function(p, R, r, f, nu = 1, eps = 1e-12) {
    # r can be a vector

    # Scarce capacity: f <= p
    if (f <= p + eps) {
        # If r > p*nu -> pure reactive (0)
        # If r < p*nu -> pure pre-emptive (1)
        return(
            ifelse(r > p * nu + eps, 0,
                ifelse(r < p * nu - eps, 1, NA_real_)
            )
        )
    }

    # Abundant capacity: f > p
    alpha_c <- (f - p) / (f * (1 - p))

    alpha_star <- rep(NA_real_, length(r))

    # r < p*nu -> pure pre-emptive
    alpha_star[r < p * nu - eps] <- 1

    # r > p*nu -> compare alpha=0 vs alpha=alpha_c
    idx <- which(r > p * nu + eps)
    if (length(idx) > 0) {
        # R_thr = (1-p) / (p*(nu - r)) if nu > r
        # If nu <= r, slope is positive => pure reactive (alpha=0)

        # We set alpha_star to 0 by default for these indices
        alpha_star[idx] <- 0

        # Check where nu > r (slope could be negative)
        # Actually, R_thr is valid only if nu > r.
        # If nu <= r, then 1-p + pR(r-nu) > 0 always (since 1-p>0, p>0, R>0, r>=nu).
        # So if nu <= r, pure reactive is always better.

        # We only consider switching to alpha_c if nu > r AND R >= R_thr
        sub_idx <- idx[nu > r[idx] + eps]

        if (length(sub_idx) > 0) {
            # Calculate R_thr for these
            R_vals_sub <- if (length(R) == 1) rep(R, length(sub_idx)) else R[sub_idx]

            R_thr <- (1 - p) / (p * (nu - r[sub_idx]))

            # If R >= R_thr, switch to alpha_c
            switch_mask <- R_vals_sub >= R_thr - eps

            # Map sub_idx back to alpha_star
            alpha_star[sub_idx[switch_mask]] <- alpha_c
        }
    }

    alpha_star
}

#' Stochastic simulation of per-population cost under heterogeneous risk
#'
#' Draws n populations with risks from Beta(1, theta), simulates outbreak
#' occurrence as Bernoulli draws, allocates pre-emptive and reactive
#' vaccination, and computes realized per-population cost.
#'
#' @param alpha Numeric vector of pre-emptive fractions to evaluate.
#' @param n Number of populations (units).
#' @param p Mean outbreak probability (determines theta = (1-p)/p).
#' @param f Total vaccine capacity as fraction of n.
#' @param R Cost ratio C_I / C_V.
#' @param r Reactive vaccination effectiveness.
#' @param nu Pre-emptive vaccine efficacy (default 1).
#' @param rho Spearman rank correlation for targeting accuracy (1 = perfect).
#' @param n_sim Number of Monte Carlo replications per alpha.
#' @param seed Base random seed.
#' @return A data.frame with columns: sim, alpha, cost.
#' @export
sim_cost_heterorisk <- function(alpha, n, p, f, R, r, nu = 1, rho = 1,
                                n_sim = 50, seed = 42) {
    theta <- (1 - p) / p

    results <- vector("list", length(alpha) * n_sim)
    idx <- 1L

    for (s in seq_len(n_sim)) {
        # Draw population risks (fresh draw per simulation)
        set.seed(seed + s)
        x <- rbeta(n, shape1 = 1, shape2 = theta)

        # Targeting: rank populations by score
        if (rho >= 1) {
            score <- x
        } else {
            # Use existing calibration for imperfect targeting
            lambda <- -log(1 - x) # transform to exponential scale
            sigma <- calibrate_sigma_for_rho(lambda, rho_target = rho,
                                             seed = seed + s)
            noise <- rnorm(n, mean = 0, sd = sigma)
            score <- lambda + noise
        }
        rank_order <- order(score, decreasing = TRUE)

        for (a in alpha) {
            n_pre <- floor(n * a * f)
            n_react <- floor(n * (1 - a) * f)

            # Pre-emptive vaccination: top n_pre by score
            pre_vax <- rep(FALSE, n)
            if (n_pre > 0) {
                pre_vax[rank_order[seq_len(n_pre)]] <- TRUE
            }

            # Draw outbreaks: Bernoulli(x_i) for each population
            Y <- rbinom(n, size = 1, prob = x)

            # Reactive allocation: among non-pre-emptive with outbreak
            react_vax <- rep(FALSE, n)
            react_eligible <- which(!pre_vax & Y == 1)
            if (length(react_eligible) > 0 && n_react > 0) {
                n_to_react <- min(n_react, length(react_eligible))
                react_selected <- react_eligible[
                    sample.int(length(react_eligible), n_to_react)
                ]
                react_vax[react_selected] <- TRUE
            }

            # Per-population cost (normalized by C_V)
            cost_i <- numeric(n)
            # Pre-emptive vaccinated: pay vaccine + possible breakthrough
            cost_i[pre_vax] <- 1 + Y[pre_vax] * (1 - nu) * R
            # Reactive vaccinated: pay vaccine + reduced illness
            cost_i[react_vax] <- 1 + (1 - r) * R
            # Unmitigated outbreak: full illness cost
            unmitigated <- !pre_vax & !react_vax & Y == 1
            cost_i[unmitigated] <- R

            results[[idx]] <- data.frame(
                sim = s, alpha = a, cost = mean(cost_i)
            )
            idx <- idx + 1L
        }
    }

    do.call(rbind, results)
}


#' Fraction of area remaining under an exponential tent curve after time t
#'
#' The outbreak follows exponential growth from t0 to tp and exponential
#' decay from tp to t1 (tent-shaped on a log scale). Growth and decay rates
#' are specified directly.
#'
#' @param t Time at which to evaluate (absolute).
#' @param t0 Outbreak start time.
#' @param tp Outbreak peak time.
#' @param t1 Outbreak end time.
#' @param rg Exponential growth rate.
#' @param rd Exponential decay rate.
#' @return Fraction of total cumulative incidence occurring after time t.
frac_remaining_exp <- function(t, t0, tp, t1, rg, rd) {
    if (t >= t1) return(0)
    if (t <= t0) return(1)

    # I(t0)/Ip = exp(-rg*(tp-t0)),  I(t1)/Ip = exp(-rd*(t1-tp))
    # Normalize Ip = 1 (cancels in ratio)
    I_t0 <- exp(-rg * (tp - t0))
    I_t1 <- exp(-rd * (t1 - tp))

    # Total area: A(t0, t1) = (1 - I_t0)/rg + (1 - I_t1)/rd
    total <- (1 - I_t0) / rg + (1 - I_t1) / rd

    # Area from t0 to t
    if (t <= tp) {
        I_t <- exp(rg * (t - tp))
        area_before <- (I_t - I_t0) / rg
    } else {
        I_t <- exp(-rd * (t - tp))
        area_before <- (1 - I_t0) / rg + (1 - I_t) / rd
    }

    (total - area_before) / total
}


#' Stochastic simulation with endogenous reactive effectiveness
#'
#' Extends \code{sim_cost_heterorisk()} by replacing the fixed parameter r with
#' a per-population reactive effectiveness derived from outbreak trajectories,
#' response delay, and vaccine effectiveness (VE_d and VE_i).
#'
#' For each reactively vaccinated population i:
#' \deqn{r_i = \text{OVE} \times \text{frac\_remaining}(t_{\text{eff},i})}
#' where OVE = 1 - [f_cov(1-VE_d)(1-VE_i) + (1-f_cov)(1-VE_i)] and
#' frac_remaining is the fraction of cumulative incidence after the vaccine
#' becomes effective under an exponential tent-shaped outbreak curve.
#'
#' @param alpha Numeric vector of pre-emptive fractions.
#' @param n Number of populations.
#' @param p Mean outbreak probability (Beta(1, theta) with theta = (1-p)/p).
#' @param f Total vaccine capacity as fraction of n.
#' @param R Cost ratio C_I / C_V.
#' @param VEd Direct vaccine effectiveness.
#' @param VEi Indirect vaccine effectiveness.
#' @param f_cov Vaccination coverage within a targeted population.
#' @param nu Pre-emptive vaccine efficacy (default 1).
#' @param rho Spearman rank correlation for targeting (1 = perfect).
#' @param duration_range Range for outbreak duration (days).
#' @param peak_frac_shape Beta distribution parameters for peak fraction.
#' @param rg_range Range for exponential growth rate.
#' @param rd_range Range for exponential decay rate.
#' @param delay_range Range for response delay (days).
#' @param n_sim Number of Monte Carlo replications per alpha.
#' @param seed Base random seed.
#' @return A data.frame with columns: sim, alpha, cost, mean_r.
#' @export
sim_cost_heterorisk_endogenous <- function(alpha, n, p, f, R,
                                           VEd = 0.7, VEi = 0.3,
                                           f_cov = 0.7,
                                           nu = 1, rho = 1,
                                           duration_range = c(60, 300),
                                           peak_frac_shape = c(2, 2),
                                           rg_range = c(0.02, 0.10),
                                           rd_range = c(0.01, 0.05),
                                           delay_range = c(7, 21),
                                           n_sim = 50, seed = 42) {
    theta <- (1 - p) / p

    # Overall vaccine effectiveness
    OVE <- 1 - (1 - VEi) * (1 - f_cov * VEd)

    results <- vector("list", length(alpha) * n_sim)
    idx <- 1L

    for (s in seq_len(n_sim)) {
        set.seed(seed + s)
        x <- rbeta(n, shape1 = 1, shape2 = theta)

        # Targeting
        if (rho >= 1) {
            score <- x
        } else {
            lambda <- -log(1 - x)
            sigma <- calibrate_sigma_for_rho(lambda, rho_target = rho,
                                             seed = seed + s)
            noise <- rnorm(n, mean = 0, sd = sigma)
            score <- lambda + noise
        }
        rank_order <- order(score, decreasing = TRUE)

        for (a in alpha) {
            n_pre <- floor(n * a * f)
            n_react <- floor(n * (1 - a) * f)

            # Pre-emptive vaccination
            pre_vax <- rep(FALSE, n)
            if (n_pre > 0) {
                pre_vax[rank_order[seq_len(n_pre)]] <- TRUE
            }

            # Outbreak draws
            Y <- rbinom(n, size = 1, prob = x)

            # Reactive allocation
            react_vax <- rep(FALSE, n)
            react_eligible <- which(!pre_vax & Y == 1)
            if (length(react_eligible) > 0 && n_react > 0) {
                n_to_react <- min(n_react, length(react_eligible))
                react_selected <- react_eligible[
                    sample.int(length(react_eligible), n_to_react)
                ]
                react_vax[react_selected] <- TRUE
            }

            # Compute per-population r_i for reactive-vaccinated populations
            r_i <- numeric(n)
            react_idx <- which(react_vax)
            if (length(react_idx) > 0) {
                for (j in react_idx) {
                    # Draw outbreak trajectory
                    D <- runif(1, duration_range[1], duration_range[2])
                    phi <- rbeta(1, peak_frac_shape[1], peak_frac_shape[2])
                    t0_j <- 0
                    tp_j <- phi * D
                    t1_j <- D
                    rg_j <- runif(1, rg_range[1], rg_range[2])
                    rd_j <- runif(1, rd_range[1], rd_range[2])

                    # Response delay
                    delay <- runif(1, delay_range[1], delay_range[2])
                    t_eff <- t0_j + delay

                    frac_rem <- frac_remaining_exp(t_eff, t0_j, tp_j, t1_j,
                                                   rg_j, rd_j)
                    r_i[j] <- OVE * frac_rem
                }
            }

            # Mean r across reactively vaccinated populations
            mean_r <- if (length(react_idx) > 0) {
                mean(r_i[react_idx])
            } else {
                NA_real_
            }

            # Per-population cost
            cost_i <- numeric(n)
            cost_i[pre_vax] <- 1 + Y[pre_vax] * (1 - nu) * R
            cost_i[react_vax] <- 1 + (1 - r_i[react_vax]) * R
            unmitigated <- !pre_vax & !react_vax & Y == 1
            cost_i[unmitigated] <- R

            results[[idx]] <- data.frame(
                sim = s, alpha = a, cost = mean(cost_i), mean_r = mean_r
            )
            idx <- idx + 1L
        }
    }

    do.call(rbind, results)
}

# ── Cost-ratio utilities ───────────────────────────────────────────────────────

#' Total economic cost per cholera case.
#'
#' Sums direct medical costs (outpatient + hospital, weighted by case severity)
#' and indirect costs (morbidity DALYs, premature-death YLLs, productivity
#' losses) per reported case.
#'
#' All monetary parameters are in USD. \code{year_val} is the value of one
#' statistical life-year (e.g., GDP per capita or 3 × GDP per capita).
#'
#' @param day_ill            Mean illness duration (days).
#' @param mean_dis_wt        Disability weight (0–1).
#' @param year_val           Value of one life-year (USD).
#' @param mean_prop_workforce Fraction of population in workforce.
#' @param patient_workday_lost Patient workdays lost per episode.
#' @param caregiver_workday_lost Caregiver workdays lost per episode.
#' @param mean_cfr           Case-fatality ratio.
#' @param mean_remaining_life Remaining life expectancy at mean age of infection (years).
#' @param pr_moderate        Proportion of cases that are moderate (outpatient).
#' @param pr_severe          Proportion of cases that are severe (hospitalised).
#' @param patient_cost_outpt Patient direct cost for outpatient case (USD).
#' @param patient_cost_hosp  Patient direct cost for hospitalised case (USD).
#' @param public_cost_outpt  Public-sector cost for outpatient case (USD).
#' @param public_cost_hosp   Public-sector cost for hospitalised case (USD).
#' @return Scalar: total cost per case (USD).
#' @export
total_cost_per_case <- function(day_ill              = 2,
                                mean_dis_wt          = 0.188,
                                year_val             = 1705.644,
                                mean_prop_workforce  = 0.6805938,
                                patient_workday_lost = 6.5,
                                caregiver_workday_lost = 3.3,
                                mean_cfr             = 0.01,
                                mean_remaining_life  = 42.63315,
                                pr_moderate          = 0.78,
                                pr_severe            = 0.22,
                                patient_cost_outpt   = 5.71,
                                patient_cost_hosp    = 38.87,
                                public_cost_outpt    = 2.84,
                                public_cost_hosp     = 65.77) {
    pr_tot <- pr_moderate + pr_severe

    direct_cost <- pr_moderate / pr_tot * (patient_cost_outpt + public_cost_outpt) +
                   pr_severe   / pr_tot * (patient_cost_hosp  + public_cost_hosp)

    indirect_morbidity  <- (day_ill / 365) * mean_dis_wt * year_val
    indirect_mortality  <- mean_cfr * mean_remaining_life * year_val
    indirect_productivity <- mean_prop_workforce *
                             (patient_workday_lost + caregiver_workday_lost) / 365 * year_val

    direct_cost + indirect_morbidity + indirect_mortality + indirect_productivity
}

#' Compute total cost for one simulation realization (cell-level).
#'
#' Merges cell, trajectory, and reactive-vaccination tables, applies
#' pre-emptive and reactive fractional reductions to the attack rate
#' multiplicatively, then sums outbreak, pre-emptive vaccination, and
#' reactive vaccination costs across all cells.
#'
#' Requires \pkg{data.table}.
#'
#' @param cells_dt    data.frame/data.table with columns:
#'   \code{id}, \code{pop_size}, \code{pre_vax} (logical), \code{outbreak} (logical).
#' @param traj_dt     data.frame/data.table with columns:
#'   \code{id}, \code{attack_rate_per_capita}.
#' @param reactive_dt data.frame/data.table with columns:
#'   \code{id}, \code{reactive_vax} (logical); optionally
#'   \code{frac_reduction_ci} and \code{rv_coverage}.
#' @param C_case               Cost per case (USD).
#' @param C_vac_per_person     Vaccination cost per person (USD).
#' @param pre_coverage         Within-cell coverage for pre-emptive vaccination.
#' @param default_reactive_coverage Within-cell coverage for reactive vaccination
#'   when \code{rv_coverage} is absent.
#' @param reactive_coverage_col Column name for per-cell reactive coverage.
#' @param react_frac_red_col   Column name for reactive fractional reduction.
#' @param pre_frac_red_col     Column name for pre-emptive fractional reduction.
#' @param nu                   Vaccine efficacy (used only when
#'   \code{pre_frac_red_col} is missing/NA).
#' @return A list with elements \code{totals} (named numeric vector) and
#'   \code{breakdown} (data.table with per-cell costs).
#' @export
compute_total_cost <- function(cells_dt,
                               traj_dt,
                               reactive_dt,
                               C_case,
                               C_vac_per_person,
                               pre_coverage              = 0.9,
                               default_reactive_coverage = 0.9,
                               reactive_coverage_col     = "rv_coverage",
                               react_frac_red_col        = "frac_reduction_ci",
                               pre_frac_red_col          = "pre_frac_reduction_ci",
                               nu                        = 0.9) {
    cells <- data.table::as.data.table(cells_dt)
    traj  <- data.table::as.data.table(traj_dt)
    rv    <- data.table::as.data.table(reactive_dt)

    req_cells <- c("id", "pop_size", "pre_vax", "outbreak")
    miss <- setdiff(req_cells, names(cells))
    if (length(miss) > 0L)
        stop("cells_dt missing columns: ", paste(miss, collapse = ", "))
    if (!("id" %in% names(traj)) || !("attack_rate_per_capita" %in% names(traj)))
        stop("traj_dt must contain 'id' and 'attack_rate_per_capita'.")
    if (!("id" %in% names(rv)) || !("reactive_vax" %in% names(rv)))
        stop("reactive_dt must contain 'id' and 'reactive_vax'.")

    dt <- merge(cells, traj[, .(id, attack_rate_per_capita)], by = "id", all.x = TRUE)
    dt <- merge(dt, rv, by = "id", all.x = TRUE)

    dt[, pre_vax     := as.logical(pre_vax)]
    dt[, outbreak    := as.logical(outbreak)]
    dt[, reactive_vax := as.logical(reactive_vax)]
    dt[is.na(attack_rate_per_capita), attack_rate_per_capita := 0]

    # Pre-emptive fractional reduction
    pre_coverage    <- pmin(pmax(pre_coverage, 0), 1)
    pre_red_default <- pmin(1, nu * pre_coverage)
    if (!(pre_frac_red_col %in% names(dt))) dt[, (pre_frac_red_col) := NA_real_]
    dt[, pre_red_used := data.table::fifelse(
        pre_vax,
        data.table::fifelse(is.na(get(pre_frac_red_col)), pre_red_default,
                            get(pre_frac_red_col)),
        0.0)]
    dt[, pre_red_used := pmin(pmax(pre_red_used, 0), 1)]

    # Reactive fractional reduction
    if (!(react_frac_red_col %in% names(dt))) dt[, (react_frac_red_col) := NA_real_]
    dt[, (react_frac_red_col) := suppressWarnings(as.numeric(get(react_frac_red_col)))]
    dt[, react_red_used := 0.0]
    dt[reactive_vax == TRUE & !is.na(get(react_frac_red_col)),
       react_red_used := get(react_frac_red_col)]
    dt[, react_red_used := pmin(pmax(react_red_used, 0), 1)]

    # Effective attack rate (multiplicative reductions)
    dt[, attack_rate_eff :=
           data.table::fifelse(outbreak,
               attack_rate_per_capita * (1 - pre_red_used) * (1 - react_red_used),
               0)]

    # Per-cell costs
    dt[, C_outbreak := pop_size * attack_rate_eff * C_case]
    dt[, C_pre      := data.table::fifelse(pre_vax,
                           pop_size * pre_coverage * C_vac_per_person, 0)]

    if (!(reactive_coverage_col %in% names(dt))) dt[, (reactive_coverage_col) := NA_real_]
    dt[, rv_cov_used := data.table::fifelse(
            is.na(get(reactive_coverage_col)), default_reactive_coverage,
            get(reactive_coverage_col))]
    dt[, rv_cov_used  := pmin(pmax(rv_cov_used, 0), 1)]
    dt[, C_reactive   := data.table::fifelse(reactive_vax,
                             pop_size * rv_cov_used * C_vac_per_person, 0)]

    totals <- list(
        outbreak_cost     = sum(dt$C_outbreak,  na.rm = TRUE),
        pre_vax_cost      = sum(dt$C_pre,       na.rm = TRUE),
        reactive_vax_cost = sum(dt$C_reactive,  na.rm = TRUE)
    )
    totals$total_cost <- totals$outbreak_cost + totals$pre_vax_cost +
                         totals$reactive_vax_cost

    list(
        totals    = totals,
        breakdown = dt[, .(id, pop_size, pre_vax, reactive_vax, outbreak,
                           attack_rate_per_capita,
                           pre_red_used, react_red_used, attack_rate_eff,
                           C_outbreak, C_pre, C_reactive,
                           C_total = C_outbreak + C_pre + C_reactive)]
    )
}
