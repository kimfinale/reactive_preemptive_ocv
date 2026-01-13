#' Pre-emptive cost (normalized by C_V)
#' @export
cost_pre_one <- function(p, R, r, nu = 1) {
    # p and r unused in the original (nu=1) case but kept for interface.
    # With nu < 1: Cost = 1 + p * (1 - nu) * R
    1 + p * (1 - nu) * R
}

#' Reactive cost (normalized by C_V)
#' @export
cost_react_one <- function(p, R, r) {
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
cost_react_multi <- function(p, R, r, f) {
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
