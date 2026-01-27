#' SIR model with simple vaccination
#'
#' @param time_step Numeric, time step for integration (e.g. 0.1 days)
#' @param duration Numeric, total simulation duration (days)
#' @param N Numeric, Total population size
#' @param I0 Numeric, Initial infected individuals
#' @param R0 Numeric, Basic reproduction number
#' @param gamma Numeric, Recovery rate (1/infectious_period)
#' @param rho Numeric, Vaccination rate (proportion of S vaccinated per time unit) - specific to campaign duration
#' @param vacc_start Numeric, Start time of vaccination campaign
#' @param vacc_duration Numeric, Duration of vaccination campaign
#' @param ve Numeric, Vaccine efficacy (0 to 1) applied to Susceptibles moving to Recovered/Vaccinated
#'
#' @return Data frame with time course of compartments and cumulative cases
#' @export
#' SIR model equations for deSolve
#'
#' @param time Current time
#' @param state Named vector of state variables
#' @param parameters Named vector of parameters
#' @return List of derivatives
sir_equations <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
        # Check if within vaccination campaign
        is_campaign <- (time >= vacc_start) && (time < (vacc_start + vacc_duration))

        # Force of infection
        lambda <- beta * I / N

        dS <- -lambda * S
        dI <- lambda * S - gamma * I
        dR <- gamma * I
        dCumCases <- lambda * S

        if (is_campaign) {
            # Vaccination moves people from S to R
            # Rate rho is applied to current S
            # Effective vaccination: rho * S * ve
            rate_vax <- rho * S * ve
            dS <- dS - rate_vax
            dR <- dR + rate_vax
        }

        list(c(dS, dI, dR, dCumCases))
    })
}

#' SIR model with simple vaccination (ODE version)
#'
#' @param time_step Numeric, time step for integration (and output)
#' @param duration Numeric, total simulation duration (days)
#' @param N Numeric, Total population size
#' @param I0 Numeric, Initial infected individuals
#' @param R0 Numeric, Basic reproduction number
#' @param gamma Numeric, Recovery rate (1/infectious_period)
#' @param coverage Numeric, Target vaccine coverage (proportion of population)
#' @param vacc_start Numeric, Start time of vaccination campaign
#' @param vacc_duration Numeric, Duration of vaccination campaign
#' @param ve Numeric, Vaccine efficacy (0 to 1) applied to Susceptibles moving to Recovered/Vaccinated
#'
#' @return Data frame with time course of compartments and cumulative cases
#' @export
run_sir_vaccination <- function(time_step = 1, duration = 365, N = 100000,
                                I0 = 1, R0 = 2.0, gamma = 1 / 5,
                                coverage = 0.8, vacc_start = 20,
                                vacc_duration = 10, ve = 0.7) {
    # Ensure deSolve is available; if not, you might need to install it.
    if (!requireNamespace("deSolve", quietly = TRUE)) {
        stop("Package 'deSolve' is required for this function.")
    }

    # Helper function for SIR equations with S -> R vaccination
    sir_equations_internal <- function(time, state, parameters) {
        with(as.list(c(state, parameters)), {
            # Check if within vaccination campaign
            is_campaign <- (time >= vacc_start) && (time < (vacc_start + vacc_duration))

            # Force of infection
            lambda <- beta * I / N

            dS <- -lambda * S
            dI <- lambda * S - gamma * I
            dR <- gamma * I
            dCumCases <- lambda * S

            if (is_campaign) {
                # Vaccination moves people from S to R (effectively)
                # The rate rho accounts for both coverage and vaccine efficacy
                rate_vax <- rho * S
                dS <- dS - rate_vax
                dR <- dR + rate_vax
            }

            list(c(dS, dI, dR, dCumCases))
        })
    }

    # Calculate effective coverage (proportion of S to move to R)
    effective_coverage <- coverage * ve

    # Calculate daily vaccination rate (rho) to achieve effective_coverage
    # Formula: rho = -log(1 - effective_coverage) / vacc_duration
    if (vacc_duration > 0 && effective_coverage > 0 && effective_coverage < 1) {
        rho <- -log(1 - effective_coverage) / vacc_duration
    } else {
        rho <- 0
    }

    beta <- R0 * gamma
    times <- seq(0, duration, by = time_step)

    # Initial state
    state <- c(
        S = N - I0,
        I = I0,
        R = 0,
        CumCases = I0
    )

    parameters <- c(
        beta = beta,
        gamma = gamma,
        rho = rho,
        vacc_start = vacc_start,
        vacc_duration = vacc_duration,
        ve = ve, # Kept for parameter completeness, though absorbed into rho
        N = N
    )

    # Solve ODE using Runge-Kutta 4 (rk4) with the internal equation function
    out <- deSolve::ode(
        y = state,
        times = times,
        func = sir_equations_internal,
        parms = parameters,
        method = "rk4"
    )

    as.data.frame(out)
}

#' Calculate Proportional Reduction (r) using ODE
#'
#' @param ... Arguments passed to run_sir_vaccination
#' @return List with r value, and data frames for no-vax and vax scenarios
#' @export
calculate_r_ode <- function(N = 100000, I0 = 1, R0 = 2.0, gamma = 1 / 5,
                            vacc_start = 20, coverage = 0.8, vacc_duration = 10, ve = 0.7, ...) {
    # Run No Vaccination Baseline (coverage = 0)
    sim_novax <- run_sir_vaccination(
        N = N, I0 = I0, R0 = R0, gamma = gamma,
        coverage = 0, vacc_start = 0, vacc_duration = 0, ve = 0, ...
    )

    # Run With Vaccination
    sim_vax <- run_sir_vaccination(
        N = N, I0 = I0, R0 = R0, gamma = gamma,
        coverage = coverage, vacc_start = vacc_start,
        vacc_duration = vacc_duration, ve = ve, ...
    )

    cases_novax <- tail(sim_novax$CumCases, 1)
    cases_vax <- tail(sim_vax$CumCases, 1)

    r <- 1 - (cases_vax / cases_novax)

    list(
        r = r, cases_novax = cases_novax, cases_vax = cases_vax,
        sim_novax = sim_novax, sim_vax = sim_vax
    )
}

#' Calculate cases using overall effectiveness approach
#'
#' A simple static model: Cases = Cases_Baseline * (1 - Overallwise_Effectiveness)
#' Overall Effectiveness depends on coverage, direct vaccine effectiveness (DVE), and indirect vaccine effectiveness (IVE).
#' Formula: OVE = 1 - (coverage * (1 - DVE) * (1 - IVE) + (1 - coverage) * (1 - IVE))
#'
#' @param cases_baseline Numeric, estimated cases without vaccination
#' @param coverage Numeric, proportion of population vaccinated (f)
#' @param ve_direct Numeric, direct vaccine efficacy (DVE)
#' @param ve_indirect Numeric, indirect vaccine effectiveness (IVE)
#'
#' @return Numeric, estimated cases after vaccination
#' @export
#' Calculate cases using overall effectiveness approach
#'
#' Adjusts a baseline incidence time series by applying an overall effectiveness (OVE) reduction
#' to cases occurring after the start of a vaccination campaign.
#'
#' @param incidence_baseline Numeric vector, daily (or per-step) incidence time series without vaccination
#' @param time Numeric vector, time points corresponding to incidence_baseline
#' @param coverage Numeric, proportion of population vaccinated (f) at the end of campaign
#' @param ve_direct Numeric, direct vaccine efficacy (DVE)
#' @param ve_indirect Numeric, indirect vaccine effectiveness (IVE)
#' @param vacc_start Numeric, Time when vaccination reduction begins to apply
#'
#' @return List with updated cumulative cases time series `cases_post`, `incidence_post` and the calculated `ove`
#' @export
calculate_cases_overall_effectiveness <- function(incidence_baseline, time, coverage, ve_direct,
                                                  ve_indirect = 0, vacc_start = 0) {
    # OVE = 1 - (f * (1 - DVE) * (1 - IVE) + (1 - f) * (1 - IVE))
    term_vaccinated <- coverage * (1 - ve_direct) * (1 - ve_indirect)
    term_unvaccinated <- (1 - coverage) * (1 - ve_indirect)

    relative_risk_overall <- term_vaccinated + term_unvaccinated
    ove <- 1 - relative_risk_overall

    # Apply OVE to reducing incidence occurring AFTER vaccination start
    incidence_post <- incidence_baseline

    idx_after <- time > vacc_start
    if (any(idx_after)) {
        incidence_post[idx_after] <- incidence_baseline[idx_after] * (1 - ove)
    }

    cases_post <- cumsum(incidence_post)

    return(list(cases_post = cases_post, incidence_post = incidence_post, ove = ove))
}


#' Stratified SIR model (2 groups) - Pre-emptive (No dynamic vaccination)
#'
#' @param time Current time
#' @param state Named vector of state variables
#' @param parameters Named vector of parameters
#' @return List of derivatives
#' @export
sir_2grp_preemptive <- function(time, state, parameters) {
    # Unpack state
    S0 <- state["S0"]
    I0 <- state["I0"]
    R0 <- state["R0"]
    C0 <- state["C0"]
    S1 <- state["S1"]
    I1 <- state["I1"]
    R1 <- state["R1"]
    C1 <- state["C1"]

    # Unpack parameters
    beta <- parameters["beta"]
    gamma <- parameters["gamma"]

    # Total population (sum of all compartments)
    N_total <- S0 + I0 + R0 + S1 + I1 + R1

    # Force of infection (assuming well-mixed population)
    lambda <- beta * (I0 + I1) / N_total

    # Group 0 (Non-Target)
    dS0 <- -lambda * S0
    dI0 <- lambda * S0 - gamma * I0
    dR0 <- gamma * I0
    dC0 <- lambda * S0

    # Group 1 (Target)
    dS1 <- -lambda * S1
    dI1 <- lambda * S1 - gamma * I1
    dR1 <- gamma * I1
    dC1 <- lambda * S1

    list(c(dS0, dI0, dR0, dC0, dS1, dI1, dR1, dC1))
}

#' Calculate Indirect Vaccine Effectiveness (IVE) from ODE (Instantaneous Pre-emptive)
#'
#' @param N Numeric, Total population size
#' @param I0 Numeric, Initial infected individuals
#' @param R0 Numeric, Basic reproduction number
#' @param gamma Numeric, Recovery rate
#' @param ve Numeric, Vaccine efficacy
#' @param coverage Numeric, Vaccine coverage (Target proportion of total population)
#' @return List with time series of IVE and final IVE value
#' @export
calculate_ive_ode <- function(N = 100000, I0 = 1, R0 = 2.0, gamma = 1 / 5,
                              ve = 0.7, coverage = 0.8, time_step = 1,
                              duration = 365) {
    # 1. Setup Parameters
    beta <- R0 * gamma
    params <- c(beta = beta, gamma = gamma)

    # --- Population Groups ---
    # Group 0: Non-target (remains unvaccinated).
    # Group 1: Target (receives vaccination).
    # IVE is calculated by comparing the attack rate in the non-target group (Group 0)
    # between the vaccinated and baseline scenarios.
    N1 <- N * coverage
    N0 <- N - N1
    # I0 is assumed to be in the target group, which makes it easy to compute the IVE based on the non-target group
    I0_1 <- I0

    # --- Scenario B: Baseline (No Vax) ---
    # R initialized to 0
    y0_base <- c(
        S0 = N0, I0 = 0, R0 = 0, C0 = 0,
        S1 = N1 - I0_1, I1 = I0_1, R1 = 0, C1 = 0
    )

    times <- seq(0, duration, by = time_step)

    out_novax <- deSolve::ode(y = y0_base, times = times, func = sir_2grp_preemptive, parms = params)
    df_novax <- as.data.frame(out_novax)

    # --- Scenario A: Vaccination ---

    N1_susceptible <- N1 - I0_1
    R1_init <- N1_susceptible * ve
    S1_init <- N1_susceptible * (1 - ve)

    y0_vax <- c(
        S0 = N0, I0 = 0, R0 = 0, C0 = 0,
        S1 = S1_init, I1 = I0_1, R1 = R1_init, C1 = 0
    )

    out_vax <- deSolve::ode(y = y0_vax, times = times, func = sir_2grp_preemptive, parms = params)
    df_vax <- as.data.frame(out_vax)

    # --- Calculate IVE ---
    # Compare AR in Group 0 (Never Vaccinated / Non-Target)
    # AR0 = C0 / N0

    # Baseline AR0
    C0_novax <- df_novax$C0
    AR0_novax <- rep(0, length(C0_novax))
    if (N0 > 0) AR0_novax <- C0_novax / N0

    # Vax Scenario AR0
    C0_vax <- df_vax$C0
    AR0_vax <- rep(0, length(C0_vax))
    if (N0 > 0) AR0_vax <- C0_vax / N0

    # IVE vector
    ive_series <- numeric(length(AR0_novax))
    valid <- AR0_novax > 1e-9
    ive_series[valid] <- 1 - (AR0_vax[valid] / AR0_novax[valid])

    final_ive <- tail(ive_series, 1)

    list(
        ive = final_ive,
        ive_series = ive_series,
        time = times,
        df_novax = df_novax,
        df_vax = df_vax
    )
}
