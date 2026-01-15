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
run_sir_vaccination <- function(time_step = 0.1, duration = 365, N = 100000,
                                I0 = 1, R0 = 2.0, gamma = 1 / 5,
                                rho = 0.01, vacc_start = 20,
                                vacc_duration = 10, ve = 0.7) {
    beta <- R0 * gamma / N
    times <- seq(0, duration, by = time_step)
    n_steps <- length(times)

    S <- numeric(n_steps)
    I <- numeric(n_steps)
    R <- numeric(n_steps)
    V <- numeric(n_steps) # Vaccinated compartment (effectively Removed)
    CumCases <- numeric(n_steps)

    S[1] <- N - I0
    I[1] <- I0
    R[1] <- 0
    V[1] <- 0
    CumCases[1] <- I0

    for (t in 2:n_steps) {
        curr_time <- times[t - 1]

        # Vaccination logic
        # Assume a constant rate `rho` of the REMAINING Susceptibles are vaccinated
        # if we are within the campaign window.
        # Alternatively, `rho` could be a fraction of N per day.
        # Let's assume rho is fraction of N per day target, but limited by available S.
        # vacc_rate = rho * N (doses per day)
        # effective_vacc_rate = min(vacc_rate, S) * ve
        # Here, let's keep it simple: proportion of current S moving to V

        is_campaign <- (curr_time >= vacc_start) && (curr_time < (vacc_start + vacc_duration))

        # Force of infection
        lambda <- beta * I[t - 1]

        # Derivatives (approximate with Euler for simplicity, or use simple discrete step)
        new_inf <- lambda * S[t - 1] * time_step
        new_rec <- gamma * I[t - 1] * time_step

        new_vac <- 0
        if (is_campaign) {
            # Rate rho applied to S
            new_vac <- rho * S[t - 1] * time_step * ve
        }

        # Updates
        S[t] <- S[t - 1] - new_inf - new_vac
        I[t] <- I[t - 1] + new_inf - new_rec
        R[t] <- R[t - 1] + new_rec
        V[t] <- V[t - 1] + new_vac
        CumCases[t] <- CumCases[t - 1] + new_inf

        # Ensure non-negativity
        if (S[t] < 0) S[t] <- 0
        if (I[t] < 0) I[t] <- 0
    }

    data.frame(time = times, S = S, I = I, R = R, V = V, CumCases = CumCases)
}

#' Calculate Proportional Reduction (r) using ODE
#'
#' @param ... Arguments passed to run_sir_vaccination
#' @return List with r value, and data frames for no-vax and vax scenarios
#' @export
calculate_r_ode <- function(N = 100000, I0 = 1, R0 = 2.0, gamma = 1 / 5,
                            vacc_start = 20, rho = 0.05, vacc_duration = 10, ve = 0.7, ...) {
    # Run No Vaccination Baseline (rho = 0)
    sim_novax <- run_sir_vaccination(
        N = N, I0 = I0, R0 = R0, gamma = gamma,
        rho = 0, vacc_start = 0, vacc_duration = 0, ve = 0, ...
    )

    # Run With Vaccination
    sim_vax <- run_sir_vaccination(
        N = N, I0 = I0, R0 = R0, gamma = gamma,
        rho = rho, vacc_start = vacc_start,
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

#' Calculate Cases using Overall Effectiveness Approach
#'
#' A simple static model: Cases = Cases_Baseline * (1 - Overallwise_Effectiveness)
#' Overall Effectiveness depends on coverage and possibly indirect effects.
#'
#' @param cases_baseline Numeric, estimated cases without vaccination
#' @param coverage Numeric, proportion of population vaccinated
#' @param ve_direct Numeric, direct vaccine efficacy
#' @param ve_indirect_func Function, takes coverage and returns indirect protection factor (0-1)
#'
#' @return Numeric, estimated cases after vaccination
#' @export
calculate_cases_overall_effectiveness <- function(cases_baseline, coverage, ve_direct,
                                                  ve_indirect_func = NULL) {
    # Direct effect: vaccinated individuals are protected
    # Simple approximation: weighted average of risk
    # Risk_vax = Risk_unvax * (1 - ve_direct)
    # This simple formula assumes determining effectiveness on the WHOLE population:
    # Overall_VE = coverage * ve_direct + (1-coverage)*0 (if no indirect)
    # But usually "Overall Effectiveness" includes indirect effects.

    # If we follow the blog post logic/standard Halloran definitions:
    # Effectiveness is 1 - (Incidence_vax_pop / Incidence_unvax_pop)

    # Here we want Total Cases Reduction.
    # Let's define a simple user-provided "Overall Effectiveness" (OE)
    # OE = f(coverage, ve_direct)

    # Simple linear model (no herd immunity threshold non-linearity):
    # OE = coverage * ve_direct

    # If we have indirect effect:
    # Indirect protection for unvac: I_unvac_reduc = f(coverage)
    # Indirect protection for vac: usually assumed same or additional

    # Let's use the provided logic: TotalCases = TotalCases_0 * (1 - OE)

    # We can define OE using a standard "linear" assumption if func not provided
    oe <- coverage * ve_direct

    if (!is.null(ve_indirect_func)) {
        # If indirect function provided, it adds to the effect.
        # This is highly model dependent.
        # Let's assume the function returns the "Indirect Effectiveness"
        ie <- ve_indirect_func(coverage)
        # Combined? 1 - (1-Direct)*(1-Indirect)?
        # Or simply OE = DirectPart + IndirectPart?
        # For this exercise, let's assume the user supplies the full OE function or we use simple linear.
        oe <- ie # Assuming the function returns the full Overall Effectiveness
    }

    cases_post <- cases_baseline * (1 - oe)
    return(list(cases_post = cases_post, oe = oe))
}

simple_nonlinear_oe <- function(coverage, ve_direct, threshold = 0.8) {
    # Toy model: Linear up to threshold, then rapid increase (herd immunity)
    direct <- coverage * ve_direct
    indirect <- 0
    if (coverage > 0) {
        indirect <- (coverage^2) * 0.5 # Arbitrary non-linear boost
    }
    pmin(1, direct + indirect)
}

#' Calculate Indirect Vaccine Effectiveness (IVE) from ODE
#'
#' @param N Numeric, Total population size
#' @param I0 Numeric, Initial infected individuals
#' @param R0 Numeric, Basic reproduction number
#' @param gamma Numeric, Recovery rate
#' @param ve Numeric, Vaccine efficacy
#' @param rho Numeric, Vaccination rate
#' @param vacc_start Numeric, Start time of vaccination
#' @param vacc_duration Numeric, Duration of vaccination
#'
#' @return List with time series of IVE and final IVE value
#' @export
calculate_ive_ode <- function(N = 100000, I0 = 1, R0 = 2.0, gamma = 1 / 5,
                              ve = 0.7, rho = 0.05, vacc_start = 20, vacc_duration = 10) {
    # 1. Run Baseline (Unmitigated) Trend: No Vax
    sim_novax <- run_sir_vaccination(
        N = N, I0 = I0, R0 = R0, gamma = gamma,
        rho = 0, vacc_start = 0, vacc_duration = 0, ve = 0
    )

    # 2. Run Vaccination Scenario
    sim_vax <- run_sir_vaccination(
        N = N, I0 = I0, R0 = R0, gamma = gamma,
        rho = rho, vacc_start = vacc_start,
        vacc_duration = vacc_duration, ve = ve
    )

    # IVE(t) = 1 - (Incidence_unvax_in_vax_pop(t) / Incidence_unvax_in_novax_pop(t))
    # OR Cumulative risk version?
    # The blog post usually refers to attack rate reduction.
    # Let's use cumulative attack rate (Risk) to be consistent with "Overall Effectiveness" case calc.

    # Cumulative Risk Unvaccinated in NoVax Pop
    # AR_novax(t) = CumCases_novax(t) / N

    # Cumulative Risk Unvaccinated in Vax Pop
    # We need to track cases among UNVACCINATED separately in the Vax scenario.
    # Our run_sir_vaccination tracks Total CumCases.
    # We need to modify it or approximate.
    # Wait, run_sir_vaccination only returns Total CumCases. It doesn't split Cases by status (S->I vs V->I).
    # However, with a leaky vaccine (ve applied to S->I rate?), V individuals get infected at rate (1-ve)*lambda.
    # If our model moves S -> V, then V people have reduced susceptibility.
    # Our model:
    # new_inf = lambda * S
    # new_vac = rho * S
    # V compartment accumulates vaccinated.
    # Do V get infected?
    # Looking at run_sir_vaccination code:
    # S[t] <- S[t-1] - new_inf - new_vac
    # I[t] <- I[t-1] + new_inf - new_rec
    # It seems V compartment is "Removed" from S, but DOES NOT contribute to I?
    # If V individuals don't get infected, it's a "Perfect" vaccine (VE=1) for those who take it?
    # Wait, the code says:
    # new_vac <- rho * S[t-1] * time_step * ve
    # And V[t] <- V[t-1] + new_vac
    # This implies that `new_vac` amount of people are effectively "immune" and removed from S.
    # The remaining S are still susceptible.
    # People who are "vaccinated but failed" (1-ve) stay in S?
    # If new_vac only moves "protected" people, then "Total Population" in S+I+R+V = N.
    # Then S people are the "Unvaccinated" (plus those who failed vax).
    # The incidence in S is `new_inf`.
    # The "Unvaccinated Population" is S(t).
    # But in the NoVax scenario, "Unvaccinated Population" is S_novax(t).

    # Actually, the user's previous request was for "Overall Effectiveness approach" where Cases_post = Cases_baseline * (1 - OE).
    # And this request asks to compute IVE.
    # If the model is "All-or-Nothing" (which the current implementation suggests: only `ve` fraction move to V),
    # then those remaining in S are the "Unvaccinated" (including failures).
    # The risk in the "Unvaccinated" group in Vax scenario vs NoVax scenario.

    # Let's refine the model or interpretation.
    # Current model `run_sir_vaccination`:
    # new_vac = rho * S * ve
    # S -> V (immune)
    # So effectively, `ve` % of vaccinated become immune, (1-ve)% stay in S.
    # This is "All-or-Nothing".
    # In this case, comparing Cumulative Risk in S_vax vs S_novax measures indirect protection.

    # Risk_unvax_vax(t) = Cumulative_Infections_in_S_vax(t) / Initial_S
    # Wait, the denominator changes.
    # Force of Infection lambda(t) is lower in Vax scenario because I(t) is lower.
    # IVE = 1 - (lambda_vax(t) / lambda_novax(t))?
    # Or 1 - (CumCases_S_vax / N) / (CumCases_S_novax / N)?

    # Let's use the Cumulative Cases reduction for the whole population as "Overall Effectiveness".
    # Indirect Effectiveness (IVE) usually refers to 1 - (Risk_unprotected_vax / Risk_unprotected_novax).

    # In our model, S are the unprotected.
    # Risk_unprotected_vax(t) can be approximated by Int(lambda_vax) ?
    # Probability of infection for a susceptible: 1 - exp(- int lambda dt).

    # Let's calculate IVE based on 1 - (Cumulative Hazard Vax / Cumulative Hazard NoVax).
    # CumHazard(t) = Sum(beta * I(t) / N * dt)

    times <- sim_novax$time
    dt <- times[2] - times[1]
    beta <- R0 * gamma / N

    lambda_novax <- beta * sim_novax$I
    lambda_vax <- beta * sim_vax$I

    cum_hazard_novax <- cumsum(lambda_novax) * dt
    cum_hazard_vax <- cumsum(lambda_vax) * dt

    # Prob of infection for a purely susceptible person
    risk_novax <- 1 - exp(-cum_hazard_novax)
    risk_vax <- 1 - exp(-cum_hazard_vax)

    # IVE = 1 - (Risk_vax / Risk_novax)
    ive_series <- 1 - (risk_vax / risk_novax)

    # Handle 0/0 at t=0
    ive_series[is.nan(ive_series)] <- 0

    final_ive <- tail(ive_series, 1)

    list(
        time = times, ive_series = ive_series, final_ive = final_ive,
        risk_novax = risk_novax, risk_vax = risk_vax
    )
}
