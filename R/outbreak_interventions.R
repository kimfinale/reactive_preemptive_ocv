#' Apply Outbreak Response Immunization (ORI)
#'
#' @param incidence_curve Vector of daily incidence (No Vax).
#' @param delay Days after start until vaccination becomes effective.
#' @param coverage Proportion of the population effectively vaccinated (or reduction factor).
#' @param efficacy Vaccine efficacy (default 1 for simplicity if coverage handles effective protection).
#'
#' @return List with 'cases_averted', 'fraction_averted', and 'curve_vax'.
#' @export
apply_ori <- function(incidence_curve, delay, coverage, efficacy = 1.0) {
    # Create a reduction vector
    # Before delay: reduction = 0
    # After delay: reduction = coverage * efficacy

    n_days <- length(incidence_curve)
    reduction <- rep(0, n_days)

    if (delay < n_days) {
        reduction[(delay + 1):n_days] <- coverage * efficacy
    }

    curve_vax <- incidence_curve * (1 - reduction)

    total_cases_no_vax <- sum(incidence_curve)
    total_cases_vax <- sum(curve_vax)
    cases_averted <- total_cases_no_vax - total_cases_vax
    fraction_averted <- if (total_cases_no_vax > 0) cases_averted / total_cases_no_vax else 0

    return(list(
        cases_averted = cases_averted,
        fraction_averted = fraction_averted,
        curve_vax = curve_vax
    ))
}
