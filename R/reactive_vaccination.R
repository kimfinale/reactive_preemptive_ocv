#' Assign reactive vaccination parameters to ongoing outbreaks
#'
#' This function assigns reactive vaccination (RV) to a subset of ongoing
#' outbreaks based on their epidemic stage and severity. It operates on the
#' output of \code{assign_outbreak_trajectory()} and augments the data frame
#' with indicators and parameters describing reactive vaccination deployment.
#'
#' Reactive vaccination is considered only for outbreaks that are ongoing at
#' the current time \code{t_now} and that are either still on the rising phase
#' or have experienced fewer than three weeks of consecutive decline after
#' the peak incidence.
#'
#' Eligible outbreaks are prioritized using a continuous priority score based
#' on a weighted combination of total attack rate and cumulative case burden.
#' Reactive vaccination is then assigned to the highest-priority outbreaks,
#' subject to a maximum number of targets.
#'
#' @param outbreaks Data frame produced by \code{assign_outbreak_trajectory()},
#'   containing outbreak timing and severity parameters. Must include
#'   \code{t0}, \code{tp}, \code{t1}, \code{attack_rate_per_capita}, and
#'   \code{pop_size}.
#'
#' @param t_now Numeric scalar giving the current time (in the same units as
#'   \code{t0}, \code{tp}, and \code{t1}) at which reactive vaccination decisions
#'   are made.
#'
#' @param max_targets Maximum number of outbreaks that can be targeted for
#'   reactive vaccination. Defaults to \code{Inf}, corresponding to no capacity
#'   constraint.
#'
#' @param coverage Numeric scalar giving the effective vaccination coverage
#'   achieved in each targeted outbreak. Defaults to 0.7.
#'
#' @param delay_range Numeric vector of length 2 giving the minimum and maximum
#'   delay (in days) between the decision time \code{t_now} and the start of
#'   reactive vaccination. Delays are sampled uniformly within this range for
#'   targeted outbreaks.
#'
#' @param w_attack Weight assigned to the normalized per-capita attack rate when
#'   constructing the priority score.
#'
#' @param w_cases Weight assigned to the normalized cumulative case burden when
#'   constructing the priority score.
#'
#' @return A data frame equal to \code{outbreaks} with additional columns:
#'   \describe{
#'     \item{has_outbreak}{Logical indicator of whether an outbreak occurred.}
#'     \item{ongoing}{Logical indicator of whether the outbreak is ongoing at
#'       time \code{t_now}.}
#'     \item{eligible}{Logical indicator of whether the outbreak is eligible for
#'       reactive vaccination.}
#'     \item{priority_score}{Continuous priority score used for targeting.}
#'     \item{reactive}{Logical indicator of whether reactive vaccination is
#'       assigned.}
#'     \item{rv_coverage}{Effective vaccination coverage for reactive vaccination
#'       (0 if not assigned).}
#'     \item{rv_delay}{Delay (in days) between \code{t_now} and the start of
#'       reactive vaccination (NA if not assigned).}
#'   }
#'
#' @details
#' An outbreak is considered ongoing if \code{t0 <= t_now <= t1}. Among ongoing
#' outbreaks, eligibility for reactive vaccination requires either:
#' \itemize{
#'   \item incidence is still increasing (\code{t_now < tp}), or
#'   \item incidence is declining but has done so for fewer than 21 days.
#' }
#'
#' The priority score is computed as a weighted sum of normalized severity
#' measures:
#' \deqn{
#' \text{priority}_i =
#' w_{\text{attack}} \tilde{A}_i +
#' w_{\text{cases}} \tilde{C}_i,
#' }
#' where \eqn{\tilde{A}_i} is the normalized attack rate and \eqn{\tilde{C}_i} is
#' the normalized cumulative case burden in outbreak \eqn{i}.
#'
#' This formulation is intended as a simple and interpretable default; users may
#' replace or extend it to incorporate mortality, age structure, or other
#' severity metrics.
#'
#' @examples
#' \dontrun{
#' outbreaks <- assign_outbreak_trajectory(regions)
#'
#' outbreaks <- assign_reactive_vaccination(
#'   outbreaks,
#'   t_now = 90,
#'   max_targets = 5,
#'   coverage = 0.8,
#'   delay_range = c(7, 14)
#' )
#' }
#'
#' @export
assign_reactive_vaccination <- function(outbreaks,
                                        t_now,
                                        max_targets = Inf,
                                        coverage = 0.7,
                                        delay_range = c(7, 21),
                                        w_attack = 0.6,
                                        w_cases = 0.4) {

  stopifnot(
    is.data.frame(outbreaks),
    is.numeric(t_now),
    is.numeric(delay_range), length(delay_range) == 2
  )

  df <- outbreaks

  ## ---- identify outbreaks ----
  df <- df |>
    dplyr::mutate(
      t_now = t_now,
      has_outbreak = !is.na(t0)
    )

  ## ---- outbreak status at current time ----
  df <- df |>
    dplyr::mutate(
      ongoing = has_outbreak & (t_now >= t0) & (t_now <= t1),
      time_since_peak = ifelse(has_outbreak, t_now - tp, NA_real_),
      eligible = ongoing & (is.na(time_since_peak) | time_since_peak < 21)
    )

  ## ---- severity proxies ----
  df <- df |>
    dplyr::mutate(
      cum_cases = attack_rate_per_capita * pop_size
    )

  ## ---- normalization helper ----
  norm <- function(x) {
    if (all(is.na(x))) return(rep(0, length(x)))
    (x - min(x, na.rm = TRUE)) /
      (max(x, na.rm = TRUE) - min(x, na.rm = TRUE) + 1e-9)
  }

  ## ---- priority score ----
  df <- df |>
    dplyr::mutate(
      attack_norm = norm(attack_rate_per_capita),
      cases_norm  = norm(cum_cases),
      priority_score =
        eligible * (w_attack * attack_norm + w_cases * cases_norm)
    )

  ## ---- select outbreaks ----
  df <- df |>
    dplyr::arrange(dplyr::desc(priority_score)) |>
    dplyr::mutate(
      reactive_vax = eligible &
        (dplyr::row_number() <= max_targets) &
        (priority_score > 0)
    )

  ## ---- assign coverage and delay ----
  df <- df |>
    dplyr::mutate(
      rv_coverage = ifelse(reactive_vax, coverage, 0),
      rv_delay = ifelse(
        reactive_vax,
        runif(sum(reactive_vax), delay_range[1], delay_range[2]),
        NA_real_
      )
    )

  df
}


#' Compute fractional reduction in cumulative incidence from reactive vaccination
#'
#' This function computes the fractional reduction in outbreak size (cumulative
#' incidence) attributable to reactive vaccination. The reduction depends on the
#' timing of vaccine effectiveness (delay), vaccine coverage, and vaccine efficacy.
#'
#' The function is designed to work directly with the output of
#' \code{assign_reactive_vaccination()}, and assumes that outbreak trajectories
#' follow a tent-shaped incidence curve defined by start time \code{t0}, peak time
#' \code{tp}, and end time \code{t1}.
#'
#' Vaccine protection is assumed to begin at time \code{t_eff = t0 + delay}. The
#' fraction of cumulative incidence remaining after \code{t_eff} is computed
#' analytically under the tent-curve assumption, and multiplied by coverage and
#' vaccine efficacy to obtain the total fractional reduction.
#'
#' @param df A data frame produced by \code{assign_reactive_vaccination()}.
#'   Must contain columns:
#'   \describe{
#'     \item{t0}{Outbreak start time.}
#'     \item{tp}{Outbreak peak time.}
#'     \item{t1}{Outbreak end time.}
#'     \item{reactive_vax}{Logical indicator for reactive vaccination.}
#'     \item{coverage}{Vaccine coverage (0--1).}
#'     \item{delay}{Delay until vaccine becomes effective (same time units as t0/t1).}
#'   }
#' @param ve Vaccine efficacy against infection (scalar in [0,1]).
#'
#' @return The input data frame with an additional column:
#' \describe{
#'   \item{frac_reduction_ci}{Fractional reduction in cumulative incidence due to
#'     reactive vaccination.}
#' }
#'
#' @details
#' Let \eqn{A_{\text{tot}}} denote the total area under the tent-shaped incidence
#' curve, and \eqn{A_{\text{rem}}(t)} the area remaining after time \eqn{t}. The
#' fractional reduction is computed as:
#'
#' \deqn{
#' \text{frac\_reduction} =
#' \mathbb{1}(\text{reactive\_vax})
#' \times \text{coverage}
#' \times \nu
#' \times \frac{A_{\text{rem}}(t_{\text{eff}})}{A_{\text{tot}}}
#' }
#'
#' where \eqn{t_{\text{eff}} = t_0 + \text{delay}}.
#'
#' If vaccine effectiveness begins after \code{t1}, the reduction is zero.
#'
#' @examples
#' df <- compute_frac_reduction_ci(df, nu = 0.65)
#'
#' @export
compute_frac_reduction_ci <- function(df, nu) {

  # helper: fraction of remaining area under a tent curve after time t
  frac_remaining_tent <- function(t, t0, tp, t1) {

    if (t >= t1) return(0)
    if (t <= t0) return(1)

    total_area <- 0.5 * (t1 - t0)

    if (t <= tp) {
      # remaining area = area of right triangle + area of descending part
      left_area <- 0.5 * (t - t0)^2 / (tp - t0)
      rem_area <- total_area - left_area
    } else {
      # remaining area is a smaller right triangle
      rem_area <- 0.5 * (t1 - t)^2 / (t1 - tp)
    }

    rem_area / total_area
  }

  t_eff <- df$t0 + df$rv_delay

  frac_remain <- mapply(
    frac_remaining_tent,
    t  = t_eff,
    t0 = df$t0,
    tp = df$tp,
    t1 = df$t1
  )

  df$frac_reduction_ci <-
    ifelse(
      df$reactive_vax,
      df$rv_coverage * nu * frac_remain,
      0
    )

  df
}
