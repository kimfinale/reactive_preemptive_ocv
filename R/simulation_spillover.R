#' Simulate outbreaks with single-generation spatial spillover
#'
#' This function simulates outbreak occurrence at the grid-cell level using a
#' two-stage process: (i) primary outbreaks driven by intrinsic risk, and
#' (ii) secondary outbreaks caused by spillover from neighboring primary outbreaks.
#'
#' Pre-emptive vaccination operates at the cell level and can reduce both
#' intrinsic outbreak risk and susceptibility to spillover.
#'
#' @param p Numeric vector of intrinsic (baseline) outbreak probabilities for
#'   each grid cell. Values must lie in [0,1].
#' @param adj List of integer vectors giving neighboring cells for each cell.
#' @param vax Optional logical vector indicating whether each cell was
#'   pre-emptively vaccinated. Defaults to all FALSE.
#' @param eta Spillover strength parameter: probability that a single infected
#'   neighbor generates a secondary outbreak in a susceptible cell.
#' @param nu Vaccine efficacy against primary outbreaks (risk reduction factor).
#' @param kappa Vaccine efficacy against spillover-induced outbreaks.
#'
#' @return A list with components:
#' \describe{
#'   \item{Y}{Logical vector indicating whether each cell experienced an outbreak.}
#'   \item{Y0}{Logical vector of primary (intrinsic) outbreaks.}
#'   \item{Y1}{Logical vector of secondary (spillover) outbreaks.}
#'   \item{K}{Integer vector giving the number of neighboring primary outbreaks
#'            for each cell.}
#' }
#'
#' @details
#' Primary outbreaks are drawn as:
#' \deqn{ Y_i^{(0)} \sim \text{Bernoulli}(p_i (1 - \nu vax_i)) }
#'
#' For cells without primary outbreak, secondary outbreak probability is:
#' \deqn{
#'   h_i = 1 - (1 - \eta (1 - \kappa vax_i))^{K_i},
#' }
#' where \eqn{K_i} is the number of neighboring primary outbreaks.
#'
#' Only a single generation of spillover is modeled; secondary outbreaks do not
#' generate further transmission.
#'
#' @examples
#' # p_true, adj, and vax defined elsewhere
#' out <- draw_outbreak_adj(p_true, adj, vax = vax)
#'
#' @export
simulate_outbreak_spillover <- function(p,
                                        adj,
                                        vax = NULL,
                                        eta = 0.6,
                                        nu = 0.8,
                                        kappa = 0.8) {

  ## ---- input checks ----
  stopifnot(
    is.numeric(p),
    all(is.finite(p)),
    all(p >= 0 & p <= 1),
    is.list(adj),
    length(adj) == length(p),
    is.numeric(eta), eta >= 0, eta <= 1,
    is.numeric(nu), nu >= 0, nu <= 1,
    is.numeric(kappa), kappa >= 0, kappa <= 1
  )

  n <- length(p)

  if (is.null(vax)) {
    vax <- rep(FALSE, n)
  } else {
    stopifnot(is.logical(vax), length(vax) == n)
  }

  ## ---- primary outbreaks (intrinsic risk) ----
  p_eff <- p * (1 - nu * vax)
  Y0 <- rbinom(n, size = 1, prob = p_eff) == 1L

  ## ---- secondary outbreaks (spillover from neighbors) ----
  K <- vapply(seq_len(n), function(i) sum(Y0[adj[[i]]]), integer(1))

  eta_eff <- eta * (1 - kappa * vax)
  h <- 1 - (1 - eta_eff)^K

  Y1 <- (!Y0) & (rbinom(n, size = 1, prob = h) == 1L)

  ## ---- final outbreak indicator ----
  Y <- Y0 | Y1

  list(
    Y  = Y,
    Y0 = Y0,
    Y1 = Y1,
    K  = K
  )
}
