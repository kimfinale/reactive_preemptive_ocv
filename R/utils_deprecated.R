#' Reactive vs. Preemptive Vaccination Simulation Functions
#'
#' This script contains functions to simulate outbreaks across subpopulations and evaluate
#' the relative impact of reactive and preemptive vaccination strategies. Outbreaks follow
#' tent-shaped incidence curves. Vaccination modifies the incidence based on coverage, efficacy,
#' and timing of response.
#'
#' Author: Jong-Hoon Kim
#' Date: 2025-06-27



#' Set up a spatial region with intrinsic outbreak risk and patch structure
#'
#' This function constructs a spatial region composed of grid cells, assigns
#' spatially correlated intrinsic outbreak risks, derives prioritization scores,
#' builds contiguous patches, and aggregates outputs into clean cell-level and
#' patch-level containers.
#'
#' A single random seed (\code{seed_rng}) controls all stochastic components,
#' including spatial risk generation, imperfect ranking (if used), and any
#' stochasticity inside patch construction. This ensures full reproducibility
#' and safe use inside Monte Carlo simulations.
#'
#' @param n_region Integer. Number of spatial cells.
#'
#' @param theta Positive numeric. Shape-2 parameter for the
#'   \eqn{\mathrm{Beta}(1,\theta)} marginal distribution of intrinsic outbreak
#'   risk.
#'
#' @param rho0 Numeric in \eqn{[0,1)}. Strength of spatial dependence in the
#'   Gaussian copula used to generate correlated risks.
#'
#' @param ell Positive numeric. Spatial range parameter controlling correlation
#'   decay with distance.
#'
#' @param adj_type Character. Neighborhood definition passed to
#'   \code{make_adj()}. One of \code{"distance"}, \code{"rook"}, or
#'   \code{"queen"}.
#'
#' @param d0 Optional numeric. Distance threshold used only when
#'   \code{adj_type = "distance"}.
#'
#' @param rho_score Numeric in \eqn{(-1,1)}. Target Spearman correlation between
#'   intrinsic risk and prioritization score. Used only when
#'   \code{score_mode = "ranking"}.
#'
#' @param n_seeds Integer. Number of seed cells from which patches are grown.
#'
#' @param patch_size Integer. Maximum number of cells per patch.
#'
#' @param seed_rng Optional integer. Single random seed controlling all
#'   stochastic components inside this function and all sub-functions.
#'
#' @param ... Additional arguments forwarded to \code{make_patches()}.
#'
#' @return A named list with three components:
#' \describe{
#'   \item{cells}{Cell-level container with elements:
#'     \describe{
#'       \item{id}{Integer vector of cell indices.}
#'       \item{coords}{\eqn{n \times 2} coordinate matrix.}
#'       \item{adj}{Adjacency list.}
#'       \item{risk}{Output of \code{spatial_beta_risk()}, including intrinsic
#'       outbreak risks \code{p}.}
#'       \item{score}{Prioritization score used for patch construction.}
#'       \item{rank}{Ranking object returned by \code{make_ranking()} (or
#'       \code{NULL} if \code{score_mode = "p_true"}).}
#'       \item{patch_id}{Integer vector mapping each cell to a patch
#'       (0 = unassigned).}
#'     }}
#'   \item{patches}{Patch-level container with elements:
#'     \describe{
#'       \item{patches}{List of integer vectors, each giving cell indices in a
#'       patch.}
#'       \item{summary}{Patch-level summary data frame produced by
#'       \code{summarize_patches()}, suitable for pre-emptive allocation.}
#'       \item{patch_out}{Full output returned by \code{make_patches()}.}
#'     }}
#'   \item{params}{List of key parameters used to generate the object.}
#' }
#'
#' @details
#' The returned object is designed to separate **cell-level state** from
#' **patch-level structure**, facilitating flexible comparisons of
#' cell-based versus patch-based vaccination strategies.
#'
#' Intrinsic outbreak risks are generated using a spatial Gaussian copula with
#' \eqn{\mathrm{Beta}(1,\theta)} marginals. Patches are grown deterministically
#' from the supplied score and adjacency structure, conditional on
#' \code{seed_rng}.
#'
#' @examples
#' reg <- setup_spatial_region(
#'   n_region   = 400,
#'   theta      = 3,
#'   rho0       = 0.6,
#'   ell        = 2,
#'   adj_type   = "rook",
#'   rho_score  = 0.7,
#'   n_seeds    = 10,
#'   patch_size = 40,
#'   seed_rng   = 1
#' )
#'
#' reg$patches$summary
#'
#' @export
setup_spatial_region <- function(
    n_region = 1000,
    pop_size_range = c(1000, 1000),
    prop_U5_range = c(0.15, 0.18),
    theta = 3,
    rho0 = 0.7,
    ell = 1.0,
    adj_type = c("distance", "rook", "queen"),
    d0 = NULL,
    rho_score = 0.7,
    n_seeds = 10,
    patch_size = 40,
    seed_rng = NULL,...) {

  adj_type    <- match.arg(adj_type)

  if (!is.null(seed_rng)) {
    set.seed(seed_rng)
  }

  ## ---- coordinates ----
  nx <- ceiling(sqrt(n_region))
  ny <- ceiling(n_region / nx)
  grd <- expand.grid(x = seq_len(nx), y = seq_len(ny))
  coords <- as.matrix(grd[seq_len(n_region), 1:2])
  colnames(coords) <- c("x", "y")

  ## ---- adjacency ----
  adj <- make_adj(coords, type = adj_type, d0 = d0)

  ## ---- intrinsic spatial risk ----
  risk <- spatial_beta_risk(
    coords   = coords,
    theta    = theta,
    rho0     = rho0,
    ell      = ell,
    seed_rng = seed_rng
  )

  p_true <- risk$p

  ## ---- score for prioritization ----

  rank_obj <- make_ranking(
    p          = p_true,
    rho_target = rho_score,
    seed_rng   = seed_rng
  )
  score <- rank_obj$score

  ## ---- patches ----
  patch_out <- make_patches(
    score,
    adj,
    n_seeds    = n_seeds,
    patch_size = patch_size,
    seed_rng   = seed_rng,
    ...
  )

  ## ---- cell-level container ----
  cells <- list(
    id       = seq_len(n_region),
    pop_size = runif(n_region, pop_size_range[1], pop_size_range[2]),
    prop_U5  = runif(n_region, prop_U5_range[1], prop_U5_range[2]),
    coords   = coords,
    adj      = adj,
    risk     = risk,
    score    = score,
    rank     = rank_obj,
    patch_id = patch_out$patch_id
  )

  ## ---- patch-level container ----
  patches <- list(
    patches   = patch_out$patches,
    summary   = summarize_patches(patch_obj = patch_out),
    patch_out = patch_out
  )

  ## ---- return ----
  list(
    cells   = cells,
    patches = patches,
    params  = list(
      n_region   = n_region,
      theta      = theta,
      rho0       = rho0,
      ell        = ell,
      adj_type   = adj_type,
      d0         = d0,
      rho_score  = rho_score,
      n_seeds    = n_seeds,
      patch_size = patch_size,
      seed_rng   = seed_rng
    )
  )
}


#' Simulate spatially correlated outbreak risks with Beta marginals
#'
#' This function generates spatially correlated baseline outbreak probabilities
#' using a Gaussian copula. A latent Gaussian field with exponential
#' distance-decay correlation is transformed to Uniform(0,1) via the standard
#' normal CDF and then mapped to Beta(1, theta) marginals.
#'
#' Spatial dependence is controlled by two interpretable parameters:
#' \describe{
#'   \item{rho0}{Correlation strength (off-diagonal scaling).}
#'   \item{ell}{Correlation range (larger values imply longer-range dependence).}
#' }
#'
#' The resulting outbreak risks are suitable for downstream simulation of
#' spatially clustered outbreaks and spillover effects.
#'
#' @param coords Numeric matrix or data.frame with at least two columns giving
#'   spatial coordinates (e.g., x and y) for each cell.
#' @param theta Positive shape parameter of the Beta(1, theta) marginal
#'   distribution. Larger values imply lower mean outbreak risk.
#' @param rho0 Numeric in [0, 1). Controls the overall strength of spatial
#'   correlation between distinct locations.
#' @param ell Positive numeric. Range parameter for exponential decay of
#'   correlation with distance.
#' @param seed Optional integer seed for reproducibility. If NULL, the RNG state
#'   is not modified.
#'
#' @return A list with components:
#' \describe{
#'   \item{p}{Vector of spatially correlated outbreak probabilities.}
#'   \item{Z}{Latent Gaussian field.}
#'   \item{Sigma}{Spatial correlation matrix used to generate \code{Z}.}
#'   \item{D}{Pairwise distance matrix derived from \code{coords}.}
#' }
#'
#' @details
#' The spatial correlation matrix is defined as
#' \deqn{
#'   \Sigma_{ij} =
#'   \begin{cases}
#'     1, & i = j \\
#'     \rho_0 \exp(-D_{ij} / \ell), & i \neq j
#'   \end{cases}
#' }
#' where \eqn{D_{ij}} is the Euclidean distance between locations \eqn{i} and
#' \eqn{j}. A small diagonal jitter is added for numerical stability.
#'
#' @seealso \code{\link[MASS]{mvrnorm}}
#'
#' @examples
#' coords <- expand.grid(x = 1:10, y = 1:10)
#' sim <- sim_spatial_beta_risk(coords, theta = 3, rho0 = 0.6, ell = 2)
#' hist(sim$p, main = "Simulated outbreak risks", xlab = "p")
#'
#' @export
spatial_beta_risk <- function(coords,
                              theta,
                              rho0,
                              ell,
                              seed_rng = NULL) {

  ## ---- input checks ----
  coords <- as.matrix(coords)
  stopifnot(
    is.numeric(coords),
    ncol(coords) >= 2,
    is.numeric(theta), length(theta) == 1, theta > 0,
    is.numeric(rho0), length(rho0) == 1, rho0 >= 0, rho0 < 1,
    is.numeric(ell),  length(ell)  == 1, ell > 0
  )

  n <- nrow(coords)

  ## ---- distances ----
  D <- as.matrix(dist(coords))

  ## ---- spatial correlation matrix ----
  Sigma <- rho0 * exp(-D / ell)
  diag(Sigma) <- 1
  Sigma <- Sigma + diag(1e-10, n)  # numerical jitter

  ## ---- latent Gaussian field ----
  if (!is.null(seed_rng)) set.seed(seed_rng)
  Z <- MASS::mvrnorm(n = 1, mu = rep(0, n), Sigma = Sigma)

  ## ---- Gaussian copula to Beta ----
  U <- pnorm(Z)
  p <- qbeta(U, shape1 = 1, shape2 = theta)

  list(
    p = p,
    Z = Z,
    Sigma = Sigma,
    D = D
  )
}


#' Construct adjacency lists for spatial units
#'
#' This function creates neighborhood (adjacency) lists for spatial units
#' defined by coordinates. It supports distance-based neighborhoods for
#' continuous space and rook/queen contiguity for regular grids.
#'
#' @param coords Numeric matrix or data.frame with at least two columns giving
#'   spatial coordinates (e.g., x and y) for each spatial unit.
#' @param type Character string specifying the neighborhood definition.
#'   One of \code{"distance"}, \code{"rook"}, or \code{"queen"}.
#' @param d0 Numeric distance threshold used only when
#'   \code{type = "distance"}. Two units are neighbors if their distance is
#'   greater than 0 and less than or equal to \code{d0}.
#'
#' @return A list of length \code{nrow(coords)}. The \code{i}-th element contains
#'   the indices of neighboring units of unit \code{i}.
#'
#' @details
#' Neighborhood definitions:
#' \describe{
#'   \item{distance}{
#'     Neighbors are all units within Euclidean distance \code{d0}.
#'     Suitable for continuous spatial layouts.
#'   }
#'   \item{rook}{
#'     Neighbors share an edge (Manhattan distance = 1).
#'     Assumes coordinates lie on a regular integer grid.
#'   }
#'   \item{queen}{
#'     Neighbors share an edge or a corner (Chebyshev distance = 1).
#'     Assumes coordinates lie on a regular integer grid.
#'   }
#' }
#'
#' Self-neighbors are excluded in all cases.
#'
#' @examples
#' # Regular grid
#' coords <- expand.grid(x = 1:5, y = 1:5)
#' adj_rook  <- make_adj(coords, type = "rook")
#' adj_queen <- make_adj(coords, type = "queen")
#'
#' # Distance-based
#' set.seed(1)
#' coords2 <- cbind(runif(20), runif(20))
#' adj_dist <- make_adj(coords2, type = "distance", d0 = 0.3)
#'
#' @export
make_adj <- function(coords,
                     type = c("distance", "rook", "queen"),
                     d0 = NULL) {

  ## ---- input checks ----
  coords <- as.matrix(coords)
  stopifnot(
    is.numeric(coords),
    ncol(coords) >= 2
  )

  type <- match.arg(type)
  n <- nrow(coords)

  adj <- vector("list", n)

  ## ---- distance-based adjacency ----
  if (type == "distance") {

    if (is.null(d0) || !is.numeric(d0) || length(d0) != 1 || d0 <= 0)
      stop("For type = 'distance', d0 must be a positive numeric scalar.")

    D <- as.matrix(dist(coords))

    for (i in seq_len(n)) {
      adj[[i]] <- which(D[i, ] > 0 & D[i, ] <= d0)
    }

  } else {

    ## ---- grid-based adjacency (rook / queen) ----
    x <- coords[, 1]
    y <- coords[, 2]

    for (i in seq_len(n)) {

      dx <- abs(x - x[i])
      dy <- abs(y - y[i])

      if (type == "rook") {
        adj[[i]] <- which(dx + dy == 1)
      } else {  # queen
        adj[[i]] <- which(pmax(dx, dy) == 1 & (dx + dy > 0))
      }
    }
  }

  adj
}



#' Compute expected outbreak risk augmented by spatial spillover
#'
#' This function computes the expected (marginal) probability that each spatial
#' unit experiences an outbreak, combining its baseline risk with the risk of
#' secondary outbreaks arising from neighboring units via spillover.
#'
#' Spillover is treated probabilistically: the function integrates over the
#' (unobserved) outbreak states of neighboring units rather than conditioning
#' on realized outbreaks.
#'
#' @param prob Numeric vector of baseline (primary) outbreak probabilities for
#'   each spatial unit.
#' @param adj List of integer vectors giving neighbors for each unit (as produced
#'   by \code{make_adj}).
#' @param eta Numeric scalar in [0,1]. Spillover strength: probability that an
#'   outbreak in a neighboring unit generates a secondary outbreak.
#'
#' @return Numeric vector of expected outbreak probabilities, augmented by
#'   spatial spillover and constrained to [0,1].
#'
#' @details
#' Let \eqn{p_i} denote the baseline outbreak probability of unit \eqn{i}.
#' The augmented risk is computed as
#' \deqn{
#'   p_i^{\mathrm{aug}}
#'   =
#'   p_i
#'   +
#'   (1 - p_i)
#'   \left[
#'     1 - \prod_{j \in \mathcal{N}(i)} (1 - \eta p_j)
#'   \right].
#' }
#'
#' This corresponds to the marginal probability that unit \eqn{i} experiences
#' either a primary outbreak or a secondary outbreak due to spillover from at
#' least one neighbor, under the assumptions that:
#' \itemize{
#'   \item primary outbreaks occur independently across units,
#'   \item spillover from each neighbor occurs independently,
#'   \item secondary outbreaks do not generate further spillover.
#' }
#'
#' @examples
#' p <- c(0.1, 0.2, 0.05)
#' adj <- list(c(2), c(1,3), c(2))
#' expected_risk_spillover(p, adj, eta = 0.3)
#'
#' @export
augment_risk_spillover <- function(prob, adj, eta) {

  ## ---- input checks ----
  stopifnot(
    is.numeric(prob),
    all(is.finite(prob)),
    all(prob >= 0 & prob <= 1),
    is.list(adj),
    length(adj) == length(prob),
    is.numeric(eta),
    length(eta) == 1,
    eta >= 0,
    eta <= 1
  )

  n <- length(prob)
  p_aug <- numeric(n)

  for (i in seq_len(n)) {
    nb <- adj[[i]]
    if (length(nb) == 0L) {
      p_aug[i] <- prob[i]
    } else {
      p_aug[i] <- prob[i] +
        (1 - prob[i]) *
        (1 - prod(1 - eta * prob[nb]))
    }
  }

  pmin(pmax(p_aug, 0), 1)
}



#' Construct an imperfect risk ranking with controlled Spearman correlation
#'
#' This function generates a noisy risk score whose ranking has a target
#' Spearman correlation with the true underlying risk. It is designed to model
#' imperfect risk assessment or prioritization in vaccination and outbreak
#' response strategies.
#'
#' The method uses a Gaussian copula construction:
#' \deqn{
#'   \text{score}_i = a z_i + \sqrt{1 - a^2}\,\varepsilon_i,
#' }
#' where \eqn{z_i} is a monotone transform of the true risk and \eqn{\varepsilon_i}
#' is fixed noise. The parameter \eqn{a} is chosen numerically so that the
#' resulting ranking achieves the target Spearman correlation.
#'
#' @param prob Numeric vector of true underlying risks.
#' @param rho_target Target Spearman rank correlation between \code{prob}
#'   and the generated score.
#' @param seed Optional integer seed. If NULL, the RNG state is not modified.
#' @param ties_method Method passed to \code{\link{rank}} for handling ties in
#'   \code{prob}.
#' @param tol Numeric tolerance for root finding.
#' @param maxiter Maximum number of iterations for root finding.
#'
#' @return A list with components:
#' \describe{
#'   \item{score}{Numeric vector of noisy scores used for ranking.}
#'   \item{ord}{Indices that sort \code{score} in decreasing order.}
#'   \item{rho_target}{Requested Spearman correlation.}
#'   \item{rho_achieved}{Achieved Spearman correlation.}
#'   \item{a_used}{Mixing parameter used in the Gaussian copula.}
#'   \item{seed}{Seed used (if any).}
#' }
#'
#' @details
#' Internally, Spearman correlation is computed using a monotone normal score
#' transform of \code{prob}. This stabilizes the mapping between the mixing
#' parameter \eqn{a} and the achieved rank correlation.
#'
#' When root-finding fails to bracket a solution (e.g., due to many ties),
#' the function falls back to minimizing absolute correlation error.
#'
#' @examples
#' set.seed(1)
#' p <- runif(50)
#' rk <- make_ranking(p, rho_target = 0.6)
#' cor(p, rk$score, method = "spearman")
#'
#' @export
make_ranking <- function(p,
                         rho_target = 0.7,
                         ties_method = "average",
                         tol = 1e-6,
                         maxiter = 200,
                         seed_rng = NULL) {

  ## ---- input checks ----
  stopifnot(
    is.numeric(p),
    all(is.finite(p)),
    length(p) >= 3,
    is.numeric(rho_target),
    length(rho_target) == 1
  )

  rho_target <- max(min(rho_target, 0.999999), -0.999999)
  n <- length(p)

  ## ---- monotone transform of true risk ----
  u <- rank(p, ties.method = ties_method) / (n + 1)
  z <- qnorm(u)

  ## ---- fixed noise for deterministic a -> rho mapping ----
  if (!is.null(seed_rng)) set.seed(seed_rng)
  eps <- rnorm(n)

  ## ---- objective: Spearman(z, score) - rho_target ----
  rho_of_a <- function(a) {
    a <- max(min(a, 0.999999), -0.999999)
    score <- a * z + sqrt(1 - a^2) * eps
    cor(z, score, method = "spearman")
  }

  f <- function(a) rho_of_a(a) - rho_target

  ## ---- root finding ----
  lo <- -0.9999
  hi <-  0.9999
  flo <- f(lo)
  fhi <- f(hi)

  if (!is.finite(flo) || !is.finite(fhi) || flo * fhi > 0) {
    opt <- optimize(function(a) abs(f(a)), interval = c(lo, hi))
    a_star <- opt$minimum
  } else {
    a_star <- uniroot(f, interval = c(lo, hi),
                      tol = tol, maxiter = maxiter)$root
  }

  ## ---- final score and ranking ----
  score <- a_star * z + sqrt(1 - a_star^2) * eps
  ord <- order(score, decreasing = TRUE)
  rho_achieved <- cor(p, score, method = "spearman")

  list(
    score = score,
    ord = ord,
    rho_target = rho_target,
    rho_achieved = rho_achieved,
    a_used = a_star,
    seed_rng = seed_rng
  )
}



#' Construct spatial patches using score-seeded greedy growth
#'
#' This function partitions spatial units into contiguous patches by first
#' generating an imperfect risk score (via \code{make_ranking()}) and then
#' selecting high-scoring seed units and greedily growing each patch through
#' neighboring units according to a specified growth rule.
#'
#' In addition to patch membership, the function returns the probability vector
#' used to create the score and the realized score itself. This allows
#' downstream functions (e.g., \code{summarize_patches()}) to work directly with
#' the output of \code{make_patches()}.
#'
#' @param p Numeric vector of baseline (cell-level) outbreak probabilities.
#'   Values should lie in \eqn{[0,1]}.
#' @param adj List of integer vectors giving neighbors for each unit.
#' @param rho_target Target Spearman rank correlation between \code{p} and
#'   the generated score used for patch seeding/growth.
#' @param n_seeds Integer number of patches (seeds) to create.
#' @param patch_size Integer maximum size of each patch.
#' @param grow_mode Character string specifying how new units are added to a
#'   growing patch:
#'   \describe{
#'     \item{"best"}{Add the neighboring unit with the highest score (default).}
#'     \item{"random"}{Add a random neighboring unit.}
#'   }
#' @param seed_rng Optional integer seed used for \code{make_ranking()}, and also
#'   used for stochastic growth when \code{grow_mode="random"}.
#'
#' @return A list with components:
#' \describe{
#'   \item{p}{The input probability vector used to generate the score.}
#'   \item{score}{The realized noisy score used for ranking/patch growth.}
#'   \item{rho_target}{Requested Spearman rank correlation used in scoring.}
#'   \item{patch_id}{Integer vector assigning each unit to a patch (0 = unassigned).}
#'   \item{patches}{List of integer vectors, each containing the units in a patch.}
#' }
#'
#' @details
#' The score is generated by \code{make_ranking(p, rho_target=...)} and is
#' used only for ordering (seed selection and growth). Patch construction is a
#' greedy heuristic and does not guarantee globally optimal patches.
#'
#' @examples
#' # score comes from make_ranking() internally
#' prob <- runif(100, 0, 0.1)
#' coords <- expand.grid(x = 1:10, y = 1:10)
#' adj <- make_adj(coords, type = "rook")
#' out <- make_patches(prob, adj, rho_target = 0.7, n_seeds = 5, patch_size = 20)
#' patch_df <- summarize_patches(out)
#'
#' @export
make_patches <- function(score,
                         adj,
                         n_seeds = 10,
                         patch_size = 25,
                         grow_mode = c("best", "random"),
                         seed_rng = NULL) {

  ## ---- input checks ----
  stopifnot(
    is.numeric(score),
    all(is.finite(score)),
    is.list(adj),
    length(adj) == length(score),
    is.numeric(n_seeds), length(n_seeds) == 1, n_seeds >= 1,
    is.numeric(patch_size), length(patch_size) == 1, patch_size >= 1
  )

  grow_mode <- match.arg(grow_mode)

  ## ---- if stochastic growth, set RNG (optional) ----
  if (grow_mode == "random" && !is.null(seed_rng)) {
    set.seed(seed_rng)
  }

  n <- length(score)
  assigned <- integer(n)  # 0 = unassigned
  order_score <- order(score, decreasing = TRUE)
  patches <- vector("list", n_seeds)

  for (k in seq_len(n_seeds)) {

    ## choose seed: highest-score unassigned
    seed <- NA_integer_
    for (i in order_score) {
      if (assigned[i] == 0L) { seed <- i; break }
    }
    if (is.na(seed)) break

    nodes <- seed
    assigned[seed] <- k
    frontier <- adj[[seed]]

    ## grow patch
    while (length(nodes) < patch_size && length(frontier) > 0L) {

      cand <- frontier[assigned[frontier] == 0L]
      if (length(cand) == 0L) break

      if (grow_mode == "best") {
        j <- cand[which.max(score[cand])]
      } else {
        j <- sample(cand, 1L)
      }

      nodes <- c(nodes, j)
      assigned[j] <- k
      frontier <- unique(c(frontier, adj[[j]]))
    }

    patches[[k]] <- nodes
  }

  list(
    score = score,
    patch_id = assigned,
    patches = patches
  )
}



#' Summarize patch-level quantities from a patch object
#'
#' This function computes patch-level summaries (size, patch outbreak probability,
#' and aggregated score). It is designed to work directly with the list returned
#' by \code{make_patches()}, which already contains \code{prob}, \code{score}, and
#' \code{patches}.
#'
#' For backward compatibility, \code{summarize_patches(prob, score_cell, patches)}
#' is also supported.
#'
#' @param patch_obj Either:
#'   \itemize{
#'     \item a list returned by \code{make_patches()} containing \code{prob},
#'     \code{score}, and \code{patches}, or
#'     \item a numeric vector \code{prob} (old calling style).
#'   }
#' @param score_cell (Old calling style only) Numeric vector of cell-level scores.
#' @param patches (Old calling style only) List of integer vectors of patch indices.
#'
#' @return A data.frame with one row per patch and columns:
#' \describe{
#'   \item{patch_id}{Patch identifier (1..K).}
#'   \item{m}{Number of cells in the patch.}
#'   \item{P_true}{Patch outbreak probability: \eqn{1-\prod(1-p_i)}.}
#'   \item{S_hat}{Mean score within the patch.}
#' }
#'
#' @details
#' Patch outbreak probability is computed as
#' \deqn{ 1 - \prod_{i\in patch}(1 - p_i), }
#' assuming conditional independence of cell-level outbreak events given \code{prob}.
#' The computation uses a numerically stable log form.
#'
#' @examples
#' # New style (recommended)
#' out <- make_patches(prob, adj, rho_target = 0.7)
#' patch_df <- summarize_patches(out)
#'
#' # Old style (still works)
#' patch_df2 <- summarize_patches(prob, out$score, out$patches)
#'
#' @export
summarize_patches <- function(patch_obj) {

  score_cell <- patch_obj$score
  patches <- patch_obj$patches

  ## ---- input checks ----
  stopifnot(
    is.numeric(score_cell),
    all(is.finite(score_cell)),
    is.list(patches)
  )

  n <- length(score_cell)
  K <- length(patches)

  m <- integer(K)
  P_true <- numeric(K)
  S_hat <- numeric(K)

  for (k in seq_len(K)) {
    idx <- patches[[k]]

    if (is.null(idx) || length(idx) == 0L) {
      stop(sprintf("Patch %d is empty or NULL.", k))
    }

    idx <- unique(as.integer(idx))
    if (any(idx < 1L | idx > n)) {
      stop(sprintf("Patch %d contains invalid cell indices.", k))
    }

    m[k] <- length(idx)

    ## stable: 1 - prod(1 - prob_cell[idx])
    # P_true[k] <- 1 - exp(sum(log1p(-prob_cell[idx])))

    ## aggregate score (mean)
    S_hat[k] <- mean(score_cell[idx])
  }

  data.frame(
    patch_id = seq_len(K),
    m = m,
    # P_true = P_true,
    S_hat = S_hat
  )
}


#' Allocate pre-emptive vaccination budget to top-ranked patches
#'
#' This function allocates a fixed pre-emptive vaccination budget across
#' spatial patches using a greedy, score-based targeting rule. Patches are
#' ranked by their prioritization score, and whole patches are selected
#' sequentially until the budget is exhausted.
#'
#' @param patch_df Data frame of patch-level summaries, as produced by
#'   \code{summarize_patches()}. Must contain columns:
#'   \describe{
#'     \item{m}{Patch size (number of cells).}
#'     \item{S_hat}{Patch prioritization score.}
#'   }
#' @param B_pre Non-negative integer giving the total number of cells that can
#'   be vaccinated pre-emptively.
#' @param patch_id Optional integer vector mapping each cell to a patch
#'   (0 = unassigned). If provided, a cell-level vaccination indicator is returned.
#'
#' @return A list with components:
#' \describe{
#'   \item{pre_vax}{Logical vector indicating which patches are selected for
#'     pre-emptive vaccination.}
#'   \item{used}{Total number of vaccine units used.}
#'   \item{V}{(Optional) Logical vector indicating which cells are vaccinated
#'     pre-emptively. Returned only if \code{patch_id} is provided.}
#' }
#'
#' @details
#' Patches are sorted in decreasing order of \code{S_hat}. Each patch is selected
#' if vaccinating all its cells does not exceed the remaining budget.
#'
#' This is a greedy heuristic and does not guarantee a globally optimal solution
#' to the knapsack problem. However, it is fast, interpretable, and consistent
#' with prioritization-based decision rules.
#'
#' @examples
#' patch_df <- data.frame(
#'   m = c(10, 20, 15),
#'   S_hat = c(0.8, 0.6, 0.9)
#' )
#'
#' # Patch-level only
#' allocate_preemptive_patches(patch_df, B_pre = 25)
#'
#' # With cell-level vaccination vector
#' patch_id <- c(1,1,2,2,3,3,3)
#' allocate_preemptive_patches(patch_df, B_pre = 25, patch_id = patch_id)
#'
#' @export
allocate_preemptive_patches <- function(patch_df, B_pre, patch_id = NULL) {

  ## ---- input checks ----
  stopifnot(
    is.data.frame(patch_df),
    all(c("m", "S_hat") %in% names(patch_df)),
    is.numeric(patch_df$m),
    is.numeric(patch_df$S_hat),
    is.numeric(B_pre),
    length(B_pre) == 1,
    B_pre >= 0
  )

  n_patch <- nrow(patch_df)

  ord <- order(patch_df$S_hat, decreasing = TRUE)
  pre_vax <- logical(n_patch)
  used <- 0

  for (k in ord) {

    if (used >= B_pre) break

    cost_k <- patch_df$m[k]

    if (used + cost_k <= B_pre) {
      pre_vax[k] <- TRUE
      used <- used + cost_k
    }
  }

  out <- list(
    pre_vax = pre_vax,
    used = used
  )

  ## ---- optional: cell-level vaccination indicator ----
  if (!is.null(patch_id)) {
    stopifnot(is.integer(patch_id) || is.numeric(patch_id))
    V <- patch_id > 0 & pre_vax[patch_id]
    out$V <- as.logical(V)
  }

  out
}

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














#' Compute instantaneous infection rate for tent-shaped outbreak
#'
#' Given a tent-shaped outbreak curve defined by start ($t_0$), peak ($t_p$), and end ($t_1$) times,
#' this function returns the infection rate at a specific time point $t$.
#'
#' @param t Time at which to compute the infection rate
#' @param t0 Time of outbreak start
#' @param tp Time of outbreak peak
#' @param t1 Time of outbreak end
#' @param r0 Baseline incidence rate (typically 0)
#' @param rp Peak incidence rate per capita
#'
#' @return Numeric value of the infection rate at time \code{t}
#' @examples
#' get_IR(t = 5, t0 = 0, tp = 10, t1 = 30, r0 = 0, rp = 0.02)
get_IR <- function(t, t0, tp, t1, r0, rp) {
  if (t <= tp) {
    r0 + (rp - r0) * (t - t0) / (tp - t0)
  } else {
    rp - (rp - r0) * (t - tp) / (t1 - tp)
  }
}

#' Compute area up to time t
#'
#' Integrates the incidence rate curve from \code{t0} to \code{t} assuming a triangular
#' (tent-shaped) outbreak with peak at \code{tp}. Uses piecewise linear incidence.
#'
#' @param t Time up to which area is calculated
#' @param t0 Time of outbreak start
#' @param tp Time of outbreak peak
#' @param t1 Time of outbreak end
#' @param r0 Baseline incidence rate
#' @param rp Peak incidence rate
#'
#' @return Area under incidence rate curve up to time t
#' @examples
#' get_area(t = 15, t0 = 0, tp = 10, t1 = 30, r0 = 0, rp = 0.02)
get_area <- function(t, t0, tp, t1, r0, rp) {
  if (t <= tp) {
    # Rising slope: triangle from t0 to t
    h <- r0 + (rp - r0) * (t - t0) / (tp - t0) - r0
    return(0.5 * (t - t0) * h)
  } else if (t <= t1) {
    # Rising triangle: t0 to tp
    A1 <- 0.5 * (tp - t0) * (rp - r0)

    # Falling trapezoid: tp to t
    r_t <- rp - (rp - r0) * (t - tp) / (t1 - tp)  # incidence at time t
    A2 <- 0.5 * (t - tp) * (rp + r_t - 2 * r0)    # subtract r0 for area above baseline

    return(A1 + A2)
  } else {
    # Total area under entire outbreak curve (from t0 to t1)
    A1 <- 0.5 * (tp - t0) * (rp - r0)
    A2 <- 0.5 * (t1 - tp) * (rp - r0)
    return(A1 + A2)
  }
}



#' #' Find the time at which cumulative incidence equals a given target number of cases
#' #'
#' #' @param C Target cumulative number of cases
#' #' @param N Population size
#' #' @param t0, tp, t1 Timepoints defining tent-shaped outbreak
#' #' @param r0, rp Baseline and peak incidence rate (per capita)
#' #'
#' #' @return Time t such that cumulative incidence equals C / N
#' get_time_for_CI <- function(C, N, t0, tp, t1, r0, rp) {
#'   area_target <- C / N # normalize with population size
#'
#'   obj_fun <- function(t) {
#'     get_area(t, t0, tp, t1, r0, rp) - area_target
#'   }
#'   uniroot(obj_fun, lower = t0, upper = t1)$root
#' }
#'
#' #' Generate a list of synthetic regions with population and outbreak probabilities
#' #'
#' #' @param n_region Number of regions to generate
#' #' @param pop_range Range of population sizes (min, max)
#' #' @param prop_U5_range Range of proportion under five in the population
#' #' @param outbreak_prob_range Range of outbreak probabilities
#' #' @return A data frame with one row per region and columns: id, pop_size, outbreak_prob
#' #' population and outbreak probability will be controlled in a more
#' #' sophisticated manner later
#' setup_region <- function(n_region = 1000,
#'                          pop_range = c(1000, 100000),
#'                          prop_U5_range = c(0.15, 0.18),
#'                          outbreak_prob_range = c(0.1, 0.9)) {
#'   tibble::tibble(
#'     id = 1:n_region,
#'     pop_size = sample(seq(pop_range[1], pop_range[2]), n_region, replace = TRUE),
#'     prop_U5 = runif(n_region, prop_U5_range[1], prop_U5_range[2]),
#'     outbreak_prob = runif(n_region, outbreak_prob_range[1], outbreak_prob_range[2])
#'   )
#' }
#'
#'
#' #' Generate outbreaks for each region
#' #'
#' #' This function simulates key outbreak parameters for each region:
#' #' start time \code{t0}, peak time \code{tp}, end time \code{t1}, baseline and peak
#' #' incidence rates \code{r0} and \code{rp}. It ensures that the total per-capita
#' #' incidence (area under the curve) remains less than 1.
#' #'
#' #' \deqn{ \text{Area} = \frac{1}{2}(t_1 - t_0)(r_p - r_0) \le 1 }
#' #'
#' #' Solving for the peak incidence rate:
#' #'
#' #' \deqn{ r_p \le r_0 + \frac{2}{t_1 - t_0} }
#' #'
#' #' @param regions Data frame with `id`, `pop_size`, `prop_U5`, `outbreak_prob`
#' #' @param peak_time_range Range to sample peak times (in days from t0)
#' #' @param duration_range Range to sample time from peak to end
#' #' @param peak_incidence_range Range of candidate values for peak incidence rate
#' #' @param max_attack_rate maximum attack rate allowed
#' #' @param r0 Baseline incidence rate (typically 0)
#' #'
#' #' @return Data frame with outbreak parameters and total cases
#' #' @export
#' generate_outbreaks <- function(regions,
#'                                start_time_range = c(0, 180),
#'                                peak_time_range = c(3*7, 20*7),
#'                                duration_range = c(5*7, 100*7),
#'                                peak_incidence_range = c(1e-6, 0.001), # peak incidence rate per capita
#'                                max_attack_rate = 0.1,
#'                                r0 = 0) {
#' df  <- purrr::map_dfr(1:nrow(regions), function(i) {
#'     region <- regions[i, ]
#'     if (runif(1) < region$outbreak_prob) {
#'       # t0 <- 0
#'       # tp <- runif(1, peak_time_range[1], peak_time_range[2])
#'       # t1 <- tp + runif(1, duration_range[1], duration_range[2])
#'
#'       t0 <- runif(1, start_time_range[1], start_time_range[2])
#'       tp <- t0 + runif(1, peak_time_range[1], peak_time_range[2])
#'       t1 <- tp + runif(1, duration_range[1], duration_range[2])
#'
#'       rp <- runif(1, peak_incidence_range[1], peak_incidence_range[2])
#'
#'       # Area of triangle: A = (1/2)*(base)*(height)
#'       base <- t1 - t0
#'       rp_max <- r0 + 2 / base * max_attack_rate
#'       rp_sample <- runif(1, peak_incidence_range[1], peak_incidence_range[2])
#'       rp <- min(rp_sample, rp_max)
#'       height <- rp - r0
#'       area <- 0.5 * base * height  # total cases per capita
#'
#'       tibble::tibble(
#'         id = region$id,
#'         t0 = t0,
#'         tp = tp,
#'         t1 = t1,
#'         r0 = r0,
#'         rp = rp,
#'         attack_rate_per_capita = area,
#'         total_cases = area * region$pop_size # ignore dynamic change of population susceptible
#'       )
#'     } else {
#'       return(NULL)
#'     }
#'   })
#'   df <- left_join(regions, df, by="id")
#' }


#' Generate stylized outbreak trajectories for each region
#'
#' This function simulates simplified outbreak dynamics at the region (grid-cell)
#' level. Outbreak occurrence is first determined by a Bernoulli trial with
#' region-specific outbreak probability. Conditional on an outbreak occurring,
#' the outbreak trajectory is represented by a triangular incidence curve
#' defined by a start time, peak time, and end time.
#'
#' Timing is parameterized using:
#' \itemize{
#'   \item a start time \eqn{t_0},
#'   \item a total outbreak duration \eqn{D = t_1 - t_0},
#'   \item a peak fraction \eqn{\phi \in (0,1)} such that
#'         \eqn{t_p = t_0 + \phi D}.
#' }
#'
#' This parameterization enforces the constraint \eqn{t_0 < t_p < t_1} by
#' construction and avoids rejection sampling.
#'
#' Incidence is assumed to increase linearly from a baseline rate \eqn{r_0}
#' at \eqn{t_0} to a peak incidence rate \eqn{r_p} at \eqn{t_p}, and then decrease
#' linearly back to \eqn{r_0} at \eqn{t_1}. The total per-capita attack rate
#' corresponds to the area under this triangular curve:
#'
#' \deqn{
#' A = \frac{1}{2} D (r_p - r_0).
#' }
#'
#' The peak incidence rate is truncated to ensure that the total attack rate
#' does not exceed \code{max_attack_rate}.
#'
#' @param regions Data frame containing at least:
#'   \code{id}, \code{pop_size}, and \code{outbreak_prob}.
#' @param start_time_range Range for outbreak start times.
#' @param total_duration_range Range for total outbreak duration.
#' @param peak_frac_shape Shape parameters (a, b) for Beta distribution of
#'   the peak fraction \eqn{\phi}.
#' @param peak_incidence_range Candidate range for peak incidence rate.
#' @param max_attack_rate Maximum per-capita attack rate.
#' @param r0 Baseline incidence rate (default = 0).
#'
#' @return Data frame with one row per region, including outbreak timing
#'   and severity parameters. Regions without outbreaks have NA values
#'   for outbreak-specific fields.
#'
#' @details
#' Outbreak shape is assumed independent of outbreak probability,
#' conditional on outbreak occurrence.
#'
#' @export
assign_outbreak_trajectory <- function(regions,
                               start_time_range = c(0, 180),
                               total_duration_range = c(10 * 7, 120 * 7),
                               peak_frac_shape = c(2, 2),
                               peak_incidence_range = c(1e-5, 5e-4),
                               cfr_range = c(0.001, 0.1),
                               max_attack_rate = 0.1,
                               r0 = 0) {

  purrr::map_dfr(seq_len(nrow(regions)), function(i) {

    region <- regions[i, ]

    ## --- outbreak occurrence ---
    # if (runif(1) >= region$outbreak_prob)
    if (!region$Y) # outbreak occurred
      return(
        tibble::tibble(
          id = region$id,
          t0 = NA_real_,
          tp = NA_real_,
          t1 = NA_real_,
          duration = NA_real_,
          peak_fraction = NA_real_,
          r0 = NA_real_,
          rp = NA_real_,
          attack_rate_per_capita = NA_real_,
          case_fatality_ratio = NA_real_
        )
      )

    ## --- timing ---
    t0 <- runif(1, start_time_range[1], start_time_range[2])
    D <- runif(1, total_duration_range[1], total_duration_range[2])
    phi <- rbeta(1, peak_frac_shape[1], peak_frac_shape[2])

    tp <- t0 + phi * D
    t1 <- t0 + D

    ## --- peak incidence with attack-rate constraint ---
    rp_max <- r0 + 2 * max_attack_rate / D

    rp <- min(
      runif(1, peak_incidence_range[1], peak_incidence_range[2]),
      rp_max
    )

    ## --- total attack rate ---
    area <- 0.5 * D * (rp - r0)
    cfr <- runif(1, cfr_range[1], cfr_range[2])

    tibble::tibble(
      id = region$id,
      t0 = t0,
      tp = tp,
      t1 = t1,
      duration = D,
      peak_fraction = phi,
      r0 = r0,
      rp = rp,
      attack_rate_per_capita = area,
      case_fatality_ratio = cfr
    )
  }) |>
    dplyr::left_join(regions, by = "id")
}



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

# ggplot2 theme for map
theme_map <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank()
    )
}

#' Generate Timestamp
#'
#' Creates a formatted timestamp with configurable components
#'
#' @param year Include year (default: TRUE)
#' @param month Include month (default: TRUE)
#' @param day Include day (default: TRUE)
#' @param hour Include hour (default: FALSE)
#' @param minute Include minute (default: FALSE)
#' @param second Include second (default: FALSE)
#' @return Formatted timestamp string
tstamp <- function(year=TRUE, month=TRUE, day=TRUE,
                   hour=FALSE, minute=FALSE, second=FALSE) {
  # Format date part
  date_format <- ""
  if (year) date_format <- paste0(date_format, "%Y")
  if (month) date_format <- paste0(date_format, "%m")
  if (day) date_format <- paste0(date_format, "%d")

  # Format time part
  time_format <- ""
  if (hour) time_format <- paste0(time_format, "%H")
  if (minute) time_format <- paste0(time_format, "%M")
  if (second) time_format <- paste0(time_format, "%S")

  # Create timestamp
  if (date_format == "" && time_format == "") {
    return("You'd better select parameters well.")
  }

  result <- if (date_format != "") format(Sys.time(), date_format) else ""

  if (time_format != "") {
    time_part <- format(Sys.time(), time_format)
    result <- if (result != "") paste0(result, "T", time_part) else time_part
  }

  return(result)
}
