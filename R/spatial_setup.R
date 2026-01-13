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
