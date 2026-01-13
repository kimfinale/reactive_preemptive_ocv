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
