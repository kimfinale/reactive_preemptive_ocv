# montecarlo.R
# This file implements a single Monte Carlo draw for testing simulations

get_sample <- function(wait, supply, Q, tat, c, f, T, L,
                       wait_exit, Q_exit, f_exit, L_exit,
                       kinetics, testing, isolation) {
  params <- kinetics()
  names(params) <- c("A", "P", "B", "tSx", "m", "M", "thr")
  with(as.list(params), {
    D <- get_D(A, P, B, m, M, L)
    first <- D[1]; last <- D[2]

    tests_scheduled <- get_scheduled_tests(tSx, wait, supply, Q, T, testing)
    tests_taken <- apply_compliance(tests_scheduled, c)
    valid_tests_taken <- apply_failure(tests_taken, f)
    hits <- get_hits(valid_tests_taken, first, last)
    tDx <- get_tDx(hits, tat)
    n_tests <- count_tests(tests_taken, tDx)
    I0 <- area_triangle(A, P, B, m, M, thr)

    if (is.infinite(tDx)) {
      Itest <- I0
      tExit <- -1
      n_tests_exit <- 0
      Iexit <- 0
    } else {
      Itest <- area_diagnosis(A, P, B, m, M, thr, tDx)
      exit_info <- compute_exit(isolation, tDx, wait_exit, Q_exit, f_exit, L_exit, T,
                                A, P, B, m, M, thr)
      tExit <- exit_info[[1]]
      n_tests_exit <- exit_info[[2]]
      Iexit <- I0 - exit_info[[3]]
    }

    return(c(tDx, tExit, n_tests, n_tests_exit, I0, Itest, Iexit,
             A, P, B, tSx, m, M, first, last, thr))
  })
}

get_D <- function(A, P, B, m, M, L) {
  if (L > M) return(c(0, 0))
  if (L < m) return(c(A, B))
  f <- (L - m) / (M - m)
  a <- f * (P - A) + A
  b <- f * (P - B) + B
  return(c(a, b))
}

get_t_trigger <- function(A, threshold) {
  if (L > M) return(c(0, 0))
  if (L < m) return(c(A, B))
  f <- (L - m) / (M - m)
  a <- f * (P - A) + A
  b <- f * (P - B) + B
  # get t at which area is larger than the threshold
  #
  return(t_trigger)
}


get_scheduled_tests <- function(tSx, wait, supply, Q, T, method) {
  if (identical(method, test_regular)) {
    return(method(Q, T))
  } else if (identical(method, test_post_exposure)) {
    return(method(wait, supply, Q))
  } else if (identical(method, test_post_symptoms)) {
    return(method(tSx, wait, supply, Q))
  }
}

test_post_symptoms <- function(tSx, wait, supply, Q) {
  phase <- Q * runif(1)
  tSx + wait + phase + seq(0, Q * (supply - 1), by = Q)
}

test_post_exposure <- function(wait, supply, Q) {
  phase <- Q * runif(1)
  wait + phase + seq(0, Q * (supply - 1), by = Q)
}

test_regular <- function(Q, T) {
  phase <- Q * runif(1)
  seq(phase, T, by = Q)
}

apply_compliance <- function(tests_scheduled, c) {
  n <- rbinom(1, length(tests_scheduled), c)
  if (n == 0) return(numeric(0))
  sort(sample(tests_scheduled, n))
}

apply_failure <- function(tests_taken, f) {
  n <- rbinom(1, length(tests_taken), 1 - f)
  if (n == 0) return(numeric(0))
  sort(sample(tests_taken, n))
}

get_hits <- function(valid_tests_taken, first, last) {
  valid_tests_taken[valid_tests_taken >= first & valid_tests_taken <= last]
}

get_tDx <- function(hits, tat) {
  if (length(hits) == 0) return(Inf)
  hits[1] + tat
}

count_tests <- function(tests_taken, tDx) {
  if (tDx == -1) return(length(tests_taken))
  sum(tests_taken <= tDx)
}

area_triangle <- function(a, b, thr) {
  if (M < thr) return(0)
  D <- get_D(A, P, B, m, M, thr)
  base <- D[2] - D[1]
  height <- M - thr
  base * height / 2
}

area_vx <- function(A, P, B, m, M, thr, tVE) {
  if (M < thr) return(0)
  D <- get_D(A, P, B, m, M, thr)
  a <- D[1]; b <- D[2]

  if (tDx < a) return(0)
  if (tDx <= P) {
    base <- tVE - a
    height <- base * (M - thr) / (P - a)
    return(base * height / 2)
  }
  if (tDx < b) {
    base <- b - tVE
    height <- base * (M - thr) / (b - P)
    return(((b - a) * (M - thr) - base * height) / 2)
  }
  return((b - a) * (M - thr) / 2)
}

get_scheduled_exit_tests <- function(tDx, wait_exit, Q_exit, T) {
  phase <- Q_exit * runif(1)
  t_start <- tDx + wait_exit + phase
  scheduled <- seq(t_start, T, by = Q_exit)
  if (length(scheduled) == 0) return(tDx + wait_exit)
  scheduled
}

get_tExit <- function(exit_tests_failed, exit_tests_scheduled, last) {
  if (length(exit_tests_failed) > 0) return(exit_tests_failed[1])
  exit_tests_scheduled[which(exit_tests_scheduled > last)[1]]
}

compute_exit <- function(isolation, tDx, wait_exit, Q_exit, f_exit, L_exit, T,
                         A, P, B, m, M, thr) {
  if (isolation == "fixed") {
    tExit <- tDx + wait_exit
    n_tests_exit <- 0
    Iexit_complement <- area_diagnosis(A, P, B, m, M, thr, tExit)
  } else if (isolation == "TTE") {
    exit_tests_scheduled <- get_scheduled_exit_tests(tDx, wait_exit, Q_exit, T)
    exit_tests_taken <- exit_tests_scheduled
    D <- get_D(A, P, B, m, M, L_exit)
    first <- D[1]; last <- D[2]
    exit_tests_positive <- get_hits(exit_tests_taken, first, last)
    positive_tests_failed <- apply_failure(exit_tests_positive, 1 - f_exit)
    tExit <- get_tExit(positive_tests_failed, exit_tests_scheduled, last)
    n_tests_exit <- count_tests(exit_tests_taken, tExit)
    Iexit_complement <- area_diagnosis(A, P, B, m, M, thr, tExit)
  }
  return(list(tExit, n_tests_exit, Iexit_complement))
}
