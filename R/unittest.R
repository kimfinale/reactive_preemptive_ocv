
# Source your get_area() function and get_IR() if not already loaded
library(testthat)
source("R/functions.R")

test_that("get_area() returns 0 before outbreak starts", {
  expect_equal(get_area(t = 0, t0 = 0, tp = 10, t1 = 30, r0 = 0, rp = 0.0), 0)
})

test_that("get_area() computes triangle area for t <= tp", {
  t0 <- 0; tp <- 10; t1 <- 30; r0 <- 0; rp <- 0.02
  t <- 5  # midway between t0 and tp
  expected_area <- 0.5 * (t - t0) * ((r0 + (rp - r0) * (t - t0)/(tp - t0)) - r0)
  expect_equal(get_area(t, t0, tp, t1, r0, rp), expected_area)
})

test_that("get_area() computes full triangle at t = tp", {
  t0 <- 0; tp <- 10; t1 <- 30; r0 <- 0; rp <- 0.02
  expected_area <- 0.5 * (tp - t0) * (rp - r0)
  expect_equal(get_area(tp, t0, tp, t1, r0, rp), expected_area)
})

test_that("get_area() computes triangle + trapezoid for t between tp and t1", {
  t0 <- 0; tp <- 10; t1 <- 30; r0 <- 0; rp <- 0.02
  t <- 20
  A1 <- 0.5 * (tp - t0) * (rp - r0)
  rt <- rp - (rp - r0) * (t - tp)/(t1 - tp)
  A2 <- 0.5 * (t - tp) * ((rp - r0) + (rt - r0))
  expect_equal(get_area(t, t0, tp, t1, r0, rp), A1 + A2)
})

test_that("get_area() returns full area for t > t1", {
  t0 <- 0; tp <- 10; t1 <- 30; r0 <- 0; rp <- 0.02
  A_total <- 0.5 * (tp - t0) * (rp - r0) + 0.5 * (t1 - tp) * (rp - r0)
  expect_equal(get_area(40, t0, tp, t1, r0, rp), A_total)
})
