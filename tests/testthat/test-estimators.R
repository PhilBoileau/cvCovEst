# Linear Shrinkage Estimator ##################################################
test_that("linear shrinkage estimator with no shrinkage is S_n", {
  expect_identical(
    linearShrinkEst(mtcars, alpha = 1),
    coop::covar(mtcars)
  )
})

test_that("linear shrinkage estimator with 100% shrinkage is identity", {
  expect_identical(
    linearShrinkEst(mtcars, alpha = 0) %>% unname(),
    diag(ncol(mtcars))
  )
})


# Ledoit Wolf Linear Shrinkage Estimator #######################################
test_that("LW LS estimator runs without issue", {
  expect_silent(
    linearShrinkLWEst(mtcars)
  )
})

# Simple Thresholding Estimator ###############################################
test_that("simple thresholing estimator without thresholding is S_n", {
  expect_identical(
    thresholdingEst(mtcars, gamma = 0),
    coop::covar(mtcars)
  )
})

test_that("simple thresholing estimator with large threshold is 0 matrix", {
  expect_identical(
    thresholdingEst(mtcars, gamma = 1000000) %>% unname(),
    matrix(data = 0, nrow = ncol(mtcars), ncol = ncol(mtcars))
  )
})

# Banding Estimator ###########################################################
test_that("banding estimator with k = 0 is diagonal of S_n", {
  expect_identical(
    bandingEst(mtcars, k = 0L) %>% unname(),
    diag( diag( coop::covar(mtcars) ) )
  )
})

test_that("banding estimator with k >> 0 is S_n", {
  expect_identical(
    bandingEst(mtcars, k = 1000000L),
    coop::covar(mtcars)
  )
})

# Tapering Estimator ##########################################################
test_that("tapering estimator with k = 0 is diagonal of S_n", {
  expect_identical(
    taperingEst(mtcars, k = 0L) %>% unname(),
    diag( diag( coop::covar(mtcars) ) )
  )
})

test_that("tapering estimator with k = 2*p is S_n", {
  expect_identical(
    taperingEst(mtcars, k = 22L) %>% unname(),
    coop::covar(mtcars)
  )
})

test_that("tapering estimator with k >> 0 is S_n", {
  expect_identical(
    taperingEst(mtcars, k = 1000000L) %>% unname(),
    coop::covar(mtcars)
  )
})

