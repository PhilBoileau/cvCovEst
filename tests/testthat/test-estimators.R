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

test_that("Linear shrinkage estimator produces shrunken estimates", {

  # aAl absolute elements of the estimate should be smaller than the
  # absolute value of the entries in the sample covariance matrix estimate
  dat <- scale(mtcars, center = TRUE, scale = TRUE)
  abs_est <- abs(linearShrinkEst(dat, alpha = 0.1))
  abs_sample_cov <- abs(cov(dat))
  expect_true(
    sum(round((abs_sample_cov - abs_est), digits = 6) >= 0) == 121
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
    diag(diag( coop::covar(mtcars)))
  )
})

test_that("tapering estimator with k = 2*p-2 is S_n", {
  expect_identical(
    taperingEst(mtcars, k = 20L) %>% unname(),
    coop::covar(mtcars) %>% unname()
  )
})

test_that("tapering estimator with k = 2*p-4 is not S_n", {
  expect_false(identical(
    taperingEst(mtcars, k = 18L) %>% unname(),
    coop::covar(mtcars) %>% unname()
  ))
})

test_that("tapering estimator with k >> 0 is S_n", {
  expect_identical(
    taperingEst(mtcars, k = 1000000L) %>% unname(),
    coop::covar(mtcars) %>% unname()
  )
})

# Ledoit Wolf Nonlinear Shrinkage Estimator ####################################
test_that("LW NLS estimator runs without issue", {
  expect_silent(
    nlShrinkLWEst(mtcars)
  )
})

# Dense linear Shrinkage Estimator ####################################
test_that("Dense linear shrinkage estimator runs without issue", {
  expect_silent(
    denseLinearShrinkEst(mtcars)
  )
})

test_that("Dense linear shrinkage estimator produces shrunken estimates", {

  # Mean covariance is -0.05 in sample covarian matrix.
  # All absolute covariance values are larger than 0.057.
  # Estimator should therefore produce smaller estimates in each entry.
  dat <- scale(mtcars, center = TRUE, scale = TRUE)
  abs_est <- abs(denseLinearShrinkEst(dat))
  abs_sample_cov <- abs(cov(dat))
  expect_true(
    sum(round((abs_sample_cov - abs_est), digits = 6) >= 0) == 121
  )
})
