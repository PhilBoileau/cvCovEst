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
    diag(diag(coop::covar(mtcars)))
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
    diag(diag(coop::covar(mtcars)))
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

  # make sure that nlShrinkLWEst can handle case where n = p

  library(MASS)
  set.seed(123)
  Sigma <- matrix(0.5, nrow = 50, ncol = 50) + diag(0.5, nrow = 50)
  dat <- mvrnorm(n = 50, mu = rep(0, 50), Sigma = Sigma)
  expect_false(
    any(is.na(nlShrinkLWEst(dat)))
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

# SCAD Thresholding Estimator ##################################################

test_that("SCAD estimator doesn't generate any errors for no reason", {
  expect_silent(
    scadEst(mtcars, lambda = 0.1)
  )
})

test_that("SCAD estimator generates zero matrix for large lambda", {
  dat <- scale(mtcars, center = TRUE, scale = TRUE)
  expect_equal(
    sum(scadEst(dat, lambda = 10) == 0),
    121
  )
})

# POET Estimator ###############################################################

test_that("Verify POET estimator's results", {
  # check that the off diagonal elemets of the POET estimator equal to the
  # off diagonal elements of the rank-1 approximation of the sample covariance
  # matrix

  # compute the rank-1 approx
  dat <- scale(mtcars, center = TRUE, scale = TRUE)
  sample_cov_mat <- coop::covar(dat)
  eig_decomp <- RSpectra::eigs_sym(sample_cov_mat, 1)
  rank_one_sample_cov <- eig_decomp$values *
    eig_decomp$vectors %*% t(eig_decomp$vectors)

  # compute the POET estimate with a large lambda
  poet_estimate <- poetEst(dat, k = 1, lambda = 10)

  # remove the diagonal
  diag(rank_one_sample_cov) <- 0
  diag(poet_estimate) <- 0
  poet_estimate <- poet_estimate %>% unname()

  # compare
  expect_identical(round(rank_one_sample_cov, 10), round(poet_estimate, 10))
})

# Robust POET Estimator #########################################################

test_that("Verify Robust POET estimator's ranks", {
  # check that the rank of robust POET equal to k when lambda is large
  # compute the robust POET estimate with a large lambda
  dat <- scale(mtcars, center = TRUE, scale = TRUE)
  k = ceiling(ncol(dat) / 5)
  robust_poet_estimate <- robustPoetEst(dat, k, lambda = 10, var_estimation = "sample")
  
  # compare
  library(Matrix)
  expect_equal(Matrix::rankMatrix(robust_poet_estimate)[1], k)
})


test_that("Verify Robust POET estimator's results", {
  # In this specific example, robust POET estimator is supposed to return a matrix of 
  # all 13 or all 19.78292 dependent on the variance estimation method
  Y = matrix(c(1:12, 2:13, 3:14, 4:15, 5:16), 12, 5)
  est = matrix(rep(13L, 25), 5, 5)
  est_mad = matrix(rep((3 * 1.4826) ** 2, 25), 5, 5)
  robust_poet_estimate = robustPoetEst(Y, 1, lambda = 10, var_estimation = "sample")
  robust_poet_estimate_mad = robustPoetEst(Y, 1, lambda = 10, var_estimation = "mad")
  # compare
  expect_equal(est, robust_poet_estimate)
  expect_equal(est_mad, robust_poet_estimate_mad)
})

# Adaptive Lasso Estimator #####################################################

test_that("adaptive Lasso estimator with no penalty is Sn", {
  expect_identical(
    adaptiveLassoEst(mtcars, lambda = 0, n = 0) %>% unname(),
    coop::covar(mtcars) %>% unname()
  )
})

test_that("adaptive Lasso estimator with no penalty and non-zero n is Sn", {
  expect_identical(
    adaptiveLassoEst(mtcars, lambda = 0, n = 1) %>% unname(),
    coop::covar(mtcars) %>% unname()
  )
})

test_that("adaptive Lasso estimator with large threshold is 0 matrix", {
  expect_identical(
    adaptiveLassoEst(mtcars, lambda = 1000000, n = 0) %>% unname(),
    matrix(data = 0, nrow = ncol(mtcars), ncol = ncol(mtcars))
  )
})
