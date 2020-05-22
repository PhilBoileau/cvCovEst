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
