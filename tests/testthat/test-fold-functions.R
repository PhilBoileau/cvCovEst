library(MASS)
library(origami)
set.seed(123)

# generate 50x50 covariance matrix with unit variances and off-diagonal
# elements equal to 0.5
Sigma <- matrix(0.5, nrow = 50, ncol = 50) + diag(0.5, nrow = 50)

# sample 200 observations from multivariate normal with mean = 0, var = Sigma
dat <- mvrnorm(n = 200, mu = rep(0, 50), Sigma = Sigma)

# generate a single fold using MC-cv
resub <- make_folds(dat,
  fold_fun = folds_vfold,
  V = 5
)[[1]]

test_that("linearShrinkEst does not throw error", {

  # get estimators expression call
  expect_silent(cvFrobeniusLoss(
    fold = resub,
    dat = dat,
    estimator_funs = rlang::expr(c(linearShrinkEst)),
    estimator_params = list(linearShrinkEst = list(alpha = c(0, 1)))
  ))
})

test_that("linearShrinkLWEst does not throw error", {
  expect_silent(cvFrobeniusLoss(
    fold = resub,
    dat = dat,
    estimator_funs = rlang::expr(c(linearShrinkLWEst)),
  ))
})

test_that("thresholdEst does not throw error", {
  expect_silent(cvFrobeniusLoss(
    fold = resub,
    dat = dat,
    estimator_funs = rlang::expr(c(thresholdingEst)),
    estimator_params = list(thresholdingEst = list(gamma = c(0, 1)))
  ))
})

test_that("sampleCovEst does not throw error", {
  expect_silent(cvFrobeniusLoss(
    fold = resub,
    dat = dat,
    estimator_funs = rlang::expr(c(sampleCovEst)),
  ))
})

test_that("trueFrobeniusLoss computes the correct validation set loss", {
  est <- matrix(c(1, 2, 3, 2, 1, 2, 3, 2, 1), nrow = 3)
  true_cov <- matrix(c(1, 2, 2, 2, 1, 2, 2, 2, 1), nrow = 3)
  expect_equal(trueFrobeniusLoss(est, true_cov), 38)
})

test_that("the true covariance matrix can be passed in", {
  expect_silent(cvFrobeniusLoss(
    fold = resub,
    dat = dat,
    estimator_funs = rlang::expr(c(sampleCovEst, poetEst)),
    estimator_params = list(poetEst = list(k = 5L, lambda = 0.2)),
    true_cov_mat = Sigma
  ))
  output <- cvFrobeniusLoss(
    fold = resub,
    dat = dat,
    estimator_funs = rlang::expr(c(sampleCovEst, poetEst)),
    estimator_params = list(poetEst = list(k = 5L, lambda = 0.2)),
    true_cov_mat = Sigma
  )[[1]]
  expect_equal(nrow(output), 2)
  expect_equal(ncol(output), 5)
  expect_true(!is.null(output$true_loss))
})
