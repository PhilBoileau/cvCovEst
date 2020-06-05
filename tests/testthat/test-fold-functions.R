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
  fold_fun = folds_montecarlo,
  V = 1, pvalidation = 0.3
)[[1]]

test_that("linearShrinkEst does not throw error", {
  expect_silent(cvFrobeniusLoss(
    fold = resub,
    dat = dat,
    estimator_funs = list("linearShrinkEst"),
    estimator_params =
      list(
        "linearShrinkEst" =
          list(alpha = c(0, 1))
      )
  ))
  expect_silent(cvPenFrobeniusLoss(
    fold = resub,
    dat = dat,
    resample_iter = 10,
    estimator_funs = list("linearShrinkEst"),
    estimator_params =
      list(
        "linearShrinkEst" =
          list(alpha = c(0, 1))
      )
  ))
})

test_that("thresholdEst does not throw error", {
  expect_silent(cvFrobeniusLoss(
    fold = resub,
    dat = dat,
    estimator_funs = list("thresholdingEst"),
    estimator_params =
      list(
        "thresholdingEst" =
          list(gamma = c(0, 1))
      )
  ))
  expect_silent(cvPenFrobeniusLoss(
    fold = resub,
    dat = dat,
    resample_iter = 10,
    estimator_funs = list("thresholdingEst"),
    estimator_params =
      list(
        "thresholdingEst" =
          list(alpha = c(0, 1))
      )
  ))
})

