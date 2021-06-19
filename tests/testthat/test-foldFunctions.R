library(MASS)
library(origami)
set.seed(123)

# generate 10x10 covariance matrix with unit variances and off-diagonal
# elements equal to 0.5
Sigma <- matrix(0.5, nrow = 10, ncol = 10) + diag(0.5, nrow = 10)

# sample 50 observations from multivariate normal with mean = 0, var = Sigma
dat <- mvrnorm(n = 50, mu = rep(0, 10), Sigma = Sigma)

# generate a single fold using MC-cv
resub <- make_folds(dat,
  fold_fun = folds_vfold,
  V = 5
)[[1]]

test_that("Estimators not throw error", {

  # get estimators expression call
  expect_silent(cvFrobeniusLoss(
    fold = resub,
    dat = dat,
    estimator_funs = rlang::quo(c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst,
      robustPoetEst
    )),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0, 1)),
      thresholdingEst = list(gamma = c(0, 1)),
      robustPoetEst = list(
        lambda = 0.1, k = 1L,
        var_est = c("mad")
      )
    )
  ))
  expect_silent(cvMatrixFrobeniusLoss(
    fold = resub,
    dat = dat,
    estimator_funs = rlang::quo(c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst
    )),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0, 1)),
      thresholdingEst = list(gamma = c(0, 1))
    )
  ))
  expect_silent(cvScaledMatrixFrobeniusLoss(
    fold = resub,
    dat = dat,
    estimator_funs = rlang::quo(c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst
    )),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0, 1)),
      thresholdingEst = list(gamma = c(0, 1))
    )
  ))
})

test_that("trueFrobeniusLoss computes the correct validation set loss", {
  est <- matrix(c(1, 2, 3, 2, 1, 2, 3, 2, 1), nrow = 3)
  true_cov <- matrix(c(1, 2, 2, 2, 1, 2, 2, 2, 1), nrow = 3)
  expect_equal(trueFrobeniusLoss(est, true_cov), 38)
})

test_that("cvFrobeniusLoss returns the expected error for the linear shrinkage
          estimator with alpha 0.5", {

  # compute the error for the given fold using the fold function
  computed_error <- cvFrobeniusLoss(
    fold = resub,
    dat = dat,
    estimator_funs = rlang::quo(c(linearShrinkEst)),
    estimator_params = list(linearShrinkEst = list(alpha = 0.5))
  )[[1]]$loss

  # now compute the error manually:

  # fit the data to the training set
  estimate_train <- linearShrinkEst(dat[resub$training_set, ], alpha = 0.5)

  # define the validation data
  valid_data <- dat[resub$validation_set, ]

  # compute the validation observation crossprod matrices
  rank_one_crossp <- lapply(
    seq_len(nrow(valid_data)),
    function(i) {
      Matrix::tcrossprod(valid_data[i, ])
    }
  )

  # compute the sum of element-wise squared outer products
  elwise_sq_list <- lapply(rank_one_crossp, `^`, 2)
  elwise_sq_sum <- Reduce(`+`, elwise_sq_list)
  elwise_sq <- matrixStats::sum2(elwise_sq_sum)

  # compute the sum of cross_products
  cross_prod <- Reduce(`+`, rank_one_crossp)

  # get the hadamard product and sum
  had_crossprod <- matrixStats::sum2(cross_prod * estimate_train)

  # get the elementwise square and sum it
  est_square <- matrixStats::sum2(estimate_train^2)

  # compute the error
  estimate_error <- 1 / nrow(valid_data) * (elwise_sq - 2 * had_crossprod) +
    est_square

  testthat::expect_identical(estimate_error, computed_error)

})

test_that("cvMatrixFrobeniusLoss returns the expected error for the hard
          thresholding estimator with gamma of 1", {

  # compute the error for the given fold using the fold function
  computed_error <- cvMatrixFrobeniusLoss(
    fold = resub,
    dat = dat,
    estimator_funs = rlang::quo(c(thresholdingEst)),
    estimator_params = list(thresholdingEst = list(gamma = 1))
  )[[1]]$loss

  # now compute the error manually:

  # fit the data to the training set
  estimate_train <- thresholdingEst(dat[resub$training_set, ], gamma = 1)

  # define the validation data
  valid_data <- dat[resub$validation_set, ]

  # compute the sample covariance over the validation set
  sample_cov_valid <- coop::covar(valid_data)

  # compute the error
  estimate_error <- matrixStats::sum2((estimate_train - sample_cov_valid)^2)

  testthat::expect_identical(estimate_error, computed_error)

})


test_that("the true covariance matrix can be passed in", {
  expect_silent(cvFrobeniusLoss(
    fold = resub,
    dat = dat,
    estimator_funs = rlang::quo(c(sampleCovEst, poetEst)),
    estimator_params = list(poetEst = list(k = 5L, lambda = 0.2)),
    true_cov_mat = Sigma
  ))
  output <- cvFrobeniusLoss(
    fold = resub,
    dat = dat,
    estimator_funs = rlang::quo(c(sampleCovEst, poetEst)),
    estimator_params = list(poetEst = list(k = 5L, lambda = 0.2)),
    true_cov_mat = Sigma
  )[[1]]
  expect_equal(nrow(output), 2)
  expect_equal(ncol(output), 6)
  expect_true(!is.null(output$true_loss))
})
