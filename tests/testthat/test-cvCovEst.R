library(MASS)
set.seed(123)

# generate a 50x50 covariance matrix with unit variances and off-diagonal
# elements equal to 0.5
Sigma <- matrix(0.5, nrow = 50, ncol = 50) + diag(0.5, nrow = 50)

# sample 200 observations from multivariate normal with mean = 0, var = Sigma
dat <- mvrnorm(n = 200, mu = rep(0, 50), Sigma = Sigma)

# simple test (TODO: improve test and add tests to cover more cases)
test_that("cross-validated covariance selector runs silently", {
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(linearShrinkEst, linearShrinkLWEst,
                   thresholdingEst, sampleCovEst, bandingEst,
                   taperingEst, nlShrinkLWEst),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L))
    ),
    cv_scheme = "v_fold", mc_split = 0.5,
    v_folds = 5, boot_iter = 20,
    center = TRUE, scale = FALSE, parallel = FALSE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(linearShrinkEst, linearShrinkLWEst,
                   thresholdingEst, sampleCovEst, bandingEst,
                   taperingEst, nlShrinkLWEst),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L))
    ),
    cv_scheme = "mc", mc_split = 0.5,
    v_folds = 5, boot_iter = 20,
    center = TRUE, scale = FALSE, parallel = FALSE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(linearShrinkEst, linearShrinkLWEst,
                   thresholdingEst, sampleCovEst, bandingEst,
                   taperingEst, nlShrinkLWEst),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L))
    ),
    cv_scheme = "v_fold", mc_split = 0.5,
    v_folds = 5, boot_iter = 10,
    cv_loss = cvPenFrobeniusLoss,
    center = TRUE, scale = FALSE, parallel = FALSE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(linearShrinkEst, linearShrinkLWEst,
                   thresholdingEst, sampleCovEst, bandingEst,
                   taperingEst, nlShrinkLWEst),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L))
    ),
    cv_scheme = "mc", mc_split = 0.5,
    v_folds = 5, boot_iter = 10,
    cv_loss = cvPenFrobeniusLoss,
    center = TRUE, scale = FALSE, parallel = FALSE
  ))
})
