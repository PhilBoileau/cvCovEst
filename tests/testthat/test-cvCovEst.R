library(future.apply)
library(MASS)
set.seed(123)

# generate a 10x10 covariance matrix with unit variances and off-diagonal
# elements equal to 0.5
Sigma <- matrix(0.5, nrow = 10, ncol = 10) + diag(0.5, nrow = 10)

# sample 50 observations from multivariate normal with mean = 0, var = Sigma
dat <- mvrnorm(n = 50, mu = rep(0, 10), Sigma = Sigma)

test_that("cross-validated covariance selector runs silently", {
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst
    ),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 2L),
        var_est = c("sample")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5))
    ),
    cv_scheme = "v_fold", mc_split = 0.5, v_folds = 5,
    center = TRUE, scale = FALSE, parallel = FALSE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst
    ),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 2L),
        var_est = c("sample")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5))
    ),
    cv_scheme = "mc", mc_split = 0.5, v_folds = 5,
    center = TRUE, scale = FALSE, parallel = FALSE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst
    ),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 2L),
        var_est = c("sample")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5))
    ),
    cv_loss = cvFrobeniusLoss, cv_scheme = "v_fold",
    mc_split = 0.5, v_folds = 5,
    center = TRUE, scale = FALSE, parallel = FALSE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst
    ),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 2L),
        var_est = c("sample")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5))
    ),
    cv_loss = cvFrobeniusLoss, cv_scheme = "mc",
    mc_split = 0.5, v_folds = 5,
    center = TRUE, scale = FALSE, parallel = FALSE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst
    ),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 2L),
        var_est = c("sample")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5))
    ),
    cv_loss = cvScaledMatrixFrobeniusLoss, cv_scheme = "v_fold",
    mc_split = 0.5, v_folds = 5,
    center = TRUE, scale = FALSE, parallel = FALSE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst
    ),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 2L),
        var_est = c("sample")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5))
    ),
    cv_loss = cvScaledMatrixFrobeniusLoss, cv_scheme = "mc",
    mc_split = 0.5, v_folds = 5,
    center = TRUE, scale = FALSE, parallel = FALSE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst
    ),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 2L),
        var_est = c("sample")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5))
    ),
    cv_loss = cvFrobeniusLoss,
    cv_scheme = "v_fold", mc_split = 0.5, v_folds = 5,
    center = TRUE, scale = FALSE, parallel = FALSE,
    true_cov_mat = Sigma
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(poetEst),
    estimator_params = list(
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L))
    ),
    cv_scheme = "v_fold", mc_split = 0.5, v_folds = 5,
    center = TRUE, scale = FALSE, parallel = FALSE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(robustPoetEst),
    estimator_params = list(
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 2L),
        var_est = "sample"
      )
    ),
    cv_scheme = "v_fold", mc_split = 0.5, v_folds = 5,
    center = TRUE, scale = FALSE, parallel = FALSE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(adaptiveLassoEst),
    estimator_params = list(
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5))
    ),
    cv_loss = cvFrobeniusLoss,
    cv_scheme = "v_fold", mc_split = 0.5, v_folds = 5,
    center = TRUE, scale = FALSE, parallel = FALSE,
    true_cov_mat = Sigma
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkLWEst, thresholdingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, adaptiveLassoEst
    ),
    estimator_params = list(
      thresholdingEst = list(gamma = c(0.2, 2)),
      scadEst = list(lambda = c(0.1, 0.2)),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5))
    ),
    cv_loss = cvFrobeniusLoss,
    cv_scheme = "v_fold", mc_split = 0.5, v_folds = 5,
    center = TRUE, scale = FALSE, parallel = FALSE,
    true_cov_mat = Sigma
  ))
})

test_that("cvCovEst automatically centers non-centered data", {
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst
    ),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 2L),
        var_est = c("sample")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5))
    ),
    cv_scheme = "v_fold", mc_split = 0.5, v_folds = 5,
    center = TRUE, scale = FALSE, parallel = FALSE
  ))
  expect_message(
    cvCovEst(
      dat = dat,
      estimators = c(
        linearShrinkEst, linearShrinkLWEst,
        thresholdingEst, sampleCovEst, bandingEst,
        taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
        scadEst, poetEst, robustPoetEst, adaptiveLassoEst
      ),
      estimator_params = list(
        linearShrinkEst = list(alpha = c(0.1, 0.9)),
        thresholdingEst = list(gamma = c(0.2, 2)),
        bandingEst = list(k = c(1L, 5L)),
        taperingEst = list(k = c(2L, 6L)),
        scadEst = list(lambda = c(0.1, 0.2)),
        poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L)),
        robustPoetEst = list(
          lambda = c(0.1, 0.2), k = c(1L, 2L),
          var_est = c("sample")
        ),
        adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5))
      ),
      cv_scheme = "v_fold", mc_split = 0.5, v_folds = 5,
      center = FALSE, scale = FALSE, parallel = FALSE
    ),
    "The columns of argument `dat` have been centered automatically"
  )
})

test_that("cvCovEst's outputs are of the correct dimensions", {
  with_true_covmat <- cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst
    ),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 2L),
        var_est = c("sample", "mad")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5))
    ),
    cv_loss = cvFrobeniusLoss,
    cv_scheme = "v_fold", mc_split = 0.5, v_folds = 5,
    center = TRUE, scale = FALSE, parallel = FALSE,
    true_cov_mat = Sigma
  )
  expect_true(length(with_true_covmat) == 9)
  expect_true(ncol(with_true_covmat$risk_df) == 5)
  expect_true(ncol(with_true_covmat$cv_df) == 6)
  expect_true(is.numeric(with_true_covmat$cv_cv_riskdiff))
  expect_true(is.numeric(with_true_covmat$oracle_cv_riskdiff))
  expect_true(is.numeric(with_true_covmat$cv_oracle_riskdiff_ratio))
  expect_true(is.numeric(with_true_covmat$full_oracle_riskdiff_ratio))
  expect_equal(
    with_true_covmat$cv_cv_riskdiff / with_true_covmat$oracle_cv_riskdiff,
    with_true_covmat$cv_oracle_riskdiff_ratio
  )

  without_true_covmat <- cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst
    ),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 2L),
        var_est = c("sample")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5))
    ),
    cv_scheme = "v_fold", mc_split = 0.5, v_folds = 5,
    center = TRUE, scale = FALSE, parallel = FALSE
  )
  expect_true(length(without_true_covmat) == 5)
  expect_true(ncol(without_true_covmat$risk_df) == 3)
  expect_true(ncol(without_true_covmat$cv_df) == 4)
  expect_true(is.null(without_true_covmat$cv_oracle_riskdiff_ratio))
})

test_that("cross-validated covariance selector handles sparse, true cov mat", {
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst
    ),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 2L),
        var_est = c("sample")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5))
    ),
    cv_loss = cvFrobeniusLoss,
    cv_scheme = "v_fold", mc_split = 0.5, v_folds = 5,
    center = TRUE, scale = FALSE, parallel = FALSE,
    true_cov_mat = as(Sigma, "dgeMatrix")
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst
    ),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 2L),
        var_est = c("sample")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5))
    ),
    cv_loss = cvFrobeniusLoss, cv_scheme = "v_fold",
    mc_split = 0.5, v_folds = 5,
    center = TRUE, scale = FALSE, parallel = FALSE,
    true_cov_mat = as(Sigma, "dgCMatrix")
  ))
})

test_that("Parallelization works", {
  plan(sequential)
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst
    ),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 2L),
        var_est = c("sample")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5))
    ),
    cv_scheme = "v_fold", mc_split = 0.5, v_folds = 5,
    center = TRUE, scale = FALSE, parallel = TRUE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst
    ),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 2L),
        var_est = c("sample")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5))
    ),
    cv_scheme = "mc", mc_split = 0.5, v_folds = 5,
    center = TRUE, scale = FALSE, parallel = TRUE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst
    ),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 2L),
        var_est = c("sample")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5))
    ),
    cv_loss = cvFrobeniusLoss, cv_scheme = "v_fold",
    mc_split = 0.5, v_folds = 5,
    center = TRUE, scale = FALSE, parallel = TRUE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst
    ),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 2L),
        var_est = c("sample")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5))
    ),
    cv_loss = cvFrobeniusLoss, cv_scheme = "mc",
    mc_split = 0.5, v_folds = 5,
    center = TRUE, scale = FALSE, parallel = TRUE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst
    ),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 2L),
        var_est = c("sample")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5))
    ),
    cv_loss = cvScaledMatrixFrobeniusLoss, cv_scheme = "mc",
    mc_split = 0.5, v_folds = 5,
    center = TRUE, scale = FALSE, parallel = TRUE
  ))
})
