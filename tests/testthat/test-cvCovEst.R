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
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst,
      spikedOperatorShrinkEst, spikedFrobeniusShrinkEst,
      spikedSteinShrinkEst
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
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L), noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L), noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L), noise = c(NULL, 0.5)
      )
    ),
    cv_scheme = "v_fold", mc_split = 0.5, v_folds = 5,
    parallel = FALSE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst,
      spikedOperatorShrinkEst, spikedFrobeniusShrinkEst,
      spikedSteinShrinkEst
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
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L), noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L), noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L), noise = c(NULL, 0.5)
      )
    ),
    cv_scheme = "mc", mc_split = 0.5, v_folds = 5,
    parallel = FALSE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst,
      spikedOperatorShrinkEst, spikedFrobeniusShrinkEst,
      spikedSteinShrinkEst
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
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L), noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L), noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L), noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cvFrobeniusLoss, cv_scheme = "v_fold",
    mc_split = 0.5, v_folds = 5,
    parallel = FALSE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst,
      spikedOperatorShrinkEst, spikedFrobeniusShrinkEst,
      spikedSteinShrinkEst
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
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L), noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L), noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L), noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cvFrobeniusLoss, cv_scheme = "mc",
    mc_split = 0.5, v_folds = 5,
    parallel = FALSE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst,
      spikedOperatorShrinkEst, spikedFrobeniusShrinkEst,
      spikedSteinShrinkEst
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
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L), noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L), noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L), noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cvScaledMatrixFrobeniusLoss, cv_scheme = "v_fold",
    mc_split = 0.5, v_folds = 5,
    parallel = FALSE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, robustPoetEst, adaptiveLassoEst,
      spikedOperatorShrinkEst, spikedFrobeniusShrinkEst,
      spikedSteinShrinkEst
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
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L), noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L), noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L), noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cvScaledMatrixFrobeniusLoss, cv_scheme = "mc",
    mc_split = 0.5, v_folds = 5,
    parallel = FALSE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(poetEst),
    estimator_params = list(
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L))
    ),
    cv_scheme = "v_fold", mc_split = 0.5, v_folds = 5,
    parallel = FALSE
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
    parallel = FALSE
  ))
})

test_that("cvCovEst's outputs are of the correct dimensions", {

  cv_ov_est_fit <- cvCovEst(
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
    parallel = FALSE
  )
  expect_true(length(cv_ov_est_fit) == 5)
  expect_true(ncol(cv_ov_est_fit$risk_df) == 3)
  expect_true(ncol(cv_ov_est_fit$cv_df) == 4)
  expect_true(is.null(cv_ov_est_fit$cv_oracle_riskdiff_ratio))
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
    parallel = TRUE
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
    parallel = TRUE
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
    parallel = TRUE
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
    parallel = TRUE
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
    parallel = TRUE
  ))
})
