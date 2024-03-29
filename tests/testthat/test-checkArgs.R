library(MASS)
set.seed(123)

# generate a 50x50 covariance matrix with unit variances and off-diagonal
# elements equal to 0.5
Sigma <- matrix(0.5, nrow = 50, ncol = 50) + diag(0.5, nrow = 50)

# sample 200 observations from multivariate normal with mean = 0, var = Sigma
dat <- mvrnorm(n = 200, mu = rep(0, 50), Sigma = Sigma)

# define the arguments as they appear inside cvCovEst
estimators <- rlang::expr(c(
  linearShrinkEst, thresholdingEst, sampleCovEst,
  linearShrinkLWEst, bandingEst, taperingEst,
  nlShrinkLWEst, denseLinearShrinkEst, scadEst,
  poetEst, robustPoetEst, adaptiveLassoEst,
  spikedOperatorShrinkEst, spikedFrobeniusShrinkEst,
  spikedSteinShrinkEst
))
estimator_params <- list(
  linearShrinkEst = list(alpha = c(0.1, 0.9)),
  thresholdingEst = list(gamma = c(0.2, 2)),
  bandingEst = list(k = c(1L, 5L)),
  taperingEst = list(k = c(2L, 6L)),
  scadEst = list(lambda = c(0.1, 0.2)),
  poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L, 3L)),
  adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
  robustPoetEst = list(
    lambda = c(0.1, 0.2), k = c(1L, 3L),
    var_est = "sample"
  ),
  spikedOperatorShrinkEst = list(
    p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L), noise = c(NULL, 0.5)
  ),
  spikedFrobeniusShrinkEst = list(
    p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L), noise = c(NULL, 0.5)
  ),
  spikedSteinShrinkEst = list(
    p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L), noise = c(NULL, 0.5)
  )
)
cv_loss <- rlang::expr(cvFrobeniusLoss)
cv_scheme <- "mc"
mc_split <- 0.5
v_folds <- 10
parallel <- FALSE

test_that("Only implmented estimators pass checks", {
  expect_true(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_true(checkArgs(
    dat = dat,
    estimators = rlang::expr(c(linearShrinkEst)),
    estimator_params = estimator_params,
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = rlang::expr(c(linearShrinkEst, newEstimator)),
    estimator_params = estimator_params,
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ), "Only estimators implemented in the cvCovEst package can be used.")
})

test_that("Lone estimators without hyperparams aren't permitted", {
  expect_error(
    checkArgs(
      dat = dat,
      estimators = rlang::expr(c(sampleCovEst)),
      estimator_params = estimator_params,
      cv_loss = cv_loss,
      cv_scheme = cv_scheme,
      mc_split = mc_split,
      v_folds = v_folds,
      parallel = parallel
    ),
    "This estimator doesn't possess any hyperparameters. Run it without using cvCovEst."
  )
  expect_error(
    checkArgs(
      dat = dat,
      estimators = rlang::expr(c(nlShrinkLWEst)),
      estimator_params = estimator_params,
      cv_loss = cv_loss,
      cv_scheme = cv_scheme,
      mc_split = mc_split,
      v_folds = v_folds,
      parallel = parallel
    ),
    "This estimator doesn't possess any hyperparameters. Run it without using cvCovEst."
  )
  expect_error(
    checkArgs(
      dat = dat,
      estimators = rlang::expr(c(linearShrinkLWEst)),
      estimator_params = estimator_params,
      cv_loss = cv_loss,
      cv_scheme = cv_scheme,
      mc_split = mc_split,
      v_folds = v_folds,
      parallel = parallel
    ),
    "This estimator doesn't possess any hyperparameters. Run it without using cvCovEst."
  )
  expect_error(
    checkArgs(
      dat = dat,
      estimators = rlang::expr(c(denseLinearShrinkEst)),
      estimator_params = estimator_params,
      cv_loss = cv_loss,
      cv_scheme = cv_scheme,
      mc_split = mc_split,
      v_folds = v_folds,
      parallel = parallel
    ),
    "This estimator doesn't possess any hyperparameters. Run it without using cvCovEst."
  )
})

test_that("Only reasonable hyperparameters pass checks", {
  expect_true(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_true(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(-0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L, 3L)),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 3L),
        var_est = "sample"
      ),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 1.1)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L, 3L)),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 3L),
        var_est = "sample"
      ),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(-0.1, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L, 3L)),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 3L),
        var_est = "sample"
      ),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1.1, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L, 3L)),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 3L),
        var_est = "sample"
      ),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(-2L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L, 3L)),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 3L),
        var_est = "sample"
      ),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(2L, 5L)),
      taperingEst = list(k = c(2.1, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L, 3L)),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 3L),
        var_est = "sample"
      ),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(2L, 5L)),
      taperingEst = list(k = c(-2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L, 3L)),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 3L),
        var_est = "sample"
      ),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(2L, 5L)),
      taperingEst = list(k = c(3L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L, 3L)),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 3L),
        var_est = "sample"
      ),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c("a", 0.9)),
      thresholdingEst = list(gamma = c(0.2, "b")),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L, 3L)),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 3L),
        var_est = "sample"
      ),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, "b")),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L, 3L)),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 3L),
        var_est = "sample"
      ),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat[1:11, ],
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L, 3L)),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 3L),
        var_est = "sample"
      ),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = -0.1),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L, 3L)),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 3L),
        var_est = "sample"
      ),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = "a"),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L, 3L)),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 3L),
        var_est = "sample"
      ),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = 0.1),
      poetEst = list(lambda = -0.1, k = c(1L, 2L)),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 3L),
        var_est = "sample"
      ),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = 0.1),
      poetEst = list(lambda = 0.1, k = c(1.1, 2)),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 3L),
        var_est = "sample"
      ),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = 0.1),
      poetEst = list(lambda = 0.1, k = c(1, 2)),
      adaptiveLassoEst = list(lambda = -0.1, n = 0),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 3L),
        var_est = "sample"
      ),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = 0.1),
      poetEst = list(lambda = 0.1, k = c(1, 2)),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 3L),
        var_est = "sample"
      ),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = 0.1),
      poetEst = list(lambda = 0.1, k = c(1, 2)),
      robustPoetEst = list(
        lambda = 0.1, k = c(1L, 2L),
        var_est = c("samp", "huber")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = 0.1),
      poetEst = list(lambda = 0.1, k = c(1, 2)),
      robustPoetEst = list(
        lambda = "a", k = c(1L, 2L),
        var_est = c("mad", "huber")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = 0.1),
      poetEst = list(lambda = 0.1, k = c(1, 2)),
      robustPoetEst = list(
        lambda = 0.1, k = c(1.1, 2L),
        var_est = c("mad", "huber")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = 0.1),
      poetEst = list(lambda = 0.1, k = c(1, 2)),
      adaptiveLassoEst = list(lambda = "a", n = 0),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 3L),
        var_est = "sample"
      ),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = 0.1),
      poetEst = list(lambda = 0.1, k = c(1, 2)),
      adaptiveLassoEst = list(lambda = 0.5, n = 0),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 3L),
        var_est = "sample"
      ),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(-1, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = 0.1),
      poetEst = list(lambda = 0.1, k = c(1, 2)),
      adaptiveLassoEst = list(lambda = 0.5, n = 0),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 3L),
        var_est = "sample"
      ),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.1, 0.8), num_spikes = c(NULL, 1L, 2.5),
        noise = c(NULL, 0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = 0.1),
      poetEst = list(lambda = 0.1, k = c(1, 2)),
      adaptiveLassoEst = list(lambda = 0.5, n = 0),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(1L, 3L),
        var_est = "sample"
      ),
      spikedOperatorShrinkEst = list(
        p_n_ratio = c(0.1, 0.8), num_spikes = c(NULL, 1L),
        noise = c(NULL, -0.5)
      ),
      spikedFrobeniusShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      ),
      spikedSteinShrinkEst = list(
        p_n_ratio = c(0.5, 0.8), num_spikes = c(NULL, 1L, 2L),
        noise = c(NULL, 0.5)
      )
    ),
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
})

test_that("Only reasonable CV schemes pass checks", {
  expect_true(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_loss = cv_loss,
    cv_scheme = "v_fold",
    mc_split = mc_split,
    v_folds = v_folds,


    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_loss = cv_loss,
    cv_scheme = "loo",
    mc_split = mc_split,
    v_folds = v_folds,


    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_loss = cv_loss,
    cv_scheme = "mc",
    mc_split = 0,
    v_folds = v_folds,


    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_loss = cv_loss,
    cv_scheme = "mc",
    mc_split = 1,
    v_folds = v_folds,


    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_loss = cv_loss,
    cv_scheme = "mc",
    mc_split = -0.1,
    v_folds = v_folds,


    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_loss = cv_loss,
    cv_scheme = "mc",
    mc_split = 1.2,
    v_folds = v_folds,


    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_loss = cv_loss,
    cv_scheme = "v_fold",
    mc_split = mc_split,
    v_folds = 1,


    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_loss = cv_loss,
    cv_scheme = "v_fold",
    mc_split = mc_split,
    v_folds = nrow(dat),
    parallel = parallel
  ))
  expect_true(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_loss = cv_loss,
    cv_scheme = "v_fold",
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
})

test_that("Flag elements work as expected", {
  expect_true(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_loss = cv_loss,
    cv_scheme = "v_fold",
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = TRUE
  ))
  expect_true(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_loss = cv_loss,
    cv_scheme = "v_fold",
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = FALSE
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_loss = cv_loss,
    cv_scheme = "v_fold",
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = "TRUE"
  ))
})

test_that("checkArgs only allows pre-defined loss functions", {
  expect_true(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_loss = rlang::expr(cvMatrixFrobeniusLoss),
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_true(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_loss = rlang::expr(cvScaledMatrixFrobeniusLoss),
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_loss = rlang::expr(cvOperatorNorm),
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    parallel = parallel
  ))
})

test_that("checkArgs works well in cvCovEst Function", {
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, adaptiveLassoEst, spikedOperatorShrinkEst
    ),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(3L, 4L)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(3L, 4L),
        var_est = c("sample", "mad")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      spikedOperatorShrinkEst = list(p_n_ratio = c(0.5, 0.8))
    ),
    cv_loss = cvFrobeniusLoss, cv_scheme = "v_fold",
    mc_split = 0.5, v_folds = 5, parallel = FALSE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, adaptiveLassoEst, spikedOperatorShrinkEst
    ),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(3L, 4L)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(3L, 4L),
        var_est = c("sample", "mad")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      spikedOperatorShrinkEst = list(p_n_ratio = c(0.5, 0.8))
    ),
    cv_loss = cvMatrixFrobeniusLoss, cv_scheme = "mc",
    mc_split = 0.5, v_folds = 5, parallel = FALSE
  ))
  expect_silent(cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkEst, linearShrinkLWEst,
      thresholdingEst, sampleCovEst, bandingEst,
      taperingEst, nlShrinkLWEst, denseLinearShrinkEst,
      scadEst, poetEst, adaptiveLassoEst, spikedOperatorShrinkEst
    ),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L)),
      scadEst = list(lambda = c(0.1, 0.2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(3L, 4L)),
      robustPoetEst = list(
        lambda = c(0.1, 0.2), k = c(3L, 4L),
        var_est = c("sample", "mad")
      ),
      adaptiveLassoEst = list(lambda = c(0, 0.5), n = c(0, 0.5)),
      spikedOperatorShrinkEst = list(p_n_ratio = c(0.5, 0.8))
    ),
    cv_loss = cvScaledMatrixFrobeniusLoss, cv_scheme = "mc",
    mc_split = 0.5, v_folds = 5, parallel = FALSE
  ))
})
