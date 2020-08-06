library(MASS)
set.seed(123)

# generate a 50x50 covariance matrix with unit variances and off-diagonal
# elements equal to 0.5
Sigma <- matrix(0.5, nrow = 50, ncol = 50) + diag(0.5, nrow = 50)

# sample 200 observations from multivariate normal with mean = 0, var = Sigma
dat <- mvrnorm(n = 200, mu = rep(0, 50), Sigma = Sigma)

# define the arguments as they appear inside cvCoveEst
estimators <- rlang::expr(c(linearShrinkEst, thresholdingEst, sampleCovEst,
                            linearShrinkLWEst, bandingEst, taperingEst,
                            nlShrinkLWEst))
estimator_params <- list(
  linearShrinkEst = list(alpha = c(0.1, 0.9)),
  thresholdingEst = list(gamma = c(0.2, 2)),
  bandingEst = list(k = c(1L, 5L)),
  taperingEst = list(k = c(2L, 6L))
)
cv_scheme <- "mc"
mc_split <- 0.5
v_folds <- 10
cv_loss <- rlang::expr(cvFrobeniusLoss)
boot_iter <- 25
center <- TRUE
scale <- FALSE
parallel <- FALSE

test_that("Only implmented estimators pass checks", {
  expect_true(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_true(checkArgs(
    dat = dat,
    estimators = rlang::expr(c(linearShrinkEst)),
    estimator_params = estimator_params,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = rlang::expr(c(linearShrinkEst, newEstimator)),
    estimator_params = estimator_params,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ), "Only estimators implemented in the cvCovEst package can be used.")
})

test_that("Only reasonable hyperparameters pass checks", {
  expect_true(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(-0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L))
    ),
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 1.1)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L))
    ),
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(-0.1, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L))
    ),
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1.1, 5L)),
      taperingEst = list(k = c(2L, 6L))
    ),
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(-2L, 5L)),
      taperingEst = list(k = c(2L, 6L))
    ),
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(2L, 5L)),
      taperingEst = list(k = c(2.1, 6L))
    ),
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(2L, 5L)),
      taperingEst = list(k = c(-2L, 6L))
    ),
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(2L, 5L)),
      taperingEst = list(k = c(3L, 6L))
    ),
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c("a", 0.9)),
      thresholdingEst = list(gamma = c(0.2, "b")),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L))
    ),
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, "b")),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L))
    ),
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat[1:11,],
    estimators = estimators,
    estimator_params = list(
      linearShrinkEst = list(alpha = c(0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L))
    ),
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
})

test_that("Only reasonable CV schemes pass checks",{
  expect_true(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_scheme = "v_fold",
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_scheme = "loo",
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_scheme = "mc",
    mc_split = 0,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_scheme = "mc",
    mc_split = 1,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_scheme = "mc",
    mc_split = -0.1,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_scheme = "mc",
    mc_split = 1.2,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_scheme = "v_fold",
    mc_split = mc_split,
    v_folds = 1,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_scheme = "v_fold",
    mc_split = mc_split,
    v_folds = nrow(dat),
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_true(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_scheme = "v_fold",
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = NULL,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_scheme = "v_fold",
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = rlang:expr(cvPenFrobeniusLoss),
    boot_iter = 9,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_scheme = "v_fold",
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = rlang:expr(cvPenFrobeniusLoss),
    boot_iter = NULL,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_scheme = "v_fold",
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = rlang:expr(cvPenFrobeniusLoss),
    boot_iter = 9,
    center = center,
    scale = scale,
    parallel = parallel
  ))
})

test_that("Choice of loss function is well defined", {
  expect_true(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_scheme = "v_fold",
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = rlang::expr(cvPenFrobeniusLoss),
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_true(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_scheme = "v_fold",
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = rlang::expr(cvFrobeniusLoss),
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_scheme = "v_fold",
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = rlang::expr(cvPenalizedFrobeniusLoss),
    boot_iter = boot_iter,
    center = center,
    scale = scale,
    parallel = parallel
  ))
})

test_that("Flag elements work as expected", {
  expect_true(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_scheme = "v_fold",
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = TRUE,
    scale = TRUE,
    parallel = TRUE
  ))
  expect_true(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_scheme = "v_fold",
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = FALSE,
    scale = FALSE,
    parallel = FALSE
  ))
  expect_error(checkArgs(
    dat = dat,
    estimators = estimators,
    estimator_params = estimator_params,
    cv_scheme = "v_fold",
    mc_split = mc_split,
    v_folds = v_folds,
    cv_loss = cv_loss,
    boot_iter = boot_iter,
    center = "TRUE",
    scale = "TRUE",
    parallel = "TRUE"
  ))
})

test_that("checkArgs works well in cvCovEstFunction", {
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
    cv_scheme = "v_fold", mc_split = 0.5, cv_loss = cvFrobeniusLoss,
    v_folds = 5, boot_iter = 10,
    center = TRUE, scale = FALSE, parallel = FALSE
  ))
  expect_error(cvCovEst(
    dat = dat,
    estimators = c(linearShrinkEst, linearShrinkLWEst,
                   thresholdingEst, sampleCovEst, bandingEst,
                   taperingEst, nlShrinkLWEst),
    estimator_params = list(
      linearShrinkEst = list(alpha = c(-0.1, 0.9)),
      thresholdingEst = list(gamma = c(0.2, 2)),
      bandingEst = list(k = c(1L, 5L)),
      taperingEst = list(k = c(2L, 6L))
    ),
    cv_scheme = "v_fold", mc_split = 0.5, cv_loss = cvPenFrobeniusLoss,
    v_folds = 5, boot_iter = 10,
    center = TRUE, scale = FALSE, parallel = FALSE
  ))
})
