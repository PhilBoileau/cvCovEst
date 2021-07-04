library(future.apply)
library(MASS)
set.seed(123)

# generate a 50x50 covariance matrix with unit variances and off-diagonal
# elements equal to 0.5
Sigma <- matrix(0.5, nrow = 50, ncol = 50) + diag(0.5, nrow = 50)

# sample 200 observations from multivariate normal with mean = 0, var = Sigma
dat <- mvrnorm(n = 200, mu = rep(0, 50), Sigma = Sigma)

estimator_params <- list(
  poetEst = list( # 3 x 3 different indexed hyperparameters
    lambda = c(0.01, 0.05, 0.1),
    k = c(1L, 3L, 5L)
  ),
  adaptiveLassoEst = list( # 2 x 1 indexed hyperparameters
    lambda = c(0.01, 0.1),
    n = c(2L, 3L)
  ),
  linearShrinkEst = list( # 2 different indexed hyperparameters
    alpha = c(0.1, 0.9)
  ),
  bandingEst = list( # 1 indexed hyperparameter
    k = 2L
  )
)

# All with hypers
cvTestH <- cvCovEst(
  dat = dat,
  estimators = c(
    poetEst,
    adaptiveLassoEst,
    linearShrinkEst,
    bandingEst
  ),
  estimator_params = estimator_params,
  cv_scheme = "v_fold",
  cv_loss = cvMatrixFrobeniusLoss,
  v_folds = 5,
  parallel = FALSE,
  center = TRUE,
  scale = TRUE
)

# With one no hyper
cvTestNH <- cvCovEst(
  dat = dat,
  estimators = c(
    nlShrinkLWEst,
    bandingEst
  ),
  estimator_params = estimator_params,
  cv_scheme = "v_fold",
  cv_loss = cvMatrixFrobeniusLoss,
  v_folds = 5,
  parallel = FALSE,
  center = TRUE,
  scale = TRUE
)

has_hypers <- c(
  "linearShrinkEst", "thresholdingEst", "bandingEst", "taperingEst",
  "scadEst", "poetEst", "robustPoetEst", "adaptiveLassoEst"
)

# Class Test
test_that("Objects of other known classes throw an error", {
  # cvCovest class
  expect_s3_class(cvTestH, "cvCovEst")
  expect_silent(
    summary(cvTestH)
  )
  expect_silent(
    cvTestH %>% summary()
  )
  # different class
  class(cvTestH) <- "lm"
  expect_error(
    summary(cvTestH)
  )
  expect_error(
    cvTestH %>% summary()
  )
  # other object disguised as cvCovEst object
  disguise <- c("disguise")
  class(disguise) <- "cvCovEst"
  expect_error(
    summary(disguise)
  )
})

test_that("Only current implemented summary functions are allowed", {
  expect_silent(
    summary(cvTestH, summ_fun = "bestInClass")
  )
  expect_error(
    summary(cvTestH, summ_fun = "other")
  )
})

test_that("Only supported summary statistics are allowed for plotting", {
  expect_silent(
    cvMultiMelt(
      dat = cvTestH,
      estimator = c("poetEst"),
      stat = c("min"),
      dat_orig = dat,
      cv_details = "",
      has_hypers = has_hypers
    )
  )
  expect_error(
    suppressWarnings(
      cvMultiMelt(
        dat = cvTestH,
        estimator = c("poetEst"),
        stat = c("mean"),
        dat_orig = dat,
        cv_details = "",
        has_hypers = has_hypers
      )
    )
  )
})

test_that("Valid estimator arguments are passed to plotting functions", {
  # Non-cvCovEst estimator
  expect_error(
    cvMultiMelt(
      dat = cvTestH,
      estimator = c("linearShrinkEst", "other"),
      stat = c("min"),
      dat_orig = dat,
      cv_details = "",
      has_hypers = has_hypers
    )
  )
  # Estimator not originally called to cvCovEst()
  expect_error(
    cvMultiMelt(
      dat = cvTestH,
      estimator = c("poetEst", "scadEst"),
      stat = c("min"),
      dat_orig = dat,
      cv_details = "",
      has_hypers = has_hypers
    )
  )
  # Multiple plots of the same estimator
  expect_error(
    cvMultiMelt(
      dat = cvTestH,
      estimator = c("nlShrinkLWEst"),
      stat = c("min", "max"),
      dat_orig = dat,
      cv_details = "",
      has_hypers = has_hypers
    )
  )
})

test_that("Indexing by only 1 hyperparameter throws an error in risk plot", {
  expect_error(
    plot.cvCovEst(
      x = cvTestH,
      dat_orig = dat,
      estimator = "bandingEst",
      plot_type = "risk"
    )
  )
})

test_that("Calling risk plot for non-hyper estimator throws an error", {
  expect_error(
    plot.cvCovEst(
      x = cvTestNH,
      dat_orig = dat,
      estimator = "nlShrinkLWEst",
      plot_type = "risk"
    )
  )
})

test_that("Calling for multiple stats for non-hyper estimator gets a message", {
  expect_message(
    plot.cvCovEst(
      x = cvTestNH,
      dat_orig = dat,
      estimator = "nlShrinkLWEst",
      plot_type = "eigen",
      k = 50,
      stat = c("min", "max")
    )
  )
})

test_that("Plotting only works if estimator was passed to cvCovEst", {
  expect_error(
    plot.cvCovEst(
      x = cvTestH,
      dat_orig = dat,
      estimator = "nlShrinkLWEst",
      plot_type = "eigen",
      k = 50,
      stat = c("min")
    )
  )
})

test_that("Asking for more k than exist throws an error", {
  expect_error(
    plot.cvCovEst(
      x = cvTestH,
      dat_orig = dat,
      estimator = "linearShrinkEst",
      plot_type = "eigen",
      k = 51,
      stat = c("min")
    )
  )
})

test_that("Plot method throws other errors where appropriate", {
  expect_message(
    plot.cvCovEst(
      x = cvTestH,
      dat_orig = dat,
      estimator = c("linearShrinkEst"),
      plot_type = ("summary")
    )
  )
})
