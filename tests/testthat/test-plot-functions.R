library(future.apply)
library(MASS)
set.seed(123)

# generate a 50x50 covariance matrix with unit variances and off-diagonal
# elements equal to 0.5
Sigma <- matrix(0.5, nrow = 50, ncol = 50) + diag(0.5, nrow = 50)

# sample 200 observations from multivariate normal with mean = 0, var = Sigma
dat <- mvrnorm(n = 200, mu = rep(0, 50), Sigma = Sigma)

estimator_params <- list(
  poetEst = list( # 5 x 5 different indexed hyperparameters
    lambda = seq(0.2, 1, by = 0.2),
    k = as.integer(seq(1, 5, by = 1))),
  linearShrinkEst = list( # 2 different indexed hyperparameters
    alpha = seq(0, 1, by = 0.5)),
  bandingEst = list( # 1 indexed hyperparameter
    k = 2L)
)

cvTest <- cvCovEst(
  dat = dat,
  estimators = c(
    poetEst,
    linearShrinkEst,
    bandingEst,
    nlShrinkLWEst # No hyperparameters
    ),
  estimator_params = estimator_params,
  cv_scheme = "v_fold",
  cv_loss = cvMatrixFrobeniusLoss,
  v_folds = 10,
  parallel = FALSE,
  center = TRUE,
  scale = FALSE
)

# Class Test
test_that("Objects of other known classes throw an error", {
  # cvCovest class
  expect_s3_class(cvTest, "cvCovEst")
  expect_silent(
    summary(cvTest)
    )
  expect_silent(
    cvTest %>% summary()
  )
  # different class
  class(cvTest) <- "lm"
  expect_error(
    summary(cvTest)
    )
  expect_error(
    cvTest %>% summary()
  )
  # other object disguised as cvCovEst object
  disguise <- c('disguise')
  class(disguise) <- "cvCovEst"
  expect_error(
    summary(disguise)
  )
})

test_that("Only current implemented summary functions are allowed", {
  expect_silent(
    summary(cvTest, summ_fun = 'bestInClass')
  )
  expect_error(
    summary(cvTest, summ_fun = 'other')
  )
})

test_that("Only supported summary statistics are allowed for plotting", {
  expect_silent(
    cvMultiMelt(
      dat = cvTest,
      estimator = c('poetEst'),
      stat = c('min'),
      dat_orig = dat)
    )
  expect_error(
    cvMultiMelt(
      dat = cvTest,
      estimator = c('poetEst'),
      stat = c('mean'),
      dat_orig = dat)
  )

} )

test_that("Valid estimator arguments are passed to plotting functions",  {
  # Non-cvCovEst estimator
  expect_error(
    cvMultiMelt(
      dat = cvTest,
      estimator = c(
        'linearShrinkEst',
        'other'),
      stat = c('min'),
      dat_orig = dat)
  )
  # Estimator not originally called to cvCovEst()
  expect_error(
    cvMultiMelt(
      dat = cvTest,
      estimator = c('poetEst',
                    'scadEst'),
      stat = c('min'),
      dat_orig = dat)
  )
  # Multiple plots of the same estimator
  expect_error(
    cvMultiMelt(
      dat = cvTest,
      estimator = c('nlShrinkLWEst'),
      stat = c('min', 'max'),
      dat_orig = dat)
  )
})


