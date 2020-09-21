test_that("SCAD Estimator regularizes as intended", {

  # hard thresholding
  expect_equal(scadThreshold(entry = 1, lambda = 1, a = 2), 0)
  expect_equal(scadThreshold(entry = -1, lambda = 1, a = 2), 0)
  expect_equal(scadThreshold(entry = 1, lambda = 0.5, a = 2), 0.5)
  expect_equal(scadThreshold(entry = -1, lambda = 0.5, a = 2), -0.5)

  # interpolation between hard and soft thresholding
  expect_equal(scadThreshold(entry = 1, lambda = 0.3, a = 4), 0.9)
  expect_equal(scadThreshold(entry = -1, lambda = 0.3, a = 4), -0.9)

  # return to hard thresholding
  expect_equal(scadThreshold(entry = 1.2, lambda = 0.3, a = 4), 1.2)
  expect_equal(scadThreshold(entry = -1.2, lambda = 0.3, a = 4), -1.2)
  expect_equal(scadThreshold(entry = 2, lambda = 0.3, a = 4), 2)
  expect_equal(scadThreshold(entry = -2, lambda = 0.3, a = 4), -2)
})

test_that("Adaptive Lasso Estimator regularizes as intended", {

  # values of lambda
  expect_equal(adaptiveLassoThreshold(entry = 1, lambda = 1, n = 1), 0)
  expect_equal(adaptiveLassoThreshold(entry = -1, lambda = 1, n = 1), 0)
  expect_equal(adaptiveLassoThreshold(entry = 1, lambda = 0.5, n = 1), 0.75)
  expect_equal(adaptiveLassoThreshold(entry = -1, lambda = 0.5, n = 1), -0.75)

  # values of n
  expect_equal(adaptiveLassoThreshold(entry = 1, lambda = 0.1, n = 1), 0.99)
  expect_equal(adaptiveLassoThreshold(entry = -1, lambda = 0.1, n = 1), -0.99)
  expect_equal(adaptiveLassoThreshold(entry = 1, lambda = 0.1, n = 2), 0.999)
  expect_equal(adaptiveLassoThreshold(entry = -1, lambda = 0.1, n = 2), -0.999)
})



