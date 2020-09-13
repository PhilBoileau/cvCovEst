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
