test_that("Make sure that symmetric apply makes symmetric matrices", {
  cov_mat <- cov(mtcars)

  expect_equal(
    symmetricApply(dat = cov_mat, sym_fun = adaptiveLassoThreshold,
                   sym_args = c(lambda = 0.1, n = 0.1))[10, 1],
    symmetricApply(dat = cov_mat, sym_fun = adaptiveLassoThreshold,
                   sym_args = c(lambda = 0.1, n = 0.1))[1, 10])
  expect_equal(
    symmetricApply(dat = cov_mat, sym_fun = scadThreshold,
                   sym_args = c(a = 3.7, lambda = 0.1))[10, 1],
    symmetricApply(dat = cov_mat, sym_fun = scadThreshold,
                   sym_args = c(a = 3.7, lambda = 0.1))[1, 10])

  expect_equal(
    symmetricApply(dat = cov_mat, sym_fun = scadThreshold,
                   sym_args = c(a = 3.7, lambda = 0.1)),
    scadEst(dat = mtcars, lambda = 0.1)
  )
})
