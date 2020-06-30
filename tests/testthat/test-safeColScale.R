context("Test routines for safe column scaling")

# initialization of inputs for testing function
library(MASS)
set.seed(123)

# sample 200 observations from multivariate normal with mean = 0, var = Sigma
Sigma <- matrix(0.5, nrow = 50, ncol = 50) + diag(0.5, nrow = 50)
dat <- mvrnorm(n = 200, mu = rep(0, 50), Sigma = Sigma)

center <- TRUE
scale <- TRUE

# unit tests
test_that("`scale` and `safeColScale` produce the same output", {
  dat_safeColScaled <- safeColScale(dat, center, scale)
  dat_scaled <- scale(dat, center, scale)
  expect_equal(dat_scaled, dat_safeColScaled)
})

test_that("`safeColScale` avoid NA in its output even when `scale` fails to", {
  # make first column constant
  dat[, 1] <- 0.5

  # check that scale fails to avoid NAs
  dat_scaled <- scale(dat, center, scale)
  expect_true(sum(colSums(is.na(dat_scaled))) > 0)

  # check that safeColScale avoids NAs
  dat_safeColScaled <- safeColScale(dat, center, scale)
  expect_true(sum(colSums(is.na(dat_safeColScaled))) == 0)
})

# `safeColScale` is faster than `scale`
# library(ggplot2)
# library(microbenchmark)
# mb_scale <- microbenchmark(safeColScale(target), scale(target),
# times = 20, unit = "s")
# p_mb_scale <- ggplot(data = mb_scale, aes(y = time / 1e9, x = expr)) +
# geom_violin() + theme_grey(base_size = 20) +
# xlab("Method") + ylab("Time (seconds)")
# print(mb_scale)
# print(p_mb_scale)
