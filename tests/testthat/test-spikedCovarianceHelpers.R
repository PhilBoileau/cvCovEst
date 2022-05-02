test_that(paste("Value near 1 is returned when medians of eigenvalues and",
                "Marcenko Pastur distributions are identical"), {
  eigenvalues <- rep(1, 500)
  p_n_ratio <- 0.8
  noise_estimate <- estimateNoise(eigenvalues, p_n_ratio)
  expect_equal(noise_estimate, 1, tolerance = 0.5)
})

test_that(paste("Three scaled eigenvalues are returned when num_spikes equals",
                "three"), {
  eig_vals <- c(15, 12, 9, 3, 3, 3)
  noise <- 3
  p_n_ratio <- 0.8
  num_spikes <- 3
  scaled_eig_vals <- scaleEigVals(eig_vals, noise, p_n_ratio, num_spikes)
  expect_equal(scaled_eig_vals, c(5, 4, 3))
})

test_that(paste("Two scaled eigenvalues are returned when num_spikes is not",
                "defined"), {
  eig_vals <- c(15, 12, 9, 3, 3, 3)
  noise <- 3
  p_n_ratio <- 0.8
  scaled_eig_vals <- scaleEigVals(eig_vals, noise, p_n_ratio, NULL)
  # cutt off: (1 + sqrt(0.8))^2 = 3.589
  expect_equal(scaled_eig_vals, c(5, 4))
})

test_that("ell is a vector of ones when the number of spikes equals zero", {
  scaled_eig_vals <- c()
  p <- 10
  p_n_ratio <- 0.5
  ell <- computeEll(scaled_eig_vals, p, p_n_ratio)
  expect_equal(ell, rep(1, 10))
})

test_that(paste("ell is a vector of a scaled eigenvalue and ones when the",
                "number of spikes equals one"), {
  scaled_eig_vals <- c(10)
  p <- 10
  p_n_ratio <- 0.5
  ell <- computeEll(scaled_eig_vals, p, p_n_ratio)
  shrunken_spike <- ((10 + 1 - 0.5) + sqrt((10 + 1 - 0.5)^2 - 4*10)) / 2
  expect_equal(ell, c(shrunken_spike, rep(1, 9)))
})

test_that("c equals has one non-zero entry when ell has one spiked value", {
  scaled_eig_vals <- c(10)
  p <- 10
  p_n_ratio <- 0.5
  ell <- computeEll(scaled_eig_vals, p, p_n_ratio)
  shrunken_spike <- ((10 + 1 - 0.5) + sqrt((10 + 1 - 0.5)^2 - 4*10)) / 2
  c_donoho <- computeC(ell, p_n_ratio)
  leading_c <- sqrt((1 - 0.5 / (shrunken_spike - 1)^2) /
                    (1 + 0.5 / (shrunken_spike - 1)))
  expect_equal(c_donoho[2:10], rep(0, 9))
  expect_equal(c_donoho[1], leading_c)
})

