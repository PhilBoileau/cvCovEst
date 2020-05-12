# generate a 50x50 covariance matrix with variances = 1 and off diagonal elements
# equal to 0.5
Sigma <- matrix(0.5, nrow = 50, ncol = 50) + diag(0.5, nrow = 50)

# sample 200 observations from multivariate normal with mean = 0 and var = Sigma
library(MASS)
set.seed(123)
dat <- mvrnorm(n = 200, mu = rep(0, 50), Sigma = Sigma)
