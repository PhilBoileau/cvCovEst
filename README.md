
# R/`cvCovEst`

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/PhilBoileau/cvCovEst.svg?token=YL3L6rYQtszHibWHgFxU&branch=master)](https://travis-ci.com/PhilBoileau/cvCovEst)
[![codecov](https://codecov.io/gh/PhilBoileau/cvCovEst/branch/master/graph/badge.svg?token=miHiqpGXxJ)](https://codecov.io/gh/PhilBoileau/cvCovEst)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

> Cross-Validated Covariance Matrix Estimation

**Authors:** [Philippe Boileau](https://pboileau.ca), [Brian
Collica](https://www.linkedin.com/in/brian-collica-553b0b94), and [Nima
Hejazi](https://nimahejazi.org)

-----

## What’s `cvCovEst`?

`cvCovEst` implements an efficient cross-validated procedure for
covariance matrix estimation, particularly useful in high-dimensional
settings. The general methodology allows for cross-validation to be used
to data adaptively identify the optimal estimator of the covariance
matrix from a pre-specified set of candidate estimators. An overview of
of the framework is provided in the package vignette.

<!-- A suite of plotting and diagnostic tools are also included. -->

-----

## Installation

Install the *development version* of the `cvCovEst` package from GitHub
via [`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("PhilBoileau/cvCovEst")
```

-----

## Example

To illustrate how `cvCovEst` may be used to select an optimal covariance
matrix estimator via cross-validation, consider the following toy
example:

``` r
library(MASS)
library(cvCovEst)
set.seed(1584)

# generate a 50x50 covariance matrix with unit variances and off-diagonal
# elements equal to 0.5
Sigma <- matrix(0.5, nrow = 50, ncol = 50) + diag(0.5, nrow = 50)

# sample 50 observations from multivariate normal with mean = 0, var = Sigma
dat <- mvrnorm(n = 50, mu = rep(0, 50), Sigma = Sigma)

# run CV-selector
cv_cov_est_out <- cvCovEst(
    dat = dat,
    estimators = c(linearShrinkLWEst, denseLinearShrinkEst,
                   thresholdingEst, poetEst, sampleCovEst),
    estimator_params = list(
      thresholdingEst = list(gamma = c(0.2, 2)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L))
    ),
    cv_loss = cvMatrixFrobeniusLoss,
    cv_scheme = "v_fold",
    v_folds = 5,
  )

# print the table of risk estimates
# NOTE: the estimated covariance matrix is accessible via the `$estimate` slot
cv_cov_est_out$risk_df
#> # A tibble: 9 x 3
#>   estimator            hyperparameters      empirical_risk
#>   <chr>                <chr>                         <dbl>
#> 1 linearShrinkLWEst    hyperparameters = NA           357.
#> 2 poetEst              lambda = 0.2, k = 1            369.
#> 3 poetEst              lambda = 0.2, k = 2            371.
#> 4 poetEst              lambda = 0.1, k = 1            377.
#> 5 poetEst              lambda = 0.1, k = 2            377.
#> 6 denseLinearShrinkEst hyperparameters = NA           379.
#> 7 sampleCovEst         hyperparameters = NA           379.
#> 8 thresholdingEst      gamma = 0.2                    384.
#> 9 thresholdingEst      gamma = 2                      827.
```

-----

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/PhilBoileau/cvCovEst/issues).

-----

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/PhilBoileau/cvCovEst/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

-----

## License

© 2020-2021 [Philippe Boileau](https://pboileau.ca)

The contents of this repository are distributed under the MIT license.
See file
[`LICENSE.md`](https://github.com/PhilBoileau/cvCovEst/blob/master/LICENSE.md)
for details.
