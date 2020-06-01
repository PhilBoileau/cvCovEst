
# R/`cvCovEst`

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/PhilBoileau/cvCovEst.svg?branch=master)](https://travis-ci.com/PhilBoileau/cvCovEst)
<!-- badges: end -->

> Cross-Validated Covariance Matrix Estimation

**Authors:** [Philippe Boileau](https://pboileau.ca) and [Nima
Hejazi](https://nimahejazi.org)

-----

## What’s `cvCovEst`?

`cvCovEst` implements an efficient cross-validated procedure for
covariance matrix estimation, particularly useful in high-dimensional
settings. The general methodology allows for cross-validation to be used
to data adaptively identify the optimal estimator of the covariance
matrix from a pre-specified set of candidate estimators. Also included
are diagnostic tools for assessing the quality of covariance matrix
estimators.

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

# sample 200 observations from multivariate normal with mean = 0, var = Sigma
dat <- mvrnorm(n = 200, mu = rep(0, 50), Sigma = Sigma)

# run CV-selector
cv_cov_est_out <- cvCovEst(
    dat = dat,
    estimators = c("linearShrinkEst", "thresholdingEst"),
    estimator_params = list("linearShrinkEst" = list("alpha" = c(0.1, 0.9)),
                            "thresholdingEst" = list("gamma" = c(0.2, 2))),
    cv_scheme = "v_fold", mc_split = 0.5,
    v_folds = 5, boot_iter = 20,
    center = TRUE, scale = FALSE, parallel = FALSE
  )

# print the table of risk estimates
# NOTE: the estimated covariance matrix is accessible via the `$estimate` slot
cv_cov_est_out$risk_df
#> # A tibble: 4 x 3
#> # Groups:   estimator [2]
#>   estimator       hyperparameters empirical_risk
#>   <chr>           <chr>                    <dbl>
#> 1 thresholdingEst gamma = 0.2               2.49
#> 2 linearShrinkEst alpha = 0.9               2.49
#> 3 linearShrinkEst alpha = 0.1              11.3 
#> 4 thresholdingEst gamma = 2                14.5
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

© 2020 [Philippe Boileau](https://pboileau.ca)

The contents of this repository are distributed under the MIT license.
See file
[`LICENSE.md`](https://github.com/PhilBoileau/cvCovEst/blob/master/LICENSE.md)
for details.
