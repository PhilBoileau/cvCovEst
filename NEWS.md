# cvCovEst 0.0.15 (2020-10-03)

* Fixing a minor dimension error in `nlShrinkLWEst` by changing a conditional,
  as per https://github.com/PhilBoileau/cvCovEst/issues/23.
* Enforcing the [`tidyverse` code style](https://style.tidyverse.org/) via the
  first call to `styler` in this codebase (via `make style`).
* Enforcing 80-columns in `NEWS.md`.

# cvCovEst 0.0.14 (2020-09-29)

* Replacing `stats::cov` with `coop::covar` after resolving the issue on Linux
  machines, as per https://github.com/PhilBoileau/cvCovEst/issues/18.
* Removed calculation of spurious risk ratios from `cvCovEst` and from
  `cvFrobeniusLoss` when the true covariance matrix is passed in.

# cvCovEst 0.0.13 (2020-09-28)

* Changing loss function computation so that it is more computationally
  efficient.

# cvCovEst 0.0.12 (2020-09-26)

* Removing `coop::covar` due to strange parallelization issue on Linux machines.
  Hopefully we can use it again one day.
* Prevent users from including a lone estimator as input to `cvCovEst` if the
  estimator in questions doesn't have any hyperparameters.
* Coerce sparse, true covariance matrices to regular matrix objects if and when
  input to `cvCovEst`.

# cvCovEst 0.0.11 (2020-09-24)

* Adding additional risk difference ratio calculations when the true covariance
  matrix of Gaussian Multivariate data is provided.

# cvCovEst 0.0.10 (2020-09-22)

* Added adaptive LASSO estimator.
* Users now have the option to include the true covariance matrix of their
  multivariate Gaussian data, allowing them to compare cvCovEst's selection
  versus that of the cross-valdidated oracle.

# cvCovEst 0.0.9 (2020-09-13)

* Adding POET estimator
* Estimators can now take multiple hyperparameter arguments.

# cvCovEst 0.0.8 (2020-09-13)

* Adding smoothly clipped absolute deviation thresholding estimator.

# cvCovEst 0.0.7 (2020-08-22)

* Updated the loss computation; it now patches the formula used in the draft.
  Note that it vastly overestimates the true risk of an estimator, but that it
  provides an equivalent decision rule compared to a matrix-based loss. Perhaps
  we're missing a scaling factor in our calculations?
* Moved Frobenius loss calculations to cv fold loss function.
* Removed the penalized cross-validation loss. Doesn't make sense to include.
* Included check for centered data matrix.

# cvCovEst 0.0.6 (2020-08-19)

* Adding dense covariance matrix linear shrinkage estimator.
* Updating citations in estimators docs.

# cvCovEst 0.0.5 (2020-08-06)

* Adding analytical nonlinear shrinkage estimator.

# cvCovEst 0.0.4 (2020-07-05)

* Adding tapering estimator.

# cvCovEst 0.0.4 (2020-06-23)

* Adding argument checker for `cvCovEst` function.

# cvCovEst 0.0.3 (2020-06-22)

* Adding banding estimator.

# cvCovEst 0.0.2 (2020-06-06)

* Added unpenalized frobenius matrix loss.
* cvCovEst() now requires a vector of candidate estimator functions be passed to
  the estimators argument, instead of a vector of characters corresponding to
  these candidates' names.

# cvCovEst 0.0.1 (2020-05-21)

* Minor changes to core routines, including changes to use of `origami`.
* Updates to documentation, including `Roxygen` styling.
* Addition of templates for vignette and JOSS paper.

# cvCovEst 0.0.0.9000 (2020-05-01)

* Added a `NEWS.md` file to track changes to the package.
