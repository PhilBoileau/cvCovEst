# cvCovEst 0.3.7 (2021-07-03)

+ Calling `summary.cvCovEst()` when a single summary function is
  specified now immediately returns a table instead of a list of length 1 that
  contains said table.
+ Tables returned by `summary.cvCovEst()` no longer have `dplyr` groups.
+ Fixed typo in Toy Dataset Example section of paper.

# cvCovEst 0.3.6 (2021-06-19)

+ Renamed `empirical_risk` column in `risk_df` table output by `cvCovEst()` to
  `cv_risk`.
+ Added additional citations of existing R packages for covariance matrix
  estimation in our JOSS submission.
+ Added more comprehensive tests for the available loss functions.

# cvCovEst 0.3.5 (2021-02-24)

+ Setting 'LazyLoad' to 'false' in DESCRIPTION to address CRAN checks notes.

# cvCovEst 0.3.4 (2021-02-24)

+ Fixing plotting labels and table column names, along with associated
  documentation.
+ Fixing links to pass CRAN checks
+ Reducing size of toy datasets to increase testing speed

# cvCovEst 0.3.3 (2021-02-23)

+ Adding note to `robustPoetEst()` warning again its use for correlation matrix
  estimation.
+ Fixing bug in `robustPoetEst` plots.

# cvCovEst 0.3.2 (2021-02-22)

+ Adding preprint citation information.

# cvCovEst 0.3.1 (2021-02-13)

+ Edited documentation to meet CRAN specifications.

# cvCovEst 0.3.0 (2021-02-10)

+ Added information and simulated data examples of plotting and summary
  functions.
+ Made the `cvCovEst` R package a public repository GitHub. 

# cvCovEst 0.2.0 (2021-02-02)

+ `cvCovEst` now possesses a slew of diagnostic and visualization tools. A
  detailed description of these functions will be added to the vignette in the
  near future.

# cvCovEst 0.1.7 (2021-01-29)

+ Minor clarifying updates to the documentation and the vignette.
+ Updates to `NEWS.md`, adding consistency in bullet point indicator and
  enforcing the 80-column rule.
+ Tweaks to dependencies, removing reliance on `stringi` since only invoked in
  a single pipe call in `checkArgs`.

# cvCovEst 0.1.6 (2021-01-24)

+ Added basic examples to all exported functions.

# cvCovEst 0.1.5 (2021-01-23)

+ Made `cvMatrixFrobeniusLoss` the default loss function.

# cvCovEst 0.1.4 (2020-12-29)

+ Added `cvScaledMatrixFrobeniusLoss`, a new matrix-based loss function that
  scales squared error calculations associated with each entry of a covariance
  matrix estimate by the sample variances of the entry's row and column
  variables. This is particularly useful if the features of your dataset are of
  different magnitude. It's approximately equivalent to estimating the
  correlation matrix, but without the need to re-scale the estimated
  correlation matrix to be an estimated covariance matrix.

# cvCovEst 0.1.3 (2020-12-29)

+ Fixed error with `denseLinearShrinkEst`: the shrinkage parameter was often
  selected such that the dense target was returned as the estimate.

# cvCovEst 0.1.2 (2020-12-21)

+ Completed vignette.

# cvCovEst 0.1.1 (2020-12-06)

+ `robustPoetEst` has been added to the library of candidate estimators.

# cvCovEst 0.1.0 (2020-11-19)

+ cvCovEst version 0.1.0 is used in the accompanying manuscript to generate all
  results.
+ It is stored as a separate branch called 'preprint'.

# cvCovEst 0.0.18 (2020-10-19)

+ cvCovEst now accepts cvMatrixFrobeniusLoss as a loss function. This loss
  function is a matrix-based alternative to the standard loss function. Through
  Proposition 1 of the method's manuscript the resulting selections of each
  loss are identical for any fixed cross-validation scheme. However, the
  matrix-based loss is more computationally efficient. Other minor tweaks to
  testing procedures.

# cvCovEst 0.0.17 (2020-10-19)

+ cvCovEst can now be run in parallel using future.

# cvCovEst 0.0.16 (2020-10-16)

+ When provided with the true covariance matrix, cvCovEst now outputs the
  conditional cross-validated risk differences of the cross-validation
  selection and the oracle selections.

# cvCovEst 0.0.15 (2020-10-03)

+ Fixing a minor dimension error in `nlShrinkLWEst` by changing a conditional,
  as per https://github.com/PhilBoileau/cvCovEst/issues/23.
+ Enforcing the [`tidyverse` code style](https://style.tidyverse.org/) via the
  first call to `styler` in this codebase (via `make style`).
+ Enforcing 80-columns in `NEWS.md`.

# cvCovEst 0.0.14 (2020-09-29)

+ Replacing `stats::cov` with `coop::covar` after resolving the issue on Linux
  machines, as per https://github.com/PhilBoileau/cvCovEst/issues/18.
+ Removed calculation of spurious risk ratios from `cvCovEst` and from
  `cvFrobeniusLoss` when the true covariance matrix is passed in.

# cvCovEst 0.0.13 (2020-09-28)

+ Changing loss function computation so that it is more computationally
  efficient.

# cvCovEst 0.0.12 (2020-09-26)

+ Removing `coop::covar` due to strange parallelization issue on Linux machines.
  Hopefully we can use it again one day.
+ Prevent users from including a lone estimator as input to `cvCovEst` if the
  estimator in questions doesn't have any hyperparameters.
+ Coerce sparse, true covariance matrices to regular matrix objects if and when
  input to `cvCovEst`.

# cvCovEst 0.0.11 (2020-09-24)

+ Adding additional risk difference ratio calculations when the true covariance
  matrix of Gaussian Multivariate data is provided.

# cvCovEst 0.0.10 (2020-09-22)

+ Added adaptive LASSO estimator.
+ Users now have the option to include the true covariance matrix of their
  multivariate Gaussian data, allowing them to compare cvCovEst's selection
  versus that of the cross-validated oracle.

# cvCovEst 0.0.9 (2020-09-13)

+ Adding POET estimator
+ Estimators can now take multiple hyperparameter arguments.

# cvCovEst 0.0.8 (2020-09-13)

+ Adding smoothly clipped absolute deviation thresholding estimator.

# cvCovEst 0.0.7 (2020-08-22)

+ Updated the loss computation; it now patches the formula used in the draft.
  Note that it vastly overestimates the true risk of an estimator, but that it
  provides an equivalent decision rule compared to a matrix-based loss. Perhaps
  we're missing a scaling factor in our calculations?
+ Moved Frobenius loss calculations to cv fold loss function.
+ Removed the penalized cross-validation loss. Doesn't make sense to include.
+ Included check for centered data matrix.

# cvCovEst 0.0.6 (2020-08-19)

+ Adding dense covariance matrix linear shrinkage estimator.
+ Updating citations in estimators docs.

# cvCovEst 0.0.5 (2020-08-06)

+ Adding analytical nonlinear shrinkage estimator.

# cvCovEst 0.0.4 (2020-07-05)

+ Adding tapering estimator.

# cvCovEst 0.0.4 (2020-06-23)

+ Adding argument checker for `cvCovEst` function.

# cvCovEst 0.0.3 (2020-06-22)

+ Adding banding estimator.

# cvCovEst 0.0.2 (2020-06-06)

+ Added unpenalized frobenius matrix loss.
+ cvCovEst() now requires a vector of candidate estimator functions be passed
  to he estimators argument, instead of a vector of characters corresponding to
  these candidates' names.

# cvCovEst 0.0.1 (2020-05-21)

+ Minor changes to core routines, including changes to use of `origami`.
+ Updates to documentation, including `Roxygen` styling.
+ Addition of templates for vignette and JOSS paper.

# cvCovEst 0.0.0.9000 (2020-05-01)

+ Added a `NEWS.md` file to track changes to the package.
