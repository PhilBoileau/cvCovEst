
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`cvCovEst`

> Cross-Validated Covariance Matrix Estimation

**Authors:** [Philippe Boileau](https://pboileau.ca) and [Nima
Hejazi](https://nimahejazi.org)

-----

# TODO

  - \[ \] Setup repository `README.Rmd` file
      - \[X\] Draft description (for `DESCRIPTION` file also)
      - \[X\] Installation instructions
      - \[X\] Issues sections
      - \[X\] Contributions
      - \[X\] License
  - \[ \] Identify classes of matrix objects we want to work with
  - \[ \] Create loss functions
  - \[ \] Identify estimators we want to consider. Determine if they are
    already implemented

-----

## Description

`cvCovEst` implements an efficient cross-validated approach to
covariance matrix estimation in high-dimensional settings. This
procedure data-adaptively identifies the optimal estimator of the
covariance matrix from a set of candidates. Dignostic tools are also
provided.

-----

## Installation

The `cvCovEst` package can installed from GitHub via `remotes`:

    remotes::install_github("PhilBoileau/scPCA")

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

© 2019-2020 [Philippe Boileau](https://pboileau.ca/)

The contents of this repository are distributed under the MIT license.
See file `LICENSE` for details.
