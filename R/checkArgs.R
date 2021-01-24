#' Check Arguments Passed to cvCovEst
#'
#' @description The \code{checkArgs} function verifies that all arguments
#'  passed to \code{\link{cvCovEst}} function meet its specifications.
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#' @param estimators A \code{list} of estimator functions to be
#'  considered in the cross-validated selection procedure.
#' @param estimator_params A named \code{list} of arguments corresponding to
#'  the hyperparameters of covariance matrix estimators in \code{estimators}.
#'  The name of each list element should match the name of an estimator passed
#'  to \code{estimators}. Each element of the \code{estimator_params} is itself
#'  a named \code{list}, with the names corresponding to a given estimators'
#'  hyperparameter(s). These hyperparameters may be in the form of a single
#'  \code{numeric} or a \code{numeric} vector. If no hyperparameter is needed
#'  for a given estimator, then the estimator need not be listed.
#' @param cv_loss A \code{function} indicating the loss function to use.
#'  Defaults to the Frobenius loss, \code{\link{cvMatrixFrobeniusLoss}}.
#'  An observation-based version, \code{\link{cvFrobeniusLoss}}, is also
#'  available. Finally, the \code{cvScaledMatrixFrobeniusLoss} is made available
#'  for situations where \code{dat}'s variables of different scales.
#' @param cv_scheme A \code{character} indicating the cross-validation scheme
#'  to be employed. There are two options: (1) V-fold cross-validation, via
#'  \code{"v_folds"}; and (2) Monte Carlo cross-validation, via \code{"mc"}.
#'  Defaults to Monte Carlo cross-validation.
#' @param mc_split A \code{numeric} between 0 and 1 indicating the proportion
#'  of data in the validation set of each Monte Carlo cross-validation fold.
#' @param v_folds A \code{integer} larger than or equal to 1 indicating the
#'  number of folds to use during cross-validation. The default is 10,
#'  regardless of cross-validation scheme.
#' @param center A \code{logical} indicating whether or not to center the
#'  columns of \code{dat}.
#' @param scale A \code{logical} indicating whether or not to scale the
#'  columns of \code{dat} to have variance 1.
#' @param parallel A \code{logical} option indicating whether to run the main
#'  cross-validation loop with \code{\link[future.apply]{future_lapply}}. This
#'  is passed directly to \code{\link[origami]{cross_validate}}.
#' @param true_cov_mat A \code{matrix} like object representing the true
#'  covariance matrix of the data generating distribution, which is assumed to
#'  be Gaussian. This parameter is intended for use only in simulation studies,
#'  and defaults to a value of \code{NULL}. If not \code{NULL}, various
#'  conditional risk difference ratios of the estimator selected
#'  by \code{\link{cvCovEst}} are computed relative to the different oracle
#'  selectors. NOTE: This parameter will be phased out by the release of version
#'  1.0.0.
#'
#' @importFrom assertthat assert_that is.flag
#' @importFrom methods is
#' @importFrom dplyr case_when "%>%"
#' @importFrom rlang is_bare_numeric is_integer is_null expr_deparse
#' @importFrom stringr str_sub str_split
#' @importFrom stringi stri_remove_empty
#'
#' @return Whether all argument conditions are satisfied
#'
#' @keywords internal
checkArgs <- function(dat,
                      estimators, estimator_params,
                      cv_loss, cv_scheme, mc_split, v_folds,
                      center, scale, parallel,
                      true_cov_mat = NULL) {

  # turn list of estimator functions to a vector of strings
  estimators <- estimators %>%
    rlang::expr_deparse() %>%
    stringr::str_sub(3, -2) %>%
    stringr::str_split(", ", simplify = TRUE) %>%
    as.vector() %>%
    stringi::stri_remove_empty()

  # assert that the data is of the right class
  # TODO: test cvCovEst with these sparse matrix classes
  assertthat::assert_that(
    tibble::is_tibble(dat) ||
      is.data.frame(dat) ||
      is.matrix(dat) ||
      is(dat, "dgeMatrix") ||
      is(dat, "dgCMatrix")
  )

  # if only one estimator is passed in, and it doesn't have any hyperparameters,
  # then tell users not to use cvCovEst
  if (length(estimators) == 1) {
    assertthat::assert_that(
      estimators != "linearShrinkLWEst", estimators != "sampleCovEst",
      estimators != "nlShrinkLWEst", estimators != "denseLinearShrinkEst",
      msg = paste(
        "This estimator doesn't possess any hyperparameters. Run it",
        "without using cvCovEst."
      )
    )
  }

  # assert that estimators are defined in cvCovEst package
  assertthat::assert_that(
    all(
      estimators %in% c(
        "linearShrinkEst", "linearShrinkLWEst",
        "thresholdingEst", "sampleCovEst", "bandingEst",
        "taperingEst", "nlShrinkLWEst",
        "denseLinearShrinkEst", "scadEst", "poetEst", "robustPoetEst",
        "adaptiveLassoEst"
      ) == TRUE
    ),
    msg = "Only estimators implemented in the cvCovEst package can be used."
  )

  # assert that the loss is well defined
  assertthat::assert_that(
    as.character(cv_loss) %in% c("cvFrobeniusLoss", "cvMatrixFrobeniusLoss",
                                 "cvScaledMatrixFrobeniusLoss")
  )

  # assert that estimator hyperparameters are well defined
  if ("linearShrinkEst" %in% estimators) {
    assertthat::assert_that(
      all(rlang::is_bare_numeric(estimator_params$linearShrinkEst$alpha))
      == TRUE,
      all(estimator_params$linearShrinkEst$alpha >= 0) == TRUE,
      all(estimator_params$linearShrinkEst$alpha <= 1) == TRUE
    )
  }
  if ("thresholdingEst" %in% estimators) {
    assertthat::assert_that(
      all(rlang::is_bare_numeric(estimator_params$thresholdingEst$gamma))
      == TRUE,
      all(estimator_params$thresholdingEst$gamma >= 0) == TRUE
    )
  }
  if ("bandingEst" %in% estimators) {
    assertthat::assert_that(
      all(rlang::is_integer(estimator_params$bandingEst$k)) == TRUE,
      all(estimator_params$bandingEst$k >= 0) == TRUE
    )
  }
  if ("taperingEst" %in% estimators) {
    assertthat::assert_that(
      all(rlang::is_integer(estimator_params$taperingEst$k)) == TRUE,
      all(estimator_params$taperingEst$k >= 0) == TRUE,
      all(estimator_params$taperingEst$k %% 2 == 0) == TRUE
    )
  }
  if ("nlShrinkLWEst" %in% estimators) {
    assertthat::assert_that(
      nrow(dat) >= 12
    )
  }
  if ("scadEst" %in% estimators) {
    assertthat::assert_that(
      all(rlang::is_bare_numeric(estimator_params$scadEst$lambda)) == TRUE,
      all(estimator_params$scadEst$lambda >= 0) == TRUE
    )
  }
  if ("poetEst" %in% estimators) {
    assertthat::assert_that(
      all(rlang::is_integer(estimator_params$poetEst$k)) == TRUE,
      all(estimator_params$poetEst$k >= 1) == TRUE,
      all(rlang::is_bare_numeric(estimator_params$poetEst$lambda)) == TRUE,
      all(estimator_params$poetEst$lambda >= 0) == TRUE
    )
  }
  if ("robustPoetEst" %in% estimators) {
    assertthat::assert_that(
      all(rlang::is_integer(estimator_params$robustPoetEst$k)) == TRUE,
      all(estimator_params$robustPoetEst$k >= 1) == TRUE,
      all(rlang::is_bare_numeric(estimator_params$robustPoetEst$lambda))
      == TRUE,
      all(estimator_params$robustPoetEst$lambda >= 0) == TRUE,
      all(estimator_params$robustPoetEst$var_est %in% c(
        "sample", "mad", "huber"
      ) == TRUE)
    )
  }
  if ("adaptiveLassoEst" %in% estimators) {
    assertthat::assert_that(
      all(estimator_params$adaptiveLassoEst$lambda >= 0) == TRUE,
      all(rlang::is_bare_numeric(estimator_params$adaptiveLassoEst$lambda))
      == TRUE,
      all(estimator_params$adaptiveLassoEst$n >= 0) == TRUE,
      all(rlang::is_bare_numeric(estimator_params$adaptiveLassoEst$n)) == TRUE
    )
  }

  # assert that cv scheme is supported and well defined
  assertthat::assert_that(
    cv_scheme == "mc" ||
      cv_scheme == "v_fold",
    mc_split > 0,
    mc_split < 1,
    v_folds > 1,
    v_folds < nrow(dat)
  )

  # assert that center, scaling, and parallel are flags
  assertthat::assert_that(assertthat::is.flag(center))
  assertthat::assert_that(assertthat::is.flag(scale))
  assertthat::assert_that(assertthat::is.flag(parallel))

  # when not null, assert that true_cov_matis of the appropriate type and dims
  if (!is.null(true_cov_mat)) {
    assertthat::assert_that(is.matrix(true_cov_mat) ||
      is(true_cov_mat, "dgeMatrix") ||
      is(true_cov_mat, "dgCMatrix"))
    assertthat::assert_that(identical(dim(true_cov_mat)[1], ncol(dat)))
    assertthat::assert_that(identical(dim(true_cov_mat)[2], ncol(dat)))
  } else {
    TRUE
  }
}
