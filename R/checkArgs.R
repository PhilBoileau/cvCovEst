#' Check Arguments Passed to cvCovEst
#'
#' @description \code{checkArgs()} verifies that all arguments
#'  passed to \code{\link{cvCovEst}()} function meet its specifications.
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#' @param estimators A \code{list} of estimator functions to be considered in
#'  the cross-validated estimator selection procedure.
#' @param estimator_params A named \code{list} of arguments corresponding to
#'  the hyperparameters of covariance matrix estimators in \code{estimators}.
#'  The name of each list element should match the name of an estimator passed
#'  to \code{estimators}. Each element of the \code{estimator_params} is itself
#'  a named \code{list}, with the names corresponding to a given estimator's
#'  hyperparameter(s). These hyperparameters may be in the form of a single
#'  \code{numeric} or a \code{numeric} vector. If no hyperparameter is needed
#'  for a given estimator, then the estimator need not be listed.
#' @param cv_loss A \code{function} indicating the loss function to be used.
#'  This defaults to the Frobenius loss, \code{\link{cvMatrixFrobeniusLoss}()}.
#'  An observation-based version, \code{\link{cvFrobeniusLoss}()}, is also made
#'  available. Additionally, the \code{\link{cvScaledMatrixFrobeniusLoss}(())}
#'  is included for situations in which \code{dat}'s variables are of different
#'  scales.
#' @param cv_scheme A \code{character} indicating the cross-validation scheme
#'  to be employed. There are two options: (1) V-fold cross-validation, via
#'  \code{"v_folds"}; and (2) Monte Carlo cross-validation, via \code{"mc"}.
#'  Defaults to Monte Carlo cross-validation.
#' @param mc_split A \code{numeric} between 0 and 1 indicating the proportion
#'  of observations to be included in the validation set of each Monte Carlo
#'  cross-validation fold.
#' @param v_folds An \code{integer} larger than or equal to 1 indicating the
#'  number of folds to use for cross-validation. The default is 10, regardless
#'  of the choice of cross-validation scheme.
#' @param parallel A \code{logical} option indicating whether to run the main
#'  cross-validation loop with \code{\link[future.apply]{future_lapply}()}. This
#'  is passed directly to \code{\link[origami]{cross_validate}()}.
#'
#' @importFrom methods is
#' @importFrom assertthat assert_that is.flag
#' @importFrom dplyr case_when "%>%"
#' @importFrom rlang is_bare_numeric is_integer is_null expr_deparse
#' @importFrom stringr str_sub str_split
#' @importFrom purrr keep
#'
#' @return Whether all argument conditions are satisfied
#'
#' @keywords internal
checkArgs <- function(
  dat,
  estimators, estimator_params,
  cv_loss, cv_scheme, mc_split, v_folds,
  parallel
) {

  # turn list of estimator functions to a vector of strings
  estimators <- estimators %>%
    rlang::expr_deparse() %>%
    stringr::str_sub(3, -2) %>%
    stringr::str_split(", ", simplify = TRUE) %>%
    as.list() %>%
    purrr::keep(function(x) nchar(x) > 0) %>%
    unlist()

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
        "adaptiveLassoEst", "spikedOperatorShrinkEst",
        "spikedFrobeniusShrinkEst", "spikedSteinShrinkEst"
      ) == TRUE
    ),
    msg = "Only estimators implemented in the cvCovEst package can be used."
  )

  # assert that the loss is well defined
  assertthat::assert_that(
    as.character(cv_loss) %in% c(
      "cvFrobeniusLoss", "cvMatrixFrobeniusLoss",
      "cvScaledMatrixFrobeniusLoss"
    )
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
  if ("spikedOperatorShrinkEst" %in% estimators) {
    assertthat::assert_that(
      all(estimator_params$spikedOperatorShrinkEst$p_n_ratio > 0) == TRUE,
      all(estimator_params$spikedOperatorShrinkEst$p_n_ratio < 1) == TRUE,
      all(is.null(estimator_params$spikedOperatorShrinkEst$num_spikes) |
          rlang::is_integer(estimator_params$spikedOperatorShrinkEst$num_spikes)
      ) == TRUE,
      all(is.null(estimator_params$spikedOperatorShrinkEst$num_spikes) |
            estimator_params$spikedOperatorShrinkEst$num_spikes >= 0
      ) == TRUE,
      all(is.null(estimator_params$spikedOperatorShrinkEst$noise) |
          rlang::is_bare_numeric(estimator_params$spikedOperatorShrinkEst$noise)
      ) == TRUE,
      all(is.null(estimator_params$spikedOperatorShrinkEst$noise) |
          estimator_params$spikedOperatorShrinkEst$noise > 0
      ) == TRUE
    )
  }
  if ("spikedFrobeniusShrinkEst" %in% estimators) {
    assertthat::assert_that(
      all(estimator_params$spikedFrobeniusShrinkEst$p_n_ratio > 0) == TRUE,
      all(estimator_params$spikedFrobeniusShrinkEst$p_n_ratio < 1) == TRUE,
      all(
        is.null(estimator_params$spikedFrobeniusShrinkEst$num_spikes) |
        rlang::is_integer(estimator_params$spikedFrobeniusShrinkEst$num_spikes)
      ) == TRUE,
      all(is.null(estimator_params$spikedFrobeniusShrinkEst$num_spikes) |
          estimator_params$spikedFrobeniusShrinkEst$num_spikes >= 0
      ) == TRUE,
      all(
        is.null(estimator_params$spikedFrobeniusShrinkEst$noise) |
        rlang::is_bare_numeric(estimator_params$spikedFrobeniusShrinkEst$noise)
      ) == TRUE,
      all(is.null(estimator_params$spikedFrobeniusShrinkEst$noise) |
          estimator_params$spikedFrobeniusShrinkEst$noise > 0
      ) == TRUE
    )
  }
  if ("spikedSteinShrinkEst" %in% estimators) {
    assertthat::assert_that(
      all(estimator_params$spikedSteinShrinkEst$p_n_ratio > 0) == TRUE,
      all(estimator_params$spikedSteinShrinkEst$p_n_ratio < 1) == TRUE,
      all(is.null(estimator_params$spikedSteinShrinkEst$num_spikes) |
          rlang::is_integer(estimator_params$spikedSteinShrinkEst$num_spikes)
      ) == TRUE,
      all(is.null(estimator_params$spikedSteinShrinkEst$num_spikes) |
          estimator_params$spikedSteinShrinkEst$num_spikes >= 0
      ) == TRUE,
      all(is.null(estimator_params$spikedSteinShrinkEst$noise) |
          rlang::is_bare_numeric(estimator_params$spikedSteinShrinkEst$noise)
      ) == TRUE,
      all(is.null(estimator_params$spikedSteinShrinkEst$noise) |
          estimator_params$spikedSteinShrinkEst$noise > 0
      ) == TRUE
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

  # assert that parallel is a flag
  assertthat::assert_that(assertthat::is.flag(parallel))

}
