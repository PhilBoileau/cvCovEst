#' Cross-Validated Covariance Matrix Estimator Selector
#'
#' @description \code{cvCovEst()} identifies the optimal covariance matrix
#'   estimator from among a set of candidate estimators.
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#' @param estimators A \code{list} of estimator functions to be considered in
#'  the cross-validated estimator selection procedure.
#' @param estimator_params A named \code{list} of arguments corresponding to
#'  the hyperparameters of covariance matrix estimators in \code{estimators}.
#'  The name of each list element should match the name of an estimator passed
#'  to \code{estimators}. Each element of the \code{estimator_params} is itself
#'  a named \code{list}, with the names corresponding to a given estimator's
#'  hyperparameter(s). The hyperparameter(s) may be in the form of a single
#'  \code{numeric} or a \code{numeric} vector. If no hyperparameter is needed
#'  for a given estimator, then the estimator need not be listed.
#' @param cv_loss A \code{function} indicating the loss function to be used.
#'  This defaults to the Frobenius loss, \code{\link{cvMatrixFrobeniusLoss}()}.
#'  An observation-based version, \code{\link{cvFrobeniusLoss}()}, is also made
#'  available. Additionally, the \code{\link{cvScaledMatrixFrobeniusLoss}()} is
#'  included for situations in which \code{dat}'s variables are of different
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
#' @param center A \code{logical} indicating whether to center the columns of
#'  \code{dat} to have mean zero.
#' @param scale A \code{logical} indicating whether to scale the columns of
#'  \code{dat} to have unit variance.
#' @param parallel A \code{logical} option indicating whether to run the main
#'  cross-validation loop with \code{\link[future.apply]{future_lapply}()}. This
#'  is passed directly to \code{\link[origami]{cross_validate}()}.
#'
#' @importFrom origami cross_validate folds_vfold folds_montecarlo
#' @importFrom dplyr arrange summarise group_by "%>%" ungroup
#' @importFrom tibble as_tibble
#' @importFrom rlang .data enquo eval_tidy
#' @importFrom matrixStats sum2
#' @importFrom purrr flatten
#' @importFrom stringr str_split
#'
#' @return A \code{list} of results containing the following elements:
#'   \itemize{
#'     \item \code{estimate} - A \code{matrix} corresponding to the estimate of
#'       the optimal covariance matrix estimator.
#'     \item \code{estimator} - A \code{character} indicating the optimal
#'       estimator and corresponding hyperparameters, if any.
#'     \item \code{risk_df} - A \code{\link[tibble]{tibble}} providing the
#'       cross-validated risk estimates of each estimator.
#'     \item \code{cv_df} - A \code{\link[tibble]{tibble}} providing each
#'       estimators' loss over the folds of the cross-validated procedure.
#'     \item \code{args} - A named \code{list} containing arguments passed to
#'       \code{cvCovEst}.
#'   }
#'
#' @examples
#' cvCovEst(
#'   dat = mtcars,
#'   estimators = c(
#'     linearShrinkLWEst, thresholdingEst, sampleCovEst
#'   ),
#'   estimator_params = list(
#'     thresholdingEst = list(gamma = seq(0.1, 0.3, 0.1))
#'   ),
#'   center = TRUE,
#'   scale = TRUE
#' )
#' @export
cvCovEst <- function(
  dat,
  estimators = c(
   linearShrinkEst, thresholdingEst, sampleCovEst
  ),
  estimator_params = list(
   linearShrinkEst = list(alpha = 0),
   thresholdingEst = list(gamma = 0)
  ),
  cv_loss = cvMatrixFrobeniusLoss,
  cv_scheme = "v_fold", mc_split = 0.5, v_folds = 10L,
  center = TRUE,
  scale = FALSE,
  parallel = FALSE
) {

  # grab estimator expression
  estimators <- rlang::enquo(estimators)

  # grab the loss function
  cv_loss <- rlang::enquo(cv_loss)

  # check inputs
  checkArgs(
    dat,
    rlang::quo_get_expr(estimators), estimator_params,
    rlang::quo_get_expr(cv_loss), cv_scheme,
    mc_split, v_folds,
    center, scale, parallel
  )

  # create arguments list
  args <- list(
    cv_loss = cv_loss,
    cv_scheme = cv_scheme,
    mc_split = mc_split,
    v_folds = v_folds,
    center = center,
    scale = scale,
    parallel = parallel
  )

  # center and scale the data, if desired
  if (center == FALSE) {
    col_means <- colMeans(dat)
    abs_diff_zero <- abs(col_means - rep(0, length(col_means)))
    if (any(abs_diff_zero > 1e-10)) {
      message("The columns of argument `dat` have been centered automatically")
      dat <- safeColScale(X = dat, center = center, scale = scale)
    }
  } else {
    dat <- safeColScale(X = dat, center = center, scale = scale)
  }

  # define the folds based on cross-validation scheme
  if (cv_scheme == "mc") {
    folds <- origami::make_folds(dat,
      fold_fun = origami::folds_montecarlo,
      V = v_folds,
      pvalidation = mc_split
    )
  } else if (cv_scheme == "v_fold") {
    folds <- origami::make_folds(dat,
      fold_fun = origami::folds_vfold,
      V = v_folds
    )
  }

  # apply the estimators to each fold
  fold_results <- origami::cross_validate(
    dat = dat,
    cv_fun = rlang::eval_tidy(cv_loss),
    folds = folds,
    estimator_funs = estimators,
    estimator_params = estimator_params,
    use_future = parallel,
    .combine = FALSE
  )

  # convert results to tibble
  fold_results_concat <- dplyr::bind_rows(fold_results[[1]])

  # compute cv risk
  cv_results <- fold_results_concat %>%
    dplyr::group_by(.data$estimator, .data$hyperparameters) %>%
    dplyr::summarise(cv_risk = mean(.data$loss)) %>%
    dplyr::arrange(.data$cv_risk) %>%
    dplyr::ungroup()

  # compute the best estimator's estimate
  best_est_fun <- get(cv_results[1, ]$estimator)
  best_est_hparams <- cv_results[1, ]$hyperparameters
  if (best_est_hparams != "hyperparameters = NA") {
    best_est_hparams_table <- best_est_hparams %>%
      stringr::str_split(pattern = ", ") %>%
      purrr::flatten() %>%
      stringr::str_split(pattern = " = ", simplify = TRUE)
    best_hparams_list <- as.list(best_est_hparams_table[, 2])
    names(best_hparams_list) <- best_est_hparams_table[, 1]
    best_hparams_list <- lapply(best_hparams_list, strToNumber)
    estimate <- rlang::exec(
      best_est_fun,
      dat,
      !!!best_hparams_list
    )
  } else {
    estimate <- best_est_fun(dat)
  }

  # prep output and return
  out <- list(
    estimate = estimate,
    estimator = paste0(
      cv_results[1, ]$estimator, ", ",
      cv_results[1, ]$hyperparameters
    ),
    risk_df = cv_results,
    cv_df = fold_results_concat,
    args = args
  )

  class(out) <- "cvCovEst"
  return(out)
}

###############################################################################

#' Convert String to Numeric or Integer When Needed
#'
#' @param x A \code{character} representing a number or an integer.
#'
#' @return \code{x} converted to the appropriate type.
#'
#' @importFrom stringr str_sub
#'
#' @keywords internal
strToNumber <- function(x) {
  if (stringr::str_sub(x, start = -1) == "L") {
    as.integer(stringr::str_sub(x, end = -2))
  } else if (!grepl("^[[:digit:]]", x)) {
    x
  } else {
    as.numeric(x)
  }
}
