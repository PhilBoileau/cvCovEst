#' Cross-Validated Covariance Matrix Estimator Selector
#'
#' @description \code{cvCovEst} identifies the optimal covariance matrix
#'   estimator from among a set of candidate estimators.
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
#' @param cv_scheme A \code{character} indicating the cross-validation scheme
#'  to be employed. There are two options: (1) V-fold cross-validation, via
#'  \code{"v_folds"}; and (2) Monte Carlo cross-validation, via \code{"mc"}.
#'  Defaults to Monte Carlo cross-validation.
#' @param mc_split A \code{numeric} between 0 and 1 indicating the proportion
#'  of data in the validation set of each Monte Carlo cross-validation fold.
#' @param v_folds A \code{integer} larger than or equal to 1 indicating the
#'  number of folds to use during cross-validation. The default is 10,
#'  regardless of cross-validation scheme.
#' @param cv_loss A \code{function} indicating the loss function to use.
#'  Defaults to the penalized scaled Frobenius loss, \code{cvPenFrobeniusLoss}.
#'  The non-penalized version, \code{cvFrobeniusLoss} is offered as well.
#' @param boot_iter A \code{integer} dictating the number of bootstrap
#'  iterations used to compute the penalty term of the cross-validated
#'  penalized scaled Frobenius loss. The default is set to 100. If
#'  \code{cvFrobeniusLoss} is selected in place of \code{cvPenFrobeniusLoss},
#'  then this argument is ignored.
#' @param center A \code{logical} indicating whether or not to center the
#'  columns of \code{dat}.
#' @param scale A \code{logical} indicating whether or not to scale the
#'  columns of \code{dat} to have variance 1.
#' @param parallel A \code{logical} option indicating whether to run the main
#'  cross-validation loop with \code{\link[future.apply]{future_lapply}}. This
#'  is passed directly to \code{\link[origami]{cross_validate}}.
#'
#' @importFrom origami cross_validate folds_vfold folds_montecarlo
#' @importFrom dplyr arrange summarise group_by "%>%"
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#' @importFrom rlang as_string
#' @importFrom rlang expr
#'
#' @return A \code{list} of results containing the following elements:
#'   \itemize{
#'     \item \code{estimate} - A \code{matrix} corresponding to the estimate of
#'       the optimal covariance matrix estimator.
#'     \item \code{estimator} - A \code{character} indicating the optimal
#'       estimator and corresponding hyperparameters, if any.
#'     \item \code{results_df} - A \code{\link[tibble]{tibble}} providing the
#'       results of the cross-validation procedure. (TODO)
#'     \item \code{origami_output} - A \code{\link[tibble]{tibble}} providing
#'       the results of the \code{\link[origami]{cross_validate}} call.
#'   }
#'
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
  cv_scheme = "mc", mc_split = 0.5, v_folds = 10L,
  cv_loss = cvFrobeniusLoss,
  boot_iter = 100L,
  center = TRUE,
  scale = TRUE,
  parallel = FALSE
)
{

  # grab estimator expression
  estimators <- rlang::enexpr(estimators)

  #grab the cv_loss function name as astring
  cv_loss_name <- rlang::enexpr(cv_loss)

  # check inputs
  checkArgs(
    dat,
    estimators, estimator_params,
    cv_scheme, mc_split, v_folds, cv_loss_name, boot_iter,
    center, scale, parallel
  )

  # center and scale the data, if desired
  dat <- safeColScale(X = dat, center = center, scale = scale)

  # define the folds based on cross-validation scheme
  n_obs <- nrow(dat)
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
  if (rlang::as_string(rlang::enexpr(cv_loss)) == "cvPenFrobeniusLoss") {

    fold_results <- origami::cross_validate(
      dat = dat,
      cv_fun = cv_loss,
      folds = folds,
      estimator_funs = estimators,
      estimator_params = estimator_params,
      resample_iter = boot_iter,
      use_future = parallel,
      .combine = FALSE
    )

  } else if (rlang::as_string(rlang::enexpr(cv_loss)) == "cvFrobeniusLoss"){

    fold_results <- origami::cross_validate(
      dat = dat,
      cv_fun = cv_loss,
      folds = folds,
      estimator_funs = estimators,
      estimator_params = estimator_params,
      use_future = parallel,
      .combine = FALSE
    )

  }

  # convert results to tibble
  fold_results_concat <- dplyr::bind_rows(fold_results[[1]])
  fold_results_concat

  # compute empirical risk
  cv_results <- fold_results_concat %>%
    dplyr::group_by(.data$estimator, .data$hyperparameters) %>%
    dplyr::summarise(empirical_risk = mean(.data$loss)) %>%
    dplyr::arrange(.data$empirical_risk)

  # compute the best estimator's estimate
  best_est <- get(cv_results[1, ]$estimator)
  best_est_hyperparams <- parse(text = cv_results[1, ]$hyperparameters)
  if (cv_results[1, ]$hyperparameters != "hyperparameters = NA") {
    estimate <- best_est(dat, eval(best_est_hyperparams))
  } else {
    estimate <- best_est(dat)
  }

  # prep output and return
  out <- list(
    estimate = estimate,
    estimator = paste0(
      cv_results[1, ]$estimator, ", ",
      cv_results[1, ]$hyperparameters
    ),
    risk_df = cv_results,
    cv_df = fold_results_concat
  )
  return(out)
}
