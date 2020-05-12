#' Cross-Validated Covariance Matrix Estimator Selector
#'
#' @description \code{cvCovEst} identifies the optimal covariance matrix
#'   estimator from among a set of candidate estimators over a high-dimensional
#'   data set.
#'
#' @param dat A numeric \code{data.frame} or matrix.
#' @param estimators A \code{list} of estimator functions.
#' @param estimator_params A named \code{list} of arguments corresponding to the
#'   hyperparameters of the covariance matrix estimator, \code{estimator_funs}.
#'   The name of each list element should be the name of an estimator passed to
#'   \code{estimator_funs}. Each element of the \code{estimator_params} is
#'   itself a named \code{list}, where the names correspond to an estimators'
#'   hyperparameter(s). These hyperparameters may be in the form of a single
#'   \code{numeric} or a \code{numeric} vector. If no hyperparameter is needed
#'   for a given estimator, then the estimator need not be listed.
#' @param cv_scheme A \code{character} indicating the cross-validation scheme
#'   to be emplyed. There are two options: (1) V-fold cross-validation as
#'   \code{"v_folds"} and (2) Monte Carlo cross-validation as \code{"mc"}.
#'   Defaults to \code{"mc"}.
#' @param mc_split A \code{numeric} between 0 and 1 indicating the proportion of
#'   data in the validation set of each Monte Carlo cross-validation fold.
#' @param v_folds A \code{integer} larger than or equal to 1 indicating the
#'   number of folds to use during cross-validation. The default is 10,
#'   regardless of cross-validation scheme.
#' @param boot_iter A \code{integer} dictating the number of bootstrap
#'   iterations used to compute the covariance terms of the cross-validated
#'   scaled Frobenius loss. The default is set to 100.
#' @param center A \code{logical} indicating whether or not to center the
#'   columns of \code{dat}.
#' @param scale A \code{logical} indicating whether or not to scale the
#'   columns of \code{dat} to have variance 1.
#' @param parallel A \code{logical} option for whether to run the main
#'   cross-validation loop with \code{\link[future.apply]{future_lapply}}.
#'
#' @return A \code{list} of results containing the following elements:
#'   \itemize{
#'     \item \code{estimate} - A \code{matrix} corresponding to the estimate of
#'       the optimal covariance matrix estimator. (TODO)
#'     \item \code{estimator} - A \code{character} indicating the optimal
#'       estimator and corresponding hyperparameters, if any. (TODO)
#'     \item \code{results_df} - A \code{tibble} providing the results of
#'       the cross-validation procedure. (TODO)
#'     \item \code{origami_output} - A \code{tibble} providing the results of
#'       the \code{\link[origami]{cross_validate}} call.
#'   }
#' @export
#'
#' @importFrom origami folds_montecarlo
#' @importFrom origami folds_vfold
#' @importFrom origami cross_validate
#' @importFrom tidyr as_tibble
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom magrittr %>%
#' @importFrom rlang .data
cvCovEst <- function(
  dat,
  estimators = c("linearShrinkEst", "thresholdingEst", "sampleCovEst"),
  estimator_params = list("linearShrinkEst" = list("alpha" = 0),
                          "thresholdingEst" = list("gamma" = 0)),
  cv_scheme = "mc", mc_split = 0.5, v_folds = 10,
  boot_iter = 100,
  center = TRUE, scale = TRUE,
  parallel = FALSE) {

  # center and scale the data, if desired. (TODO: efficient implementation?)
  dat <- scale(dat, center = center, scale = scale)

  # define the folds based on cross-validation scheme
  n_obs <- nrow(dat)
  if (cv_scheme == "mc") {
    folds <- origami::make_folds(dat,
                                 fold_fun = folds_montecarlo,
                                 V = v_folds,
                                 pvalidation = mc_split)
  } else if (cv_scheme == "v_fold") {
    folds <- origami::make_folds(dat,
                                 fold_fun = folds_vfold,
                                 V = v_folds)
  }

  # apply the estimators to each fold
  fold_results <- origami::cross_validate(
    dat = dat,
    cv_fun = cvFrobeniusLoss,
    folds = folds,
    estimator_funs = estimators,
    estimator_params = estimator_params,
    resample_iter = boot_iter,
    use_future = parallel
  )

  # remove error list from cv_results
  errors <- fold_results$errors
  fold_results$errors <- NULL

  # fix tpes
  fold_results$loss <- as.numeric(fold_results$loss)
  fold_results$fold <- as.numeric(fold_results$fold)

  # turn results to tibble
  fold_results <- tidyr::as_tibble(fold_results)

  # compute empirical risk
  cv_results <- fold_results %>%
    dplyr::group_by(.data$estimator, .data$hyperparameters) %>%
    dplyr::summarise(empirical_risk = mean(.data$loss))

  # prep output
  out <- list(
    risk_df = cv_results,
    cv_df = fold_results
  )

  # retun dataframe of cv results
  return(out)
}
