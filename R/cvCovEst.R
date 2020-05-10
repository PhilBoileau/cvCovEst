#' Cross-Validated Covariance Matrix Estimator Selector
#'
#' @description \code{cvCovEst}
#'
#' @param dat A numeric \code{data.frame} or matrix.
#' @param estimators A \code{list} of estimator functions.
#' @param estimator_params A named \code{list} of estimator parameters.
#' @param cv_scheme A \code{character} indicating the cross-validation scheme
#'   to be emplyed. Defaults to Monte Carlo (\code{mc}).
#' @param mc_split A \code{numeric} between 0 and 1 indicating the proportion of
#'   data in the validation set for each Monte Carlo iteration.
#' @param v_folds A \code{integer} larger than or equal to 1 indicating the
#'   number of folds to use during cross-validation. The default is 10,
#'   regardless of cross-validation scheme.
#' @param boot_iter A \code{integer} dictating the number of bootstrap
#'   iterations used to compute the covariance terms of the cross-validated
#'   scaled Frobenius loss. The default is set to 100.
#' @param center A \code{logical} indicating whether or not to center the
#'   columns \code{dat}.
#' @param scale A \code{logical} indicating whether or not to scale the
#'   columns \code{dat} to have variance 1.
#' @param parallel A \code{logical} option for whether to run the main
#'   cross-validation loop with \code{\link[future.apply]{future_lapply}}.
#'
#' @return A \code{list} of results containing the following elements:
#'   \itemize{
#'     \item \code{estimate} - A \code{matrix} corresponding to the estimate of
#'       the optimal covariance matrix estimator.
#'     \item \code{estimator} - A \code{character} indicating the optimal
#'       estimator and corresponding hyperparameters, if any.
#'     \item \code{results_df} - A \code{data.frame} providing the results of
#'       the cross-validation procedure.
#'   }
#' @export
#'
#' @importFrom origami folds_montecarlo
#' @importFrom origami folds_vfold
#' @importFrom origami cross_validate
cvCovEst <- function(
  dat,
  estimators = list(linearShrinkEst, thresholdingEst),
  estimator_params = list("linearShrinkEst" = c("alpha" = 0),
                          "thresholdingEst" = c("gamma" = 0)),
  cv_scheme = "mc", mc_split = 0.5, v_folds = 10,
  boot_iter = 100,
  center = TRUE, scale = TRUE,
  parallel = FALSE) {

  # center and scale the data, if desired. (TODO: efficient implementation?)
  dat <- scale(dat, center = center, scale = scale)

  # define the folds based on cross-validation scheme
  n_obs <- nrow(dat)
  if (cv_scheme == "mc") {
    folds <- origami::folds_montecarlo(n = n_obs,
                                       V = v_folds,
                                       pvalidation = mc_split)
  } else if (cv_scheme == "v_fold") {
    folds <- origami::folds_vfold(n = n_obs,
                                  V = v_folds)
  }

  # apply the estimators to each fold
  cv_results <- origami::cross_validate(
    cv_fun = cvFrobeniusLoss,
    folds = folds,
    estimator_list = estimators,
    estimator_params = estimator_params,
    resample_iter = boot_iter,
    use_future = parallel
  )

  # return cv_results (TODO: run best estimate on data and return as well)
  return(cv_results)
}
