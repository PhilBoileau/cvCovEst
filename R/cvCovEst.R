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
#'  Defaults to V-fold cross-validation.
#' @param mc_split A \code{numeric} between 0 and 1 indicating the proportion
#'  of data in the validation set of each Monte Carlo cross-validation fold.
#' @param v_folds A \code{integer} larger than or equal to 1 indicating the
#'  number of folds to use during cross-validation. The default is 10,
#'  regardless of cross-validation scheme.
#' @param center A \code{logical} indicating whether or not to center the
#'  columns of \code{dat}. Set to \code{FALSE} only if the columns have already
#'  been centered. Defaults to \code{TRUE}.
#' @param scale A \code{logical} indicating whether or not to scale the
#'  columns of \code{dat} to have variance 1. Defaults to \code{FALSE}.
#' @param parallel A \code{logical} option indicating whether to run the main
#'  cross-validation loop with \code{\link[future.apply]{future_lapply}}. This
#'  is passed directly to \code{\link[origami]{cross_validate}}.
#' @param true_cov_mat A \code{matrix} like object representing the true
#'  covariance matrix of the data generating distribution, which is assumed to
#'  be Gaussian. This parameter is intended for use only in simulation studies,
#'  and defaults to a value of \code{NULL}. If not \code{NULL}, the
#'  cross-validated conditional risk difference ratio of the estimator selected
#'  by \code{cvCovEst} is computed relative to the cross-validated oracle
#'  selector.
#'
#' @importFrom origami cross_validate folds_vfold folds_montecarlo
#' @importFrom dplyr arrange summarise group_by "%>%" pull
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#' @importFrom rlang enexpr
#' @importFrom matrixStats colMeans2 sum2
#' @importFrom Matrix tcrossprod
#'
#' @return A \code{list} of results containing the following elements:
#'   \itemize{
#'     \item \code{estimate} - A \code{matrix} corresponding to the estimate of
#'       the optimal covariance matrix estimator.
#'     \item \code{estimator} - A \code{character} indicating the optimal
#'       estimator and corresponding hyperparameters, if any.
#'     \item \code{risk_df} - A \code{\link[tibble]{tibble}} providing the
#'       cross-validated risk estimates of each estimator.
#'     \item \code{cv_df} - A \code{\link[tibble]{tibble}} providing
#'       each estimators' loss over the various folds of the cross-validatied
#'       procedure.
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
  cv_scheme = "v_fold", mc_split = 0.5, v_folds = 10L,
  center = TRUE,
  scale = FALSE,
  parallel = FALSE,
  true_cov_mat = NULL
)
{

  # grab estimator expression
  estimators <- rlang::enexpr(estimators)

  # check inputs
  checkArgs(
    dat,
    estimators, estimator_params,
    cv_scheme, mc_split, v_folds,
    center, scale, parallel
  )

  # center and scale the data, if desired
  if (center == FALSE) {
    col_means <- matrixStats::colMeans2(dat)
    abs_diff_zero <- abs(col_means - rep(0, length(col_means)))
    if (any(abs_diff_zero > 1e-10)) {
      message("`dat` argument's columns have been centered automatically")
      dat <- safeColScale(X = dat, center = center, scale = scale)
    }
  } else {
    dat <- safeColScale(X = dat, center = center, scale = scale)
  }

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
  fold_results <- origami::cross_validate(
      dat = dat,
      cv_fun = cvFrobeniusLoss, # might provide other options at a later date
      folds = folds,
      estimator_funs = estimators,
      estimator_params = estimator_params,
      use_future = parallel,
      .combine = FALSE
  )

  # convert results to tibble
  fold_results_concat <- dplyr::bind_rows(fold_results[[1]])
  fold_results_concat

  # compute the true cross-validated risk if true_covar_mat is passed in
  if(is.null(true_cov_mat)) {

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
  } else {

    # compute the minimum Frobenius risk, assuming Gaussian data
    d_true_cov <- diag(true_cov_mat)
    combo_mat <- Matrix::tcrossprod(d_true_cov, d_true_cov) + true_cov_mat^2
    min_full_risk <- matrixStats::sum2(combo_mat)

    # compute the true cross-validated risk and the cv-estimated risk
    cv_results <- fold_results_concat %>%
      dplyr::group_by(.data$estimator, .data$hyperparameters) %>%
      dplyr::summarise(true_cv_risk = mean(.data$true_loss),
                       empirical_risk = mean(.data$loss))  %>%
      dplyr::arrange(.data$empirical_risk)

    # compute the risk distance ratio under the cross-validated risk
    # of the cross-validated oracle and the cross-validated selection
    cvCovEst_true_cv_risk <- cv_results$true_cv_risk[1]
    oracle_true_cv_risk <- cv_results %>%
      dplyr::arrange(.data$true_cv_risk) %>%
      pull(.data$true_cv_risk)[1]
    cv_oracle_riskdiff_ratio <- (cvCovEst_true_cv_risk - min_full_risk) /
      (oracle_true_cv_risk - min_full_risk)

    # compute the cvCovEst estimator's estimate
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
      cv_df = fold_results_concat,
      cv_oracle_riskdiff_ratio = cv_oracle_riskdiff_ratio
    )
  }

  return(out)
}
