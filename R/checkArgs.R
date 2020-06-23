#' Check Arguments Passed to cvCovEst
#'
#' @description The \code{checkArgs} function verifies that all arguments passed
#'  to the \code{cvCovEst} function meet its specifications.
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
#' @return Whether all argument conditions are satisfied
#'
#' @keywords internal
checkArgs <- function(dat,
                      estimators, estimator_params,
                      cv_scheme, mc_split, v_folds,
                      cv_loss, boot_iter,
                      center, scale,
                      parallel) {

}
