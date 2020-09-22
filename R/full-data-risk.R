#' Full Dataset Risk
#'
#' @description \code{fullDataRisk} computes the full dataset risk under the
#'   true data generating distribution, assuming that it is Gaussian.
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
#' @param true_cov_mat A \code{matrix} like object representing the true
#'  covariance matrix of the data generating distribution, which is assumed to
#'  be Gaussian.
#'
#' @return A \code{\link[tibble]{tibble}} providing the full dataset risk of
#'   each estimator.
#'
#' @keywords internal
fullDataRisk <- function(dat, estimators, estimator_params, true_cov_mat) {

  # fit the estimators on the full dataset

  # compute each estimator's true, full dataset Frobenius risk

  # retun the tibble

}
