#' Cross-Validation Function for Scaled Frobenius Loss
#'
#' @description \code{cvFrobeniusLoss} evaluates the scaled Frobenius loss over
#'   a \code{fold} object defined via the \pkg{origami} package
#'   \insertRef{Coyle2018}{cvCovEst}.
#'
#' @param fold A \code{fold} object over which the estimation procedure will
#'   take place.
#' @param data The full \code{data.frame} on which the cross-validated procedure
#'   is performed.
#' @param estimator_fun The covariance matrix estimator to be applied to the
#'   training data.
#' @param estimator_ctrl A \code{list} of hyperparameters to be passed to the
#'   covariance matrix estimator, \code{estimator_fun}.
#' @param resample_fun The function defining the resampling-based procedure used
#'   to estimate the covariances of entries in the covariance matrix estimator
#'   and the sample covariance matrix.
#' @param resample_ctrl A \code{list} of hyperparameters to be passed to
#'   \code{resample_fun}.
#'
#' @importFrom coop covar
#' @importFrom origami training
#' @importFrom origami validation
#' @importFrom Rdpack reprompt
#'
#' @return A \code{list} providing information on the estimator, its
#'   hyperparameters (if any), and scaled Frobenius loss for the
#'   given \code{fold}.
#'
#' @keywords internal
cvFrobeniusLoss <- function(fold, data,
                            estimator_fun, estimator_ctrl,
                            resample_cov_fun, resample_ctrl) {

  # split the data into training and validation
  train_data <- origami::training(data)
  valid_data <- origami::validation(data)

  # fit the covariance matrix estimator on the training set
  est_mat <- estimator_fun(train_data, estimator_ctrl)

  # compute the sample covariance matrix over the validation set
  sample_cov_mat <- coop::covar(valid_data)

  # estimate the sum of covariance terms of Cov(est_mat, sample_cov_mat)
  cov_sum <- resample_cov_fun(train_data, valid_data, resample_ctrl)

  # return the results over given fold
  out <- list(
    estimator = as.character(substitute(estimator_fun)),
    estimator_params = paste(names(estimator_ctrl), estimator_ctrl,
                             sep = "=", collapse = ", "),
    loss = scaledFrobeniusLoss(est_mat, sample_cov_mat, cov_sum)
  )
  return(out)
}
