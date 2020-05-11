#' Naive Bootstrap Estimator of Covariance Matrix and Sample Covariance Matrix Covariance
#'
#' @description \code{sumNaiveCovBootstrap} estimates the covariance matrix
#'   between a generic covariance matrix estimator and the sample covariance
#'   matrix fit over a training set and a validation set, respectively. The
#'   sum of the estimated covariance matrix's elements are returned.
#'
#' @param estimator_fun The covariance matrix estimator to be applied to the
#'   resampled training data.
#' @param train_data A \code{data.frame} of the training data.
#' @param valid_data A \code{data.frame} of the validation data.
#' @param num_iter A positive integer defining the number of bootstrap samples
#'   to compute.
#' @param ... Arguments corresponding to the hyperparameters of the
#'   covariance matrix estimator, \code{estimator_fun}.
#'
#' @return A \code{numeric} representing the sum of all elements in the
#'   estimated covariance matrix of the \code{estimator_fun} and the sample
#'   covariance matrix.
#'
#' @importFrom coop covar
#' @importFrom dplyr sample_frac
#' @importFrom matrixStats sum2
#' @importFrom tidyr as_tibble
#'
#' @keywords internal
sumNaiveCovBootstrap <- function(estimator_fun, train_data, valid_data,
                                 num_iter, ...) {

  # get sequence of bootstrap samples
  idx <- seq_len(num_iter)

  # convert training data and validation data to tibbles
  train_data <- suppressMessages(
    tidyr::as_tibble(train_data, .name_repair = "universal")
    )
  valid_data <- suppressMessages(
    tidyr::as_tibble(valid_data, .name_repair = "universal")
  )

  # resample training data num_iter times, estimate covariance matrix on each
  estimates_list <- lapply(
    idx,
    function(x) {
      estimator_fun(sample_frac(train_data, replace = TRUE), ...)
    }
  )

  # resample validation data num_iter times, estimate the sample cov mat on each
  sample_cov_list <- lapply(
    idx,
    function(x) {
      coop::covar(sample_frac(valid_data, replace = TRUE))
    }
  )

  # compute the bootstrap means of the estimates_list and sample_cov_list
  estimates_boot_mean <- 1/num_iter * Reduce(`+`, estimates_list)
  sample_cov_boot_mean <- 1/num_iter * Reduce(`+`, sample_cov_list)

  # "center" estimates_list and sample_cov_list
  estimates_list <- lapply(estimates_list, function(mat) mat - estimates_boot_mean)
  sample_cov_list <- lapply(sample_cov_list, function(mat) mat - sample_cov_boot_mean)

  # multiply centered matrices sharing a bootstrap index
  boot_mat_mult_list <- lapply(
    idx,
    function(x) estimates_list[[x]] * sample_cov_list[[x]]
  )

  # return the sum of all the bootstraped covariance estimates
  return(1/num_iter * matrixStats::sum2(Reduce(`+`, boot_mat_mult_list)))
}
