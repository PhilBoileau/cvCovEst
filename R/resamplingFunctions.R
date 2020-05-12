#' Naive Bootstrap Estimator of Second Term in Scaled Frobenius Loss
#'
#' @description \code{sumNaiveBootstrap} estimates the second term in scaled
#'   Frobenius loss by computing the bootstrap means of the entries in the
#'   covariance matrix  estimator and sample covariance matrix.
#'
#' @param estimator_fun The covariance matrix estimator to be applied to the
#'   resampled training data.
#' @param estimates A \code{matrix} of the estimator's covariance matrix
#'   estimate over a given fold.
#' @param sample_cov_mat A \code{matrix} of the sample covariance matrix over a
#'   given fold.
#' @param train_data A \code{data.frame} of the training data.
#' @param valid_data A \code{data.frame} of the validation data.
#' @param num_iter A positive integer defining the number of bootstrap samples
#'   to compute.
#' @param est_args Arguments corresponding to the hyperparameters of the
#'   covariance matrix estimator, \code{estimator_fun}.
#'
#' @return A \code{numeric} estimating the second term of the scaled Frobenius
#'   loss.
#'
#' @importFrom coop covar
#' @importFrom dplyr sample_frac
#' @importFrom matrixStats sum2
#' @importFrom tidyr as_tibble
#'
#' @keywords internal
sumNaiveBootstrap <- function(estimator_fun, estimates, sample_cov_mat,
                              train_data, valid_data,
                              num_iter, est_args = NULL) {

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
  if (is.null(est_args)) {
    estimates_list <- lapply(
      idx,
      function(x) {
        estimator_fun(sample_frac(train_data, replace = TRUE))
      }
    )
  } else {
    estimates_list <- lapply(
      idx,
      function(x) {
        estimator_fun(sample_frac(train_data, replace = TRUE), est_args)
      }
    )
  }

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

  # "center" estimates and sample_cov
  estimates <- estimates - estimates_boot_mean
  sample_cov_mat <- sample_cov_mat - sample_cov_boot_mean

  # compute hadamard produce of estimats and sample cov mat
  boot_mat_mult_mat <- estimates * sample_cov_mat

  # return the sum of all the bootstraped covariance estimates
  return(matrixStats::sum2(boot_mat_mult_mat))
}
