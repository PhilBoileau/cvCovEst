#' Naive Bootstrap Estimator of Cross Term in Scaled Frobenius Loss
#'
#' @description \code{sumNaiveBootstrap} estimates the cross (covariance) term
#'  in the scaled Frobenius loss by computing the bootstrap means of entries in
#'  the covariance matrix estimator and sample covariance matrix.
#'
#' @param estimator_fun A \code{character} indicating the covariance matrix
#'  estimator to be applied to the resampled training data.
#' @param estimates A \code{matrix} of the estimator's covariance matrix
#'   estimate over a given fold.
#' @param sample_cov_mat A \code{matrix} of the sample covariance matrix,
#'  estimated over a given fold.
#' @param train_data A \code{data.frame} (or similar) of the training data.
#' @param num_iter A positive integer defining the number of bootstrap samples
#'   to compute.
#' @param est_args Arguments corresponding to the hyperparameters of the
#'   covariance matrix estimator, \code{estimator_fun}.
#'
#' @importFrom coop covar
#' @importFrom dplyr sample_frac
#' @importFrom matrixStats sum2
#' @importFrom tibble as_tibble
#'
#' @return A \code{numeric} estimating the cross (covariance) term of the
#'  scaled Frobenius loss.
#'
#' @keywords internal
sumNaiveBootstrap <- function(estimator_fun,
                              estimates,
                              sample_cov_mat,
                              train_data,
                              num_iter,
                              est_args = NULL) {
  # get sequence of bootstrap samples
  idx <- seq_len(num_iter)

  # convert training data and validation data to tibbles
  train_data <- suppressMessages(
    tibble::as_tibble(train_data, .name_repair = "universal")
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

  # resample validation data num_iter times, estimating the sample covariance
  # matrix on each iteration
  sample_cov_list <- lapply(idx, function(x) {
    coop::covar(sample_frac(train_data, replace = TRUE))
  })

  # compute the bootstrap means of the estimates_list and sample_cov_list
  estimates_boot_mean <- (1 / num_iter) * Reduce(`+`, estimates_list)
  sample_cov_boot_mean <- (1 / num_iter) * Reduce(`+`, sample_cov_list)

  # "center" estimates and sample_cov
  estimates <- estimates - estimates_boot_mean
  sample_cov_mat <- sample_cov_mat - sample_cov_boot_mean

  # compute hadamard product of estimates and sample cov mat
  boot_mat_mult_mat <- estimates * sample_cov_mat

  # return the sum of all the bootstraped covariance estimates
  out <- matrixStats::sum2(boot_mat_mult_mat)
  return(out)
}
