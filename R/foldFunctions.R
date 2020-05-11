#' Cross-Validation Function for Scaled Frobenius Loss
#'
#' @description \code{cvFrobeniusLoss} evaluates the scaled Frobenius loss over
#'   a \code{fold} object defined via the \pkg{origami} package
#'   \insertRef{Coyle2018}{cvCovEst}.
#'
#' @param fold A \code{fold} object over which the estimation procedure will
#'   take place.
#' @param dat The full \code{data.frame} on which the cross-validated procedure
#'   is performed.
#' @param estimator_funs A \code{list} of covariance matrix estimator functions
#'   to be applied to the training data.
#' @param resample_fun The function defining the resampling-based procedure used
#'   to estimate the covariances of entries in the covariance matrix estimator
#'   and the sample covariance matrix.
#' @param resample_iter A \code{numeric} indicating the number of repetitions
#'   to be performed by \code{resample_fun}.
#' @param ... Arguments corresponding to the hyperparameters of the
#'   covariance matrix estimator, \code{estimator_fun}.
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
cvFrobeniusLoss <- function(fold, dat,
                            resample_cov_fun, resample_iter,
                            estimator_funs, ...) {

  # split the data into training and validation
  train_data <- origami::training(dat)
  valid_data <- origami::validation(dat)

  # compute the sample covariance matrix over the validation set
  sample_cov_mat <- coop::covar(valid_data)

  est_out <- lapply(
    estimator_funs,
    function(est_fun) {

      # fit the covariance matrix estimator on the training set
      est_name <- est_fun
      est_fun <- get(est_fun)
      est_mat <- est_fun(train_data, ...)

      # estimate the sum of covariance terms of Cov(est_mat, sample_cov_mat)
      cov_sum <- resample_cov_fun(est_fun, train_data, valid_data,
                                  resample_iter, ...)

      # get the list of estimator params
      estimator_hparams <- list(...)
      if (length(estimator_hparams) == 0)
        estimator_hparams <- list("hyperparameters" = "NA")

      # return the results # return the results from the fold
      out <- list(
        estimator = est_name,
        hyperparameters = paste(names(estimator_hparams), estimator_hparams,
                                sep = "=", collapse = ", "),
        loss = scaledFrobeniusLoss(est_mat, sample_cov_mat, cov_sum),
        fold = fold_index(fold = fold)
      )
    })

  # combine lists into datafra,e
  keys <- unique(unlist(lapply(est_out, names)))
  results <- as.data.frame(as.list(
    setNames(do.call(mapply, c(FUN = c, lapply(est_out, `[`, keys))), keys)
  ))

  return(results)
}
