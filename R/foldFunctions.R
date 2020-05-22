#' Cross-Validation Function for Scaled Frobenius Loss
#'
#' @description \code{cvFrobeniusLoss} evaluates the scaled Frobenius loss over
#'  a \code{fold} object defined via the \pkg{origami} package
#'  \insertRef{Coyle2018}{cvCovEst}.
#'
#' @param fold A \code{fold} object over which the estimation procedure will
#'  take place.
#' @param dat The full \code{data.frame} on which the cross-validated procedure
#'  is performed.
#' @param estimator_funs A \code{list} of covariance matrix estimator functions
#'  to be applied to the training data. Functions should be input as
#'  \code{character}s.
#' @param resample_fun The \code{function} defining the resampling-based
#'  procedure used to estimate the entries of the second term in the scaled
#'  Frobenius loss.
#' @param resample_iter A \code{numeric} indicating the number of repetitions
#'  to be performed by \code{resample_fun}.
#' @param estimator_params A named \code{list} of arguments corresponding to the
#'  hyperparameters of the covariance matrix estimator, \code{estimator_funs}.
#'  The name of each list element should be the name of an estimator passed to
#'  \code{estimator_funs}. Each element of the \code{estimator_params} is
#'  itself a named \code{list}, where the names correspond to an estimators'
#'  hyperparameter(s). These hyperparameters may be in the form of a single
#'  \code{numeric} or a \code{numeric} vector. If no hyperparameter is needed
#'  for a given estimator, then the estimator need not be listed.
#'
#' @importFrom coop covar
#' @importFrom dplyr bind_rows
#' @importFrom origami training validation fold_index
#' @importFrom Rdpack reprompt
#' @importFrom tibble tibble
#'
#' @return A \code{\link[tibble]{tibble}} providing information on estimators, their
#'   hyperparameters (if any), and their scaled Frobenius loss over
#'   a \code{fold}.
#'
#' @keywords internal
cvFrobeniusLoss <- function(
  fold, dat,
  resample_fun = sumNaiveBootstrap, resample_iter,
  estimator_funs, estimator_params = NULL) {

  # split the data into training and validation
  train_data <- origami::training(dat)
  valid_data <- origami::validation(dat)

  # compute the sample covariance matrix over the validation set
  sample_cov_mat <- coop::covar(valid_data)

  # loop through estimator functions
  est_out <- lapply(
    estimator_funs,
    function(est_fun) {

      est_name <- est_fun
      est_fun <- get(est_fun)

      # check if a hyperparameter is needed
      hyp_name <- names(estimator_params[[est_name]])
      if (is.null(hyp_name)) {
        # fit the covariance matrix estimator on the training set
        est_mat <- est_fun(train_data)

        # estimate the sum of covariance terms of Cov(est_mat, sample_cov_mat)
        cov_sum <- resample_fun(est_fun, est_mat, sample_cov_mat,
                                train_data, valid_data, resample_iter)

        # indicate that there are no hyperparameters
        estimator_hparams <- "hyperparameters = NA"

        # return the results from the fold
        out <- tibble::tibble(
          estimator = est_name,
          hyperparameters = estimator_hparams,
          loss = scaledFrobeniusLoss(est_mat, sample_cov_mat, cov_sum),
          fold = origami::fold_index(fold = fold)
        )
        return(out)

      } else {

        # loop through the estimator hyperparameters
        param_out <- lapply(
          estimator_params[[est_name]][[hyp_name]],
          function(param) {

            # fit the covariance matrix estimator on the training set
            est_mat <- est_fun(
              train_data,
              eval(parse(text = paste(hyp_name, "=", param)))
            )

            # estimate the sum of covariance terms of Cov(est_mat, sample_cov_mat)
            cov_sum <- resample_fun(
              est_fun, est_mat, sample_cov_mat,
              train_data, valid_data, resample_iter,
              eval(parse(text = paste(hyp_name, "=", param)))
            )

            # get the list of estimator params
            estimator_hparams <- paste(hyp_name, "=", param)
            if (length(estimator_hparams) == 0)
              estimator_hparams <- "hyperparameters = NA"

            # return the results from the fold
            out <- list(
              estimator = est_name,
              hyperparameters = estimator_hparams,
              loss = scaledFrobeniusLoss(est_mat, sample_cov_mat, cov_sum),
              fold = origami::fold_index(fold = fold)
            )
            return(out)
          }
        )

        # return data frame of estimator for all considered hyperparameters
        param_out <- dplyr::bind_rows(param_out)
      }

    }
  )

  # combine into a tibble, but wrap in list for origami before returning
  est_out <- dplyr::bind_rows(est_out)
  return(list(est_out))
}
