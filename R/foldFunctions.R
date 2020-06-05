#' Cross-Validation Function for Scaled Frobenius Loss
#'
#' @description \code{cvFrobeniusLoss} evaluates the scaled Frobenius loss over
#'  a \code{fold} object (from \pkg{origami} \insertRef{Coyle2018}{cvCovEst}).
#'
#' @param fold A \code{fold} object (from \code{\link[origami]{make_folds}})
#'  over which the estimation procedure is to be performed.
#' @param dat A \code{data.frame} containing the full (non-sample-split) data,
#'  on which the cross-validated procedure is performed.
#' @param estimator_funs An \code{expression} corresponding to a vector of
#'  covariance matrix estimator functions to be applied to the training data.
#' @param estimator_params A named \code{list} of arguments corresponding to
#'  the hyperparameters of covariance matrix estimators, \code{estimator_funs}.
#'  The name of each list element should be the name of an estimator passed to
#'  \code{estimator_funs}. Each element of the \code{estimator_params} is
#'  itself a named \code{list}, with names corresponding to an estimators'
#'  hyperparameter(s). These hyperparameters may be in the form of a single
#'  \code{numeric} or a \code{numeric} vector. If no hyperparameter is needed
#'  for a given estimator, then the estimator need not be listed.
#'
#' @importFrom coop covar
#' @importFrom dplyr bind_rows
#' @importFrom origami training validation fold_index
#' @importFrom Rdpack reprompt
#' @importFrom tibble tibble
#' @importFrom rlang enexpr
#' @importFrom rlang parse_expr
#'
#' @return A \code{\link[tibble]{tibble}} providing information on estimators,
#'  their hyperparameters (if any), and their scaled Frobenius loss evaluated
#'  on a given \code{fold}.
#'
#' @keywords internal
cvFrobeniusLoss <- function(fold, dat,
                            resample_fun = sumNaiveBootstrap, resample_iter,
                            estimator_funs, estimator_params = NULL) {

  # split the data into training and validation
  train_data <- origami::training(dat)
  valid_data <- origami::validation(dat)

  # compute the sample covariance matrix over the validation set
  sample_cov_mat <- coop::covar(valid_data)

  # get number of estimators
  num_estimators <- seq(from = 2, to  = length(estimator_funs))

  # loop through estimator functions
  est_out <- lapply(num_estimators, function(x) {

    # extract estimator function and name
    est_fun <- eval(estimator_funs[[x]])
    est_name <- as.character(estimator_funs[[x]])

    # check if a hyperparameter is needed
    hyp_name <- names(estimator_params[[est_name]])
    if (is.null(hyp_name)) {

      # fit the covariance matrix estimator on the training set
      est_mat <- est_fun(train_data)

      # indicate that there are no hyperparameters
      estimator_hparam <- "hyperparameters = NA"

      # return the results from the fold
      out <- tibble::tibble(
        estimator = est_name,
        hyperparameters = estimator_hparam,
        loss = scaledFrobeniusLoss(est_mat, sample_cov_mat),
        fold = origami::fold_index(fold = fold)
      )

      return(out)

    } else {

      # loop through the estimator hyperparameters
      param_out <- lapply(
        estimator_params[[est_name]][[hyp_name]],
        function(param) {

          # fit the covariance matrix estimator on the training set
          estimator_hparam <- paste(hyp_name, "=", param)
          est_mat <- est_fun(
            train_data,
            eval(rlang::parse_expr(estimator_hparam))
          )

          # return the results from the fold
          out <- list(
            estimator = est_name,
            hyperparameters = estimator_hparam,
            loss = scaledFrobeniusLoss(est_mat, sample_cov_mat),
            fold = origami::fold_index(fold = fold)
          )

          return(out)

        }
      )

      # return data frame of estimator for all considered hyperparameters
      param_out <- dplyr::bind_rows(param_out)
    }
  })

  # combine into a tibble, but wrap in list for origami before returning
  est_out <- dplyr::bind_rows(est_out)
  return(list(est_out))
}


#' Cross-Validation Function for Penalized Scaled Frobenius Loss
#'
#' @description \code{cvFrobeniusLoss} evaluates the penalized scaled Frobenius
#'  loss over a \code{fold} object (from \pkg{origami}
#'  \insertRef{Coyle2018}{cvCovEst}).
#'
#' @param fold A \code{fold} object (from \code{\link[origami]{make_folds}})
#'  over which the estimation procedure is to be performed.
#' @param dat A \code{data.frame} containing the full (non-sample-split) data,
#'  on which the cross-validated procedure is performed.
#' @param estimator_funs An \code{expression} corresponding to a vector of
#'  covariance matrix estimator functions to be applied to the training data.
#' @param resample_fun The \code{function} defining the resampling-based
#'  procedure used to estimate the entries of the cross (covariance) term in
#'  the scaled Frobenius loss.
#' @param resample_iter A \code{numeric} indicating the number of repetitions
#'  to be performed by \code{resample_fun}.
#' @param estimator_params A named \code{list} of arguments corresponding to
#'  the hyperparameters of covariance matrix estimators, \code{estimator_funs}.
#'  The name of each list element should be the name of an estimator passed to
#'  \code{estimator_funs}. Each element of the \code{estimator_params} is
#'  itself a named \code{list}, with names corresponding to an estimators'
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
#' @return A \code{\link[tibble]{tibble}} providing information on estimators,
#'  their hyperparameters (if any), and their scaled Frobenius loss evaluated
#'  on a given \code{fold}.
#'
#' @keywords internal
cvPenFrobeniusLoss <- function(fold, dat,
                               resample_fun = sumNaiveBootstrap, resample_iter,
                               estimator_funs, estimator_params = NULL) {

  # split the data into training and validation
  train_data <- origami::training(dat)
  valid_data <- origami::validation(dat)

  # compute the sample covariance matrix over the validation set
  sample_cov_mat <- coop::covar(valid_data)

  # get number of estimators
  num_estimators <- seq(from = 2, to  = length(estimator_funs))

  # loop through estimator functions
  est_out <- lapply(num_estimators, function(x) {

    # extract estimator function and name
    est_fun <- eval(estimator_funs[[x]])
    est_name <- as.character(estimator_funs[[x]])

    # check if a hyperparameter is needed
    hyp_name <- names(estimator_params[[est_name]])
    if (is.null(hyp_name)) {

      # fit the covariance matrix estimator on the training set
      est_mat <- est_fun(train_data)

      # estimate the sum of covariance terms of Cov(est_mat, sample_cov_mat)
      cov_sum <- resample_fun(
        est_fun, est_mat, sample_cov_mat,
        train_data, resample_iter
      )

      # return the results from the fold
      out <- tibble::tibble(
        estimator = est_name,
        hyperparameters = "hyperparameters = NA",
        loss = penScaledFrobeniusLoss(est_mat, sample_cov_mat, cov_sum),
        fold = origami::fold_index(fold = fold)
      )

      return(out)

    } else {
      # loop through the estimator hyperparameters
      param_out <- lapply(
        estimator_params[[est_name]][[hyp_name]],
        function(param) {

          # fit the covariance matrix estimator on the training set
          # fit the covariance matrix estimator on the training set
          estimator_hparam <- paste(hyp_name, "=", param)
          est_mat <- est_fun(
            train_data,
            eval(rlang::parse_expr(estimator_hparam))
          )

          # estimate the sum of covariance terms of Cov(est_mat, sample_cov_mat)
          cov_sum <- resample_fun(
            est_fun, est_mat, sample_cov_mat,
            train_data, resample_iter,
            eval(rlang::parse_expr(estimator_hparam))
          )

          # return the results from the fold
          out <- list(
            estimator = est_name,
            hyperparameters = estimator_hparam,
            loss = penScaledFrobeniusLoss(est_mat, sample_cov_mat, cov_sum),
            fold = origami::fold_index(fold = fold)
          )
          return(out)
        }
      )

      # return data frame of estimator for all considered hyperparameters
      param_out <- dplyr::bind_rows(param_out)
    }
  })

  # combine into a tibble, but wrap in list for origami before returning
  est_out <- dplyr::bind_rows(est_out)
  return(list(est_out))
}
