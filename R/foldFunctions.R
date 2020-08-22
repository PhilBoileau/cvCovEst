#' Cross-Validation Function for Aggregated Frobenius Loss
#'
#' @description \code{cvFrobeniusLoss} evaluates the aggregated Frobenius loss
#'   over a \code{fold} object (from \pkg{origami}
#'   \insertCite{Coyle2018}{cvCovEst}).
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
#' @export
cvFrobeniusLoss <- function(fold, dat, estimator_funs,
                            estimator_params = NULL) {

  # split the data into training and validation
  train_data <- origami::training(dat)
  valid_data <- origami::validation(dat)

  # compute the validation observation crossprod matrices
  rank_one_crossp <- lapply(
    seq_len(nrow(valid_data)),
    function(i) {
      Matrix::tcrossprod(valid_data[i, ])
    }
  )

  # get number of estimators
  num_estimators <- seq(from = 2, to = length(estimator_funs))

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

      # compute the loss for each observation
      valid_losses <- sapply(
        seq_len(nrow(valid_data)),
        function(x) {

          # compute the Frobenius loss
          matrixStats::sum2((est_mat - rank_one_crossp[[x]])^2)

        }
      )

      # return the results from the fold
      out <- tibble::tibble(
        estimator = est_name,
        hyperparameters = estimator_hparam,
        loss = 1/nrow(valid_data) * sum(valid_losses),
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

          # compute the loss for each observation
          valid_losses <- sapply(
            seq_len(nrow(valid_data)),
            function(x) {

              # compute the Frobenius loss
              matrixStats::sum2((est_mat - rank_one_crossp[[x]])^2)

            }
          )

          # return the results from the fold
          out <- list(
            estimator = est_name,
            hyperparameters = estimator_hparam,
            loss = 1/nrow(valid_data)*sum(valid_losses),
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
