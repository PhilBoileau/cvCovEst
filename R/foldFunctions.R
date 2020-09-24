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
#' @param true_cov_mat A \code{matrix} like object representing the true
#'  covariance matrix of the data generating distribution. This parameter
#'  is intended for use only in simulation studies, and defaults to a value
#'  of \code{NULL}. If not null, \code{cvFrobeniusLoss} returns the true
#'  cross-validated Frobenius loss in addition to the cross-validated Frobenius
#'  loss estimate for the given \code{fold}.
#'
#' @importFrom coop covar
#' @importFrom dplyr bind_rows
#' @importFrom origami training validation fold_index
#' @importFrom Rdpack reprompt
#' @importFrom tibble tibble
#' @importFrom rlang eval_tidy !!! exec
#'
#' @return A \code{\link[tibble]{tibble}} providing information on estimators,
#'  their hyperparameters (if any), and their scaled Frobenius loss evaluated
#'  on a given \code{fold}.
#'
#' @references
#'   \insertAllCited{}
#'
#' @export
cvFrobeniusLoss <- function(fold, dat, estimator_funs,
                            estimator_params = NULL,
                            true_cov_mat = NULL) {

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
    est_fun <- eval_tidy(estimator_funs[[x]])
    est_name <- as.character(estimator_funs[[x]])

    # check if a hyperparameter is needed
    hyp_name <- names(estimator_params[[est_name]])
    if (is.null(hyp_name)) {

      # fit the covariance matrix estimator on the training set
      est_mat <- est_fun(train_data)

      # indicate that there are no hyperparameters
      estimator_hparam <- "hyperparameters = NA"

      # return the results from the fold
      if (is.null(true_cov_mat)) {
        out <- tibble::tibble(
          estimator = est_name,
          hyperparameters = estimator_hparam,
          loss = frobeniusLoss(est_mat, valid_data, rank_one_crossp),
          fold = origami::fold_index(fold = fold)
        )
      } else {
        out <- tibble::tibble(
          estimator = est_name,
          hyperparameters = estimator_hparam,
          loss = frobeniusLoss(est_mat, valid_data, rank_one_crossp),
          true_loss = trueFrobeniusLoss(est_mat, true_cov_mat),
          fold = origami::fold_index(fold = fold)
        )
      }

      return(out)

    } else {

      # Compute the grid of hyperparameters
      hparam_grid <- expand.grid(estimator_params[[est_name]])

      # loop through the estimator hyperparameters
      param_out <- lapply(
        seq_len(nrow(hparam_grid)),
        function(idx) {

          # fit the covariance matrix estimator on the training set
          estimator_hparam <- paste(hyp_name, "=", hparam_grid[idx, ])
          if (length(estimator_hparam) > 1)
            estimator_hparam <- paste(estimator_hparam, collapse = ", ")
          est_mat <- rlang::exec(
            est_fun,
            train_data,
            !!!as.list(unlist(hparam_grid[idx, ]))
          )

          # fit the covariance matrix estimator on the full dataset
          # NOTE: this is located here out of convenience... not computationally
          # efficient, but not all that important since only for simulations,
          # and this will not be released in the main package.
          est_mat_full <- rlang::exec(
            est_fun,
            dat,
            !!!as.list(unlist(hparam_grid[idx, ]))
          )

          # return the results from the fold
          if (is.null(true_cov_mat)) {
            out <- tibble::tibble(
              estimator = est_name,
              hyperparameters = estimator_hparam,
              loss = frobeniusLoss(est_mat, valid_data, rank_one_crossp),
              fold = origami::fold_index(fold = fold)
            )
          } else {
            out <- tibble::tibble(
              estimator = est_name,
              hyperparameters = estimator_hparam,
              loss = frobeniusLoss(est_mat, valid_data, rank_one_crossp),
              true_loss = trueFrobeniusLoss(est_mat, true_cov_mat),
              true_full_loss = trueFrobeniusLoss(est_mat_full, true_cov_mat),
              fold = origami::fold_index(fold = fold)
            )
          }

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

################################################################################

#' Frobenius Loss
#'
#' @description \code{frobeniusLoss} computes the Frobenius
#'   loss over the given dataset for the given covariance matrix estimate.
#'
#' @param estimate A \code{matrix} corresponding to the estimate of the
#'   covariance matrix.
#' @param dat A \code{data.frame} containing the data on which the loss will be
#'   computed.
#' @param rank_one_crossp A \code{list} of \code{matrix} objects consisting
#'   of the outer products of each of \code{dat}'s rows.
#'
#' @return The average Frobenius loss over \code{dat} of \code{estimate} as
#'   a \code{numeric}.
#'
#' @keywords internal
frobeniusLoss <- function(estimate, dat, rank_one_crossp) {

  # compute the losses
  losses <- sapply(
    seq_len(nrow(dat)),
    function(x) {
      # compute the Frobenius loss
      matrixStats::sum2((estimate - rank_one_crossp[[x]])^2)
    }
  )

  # compute the fold loss
  return(1/nrow(dat)*sum(losses))
}

################################################################################

#' True Cross-Validated Frobenius Loss
#'
#' @description \code{trueFrobeniusLoss} computes the true cross-validated
#'   Frobenius loss over the validation dataset.
#'
#' @param estimate A \code{matrix} corresponding to the estimate of the
#'   covariance matrix.
#' @param true_covar A \code{matrix} corresponding to the true covariance matrix
#'   of the data generating distribution.
#'
#' @return The true average Frobenius loss over the validation dataset of
#'   \code{estimate} as a \code{numeric}.
#'
#' @import Matrix
#'
#' @keywords internal
trueFrobeniusLoss <- function(estimate, true_covar) {

  # compute the matrix of the cross products of variances
  diag_true_covar <- diag(true_covar)
  cross_prod_mat <- Matrix::tcrossprod(diag_true_covar, diag_true_covar)

  # compute the element-wise square of the true covariance matrix
  elem_square_true <- true_covar^2

  # compute the element wise multiplication of the estimate and true covariance
  elem_mult <- estimate * true_covar

  # compute the element-wise square of the estimate
  elem_square_est <- estimate^2

  # combine all loss entries
  combo_mat <- cross_prod_mat + 2*elem_square_true -
    2*elem_mult + elem_square_est

  # compute the true loss
  loss <- matrixStats::sum2(combo_mat)

  return(loss)
}
