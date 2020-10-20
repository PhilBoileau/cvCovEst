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
#' @importFrom dplyr bind_rows
#' @importFrom origami training validation fold_index
#' @importFrom Rdpack reprompt
#' @importFrom tibble tibble
#' @importFrom rlang eval_tidy !!! exec quo_get_expr
#' @importFrom matrixStats sum2
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

  # compute the sum of element-wise squared outer products
  elwise_sq_list <- lapply(rank_one_crossp, `^`, 2)
  elwise_sq_sum <- Reduce(`+`, elwise_sq_list)
  elwise_sq <- matrixStats::sum2(elwise_sq_sum)

  # compute the sum of cross_products
  cross_prod <- Reduce(`+`, rank_one_crossp)

  # get number of estimators
  estimator_funs <- rlang::quo_get_expr(estimator_funs)
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

      # get the hadamard product and sum
      had_crossprod <- matrixStats::sum2(cross_prod * est_mat)

      # get the elementwise square and sum it
      est_square <- matrixStats::sum2(est_mat^2)

      # return the results from the fold
      if (is.null(true_cov_mat)) {
        out <- tibble::tibble(
          estimator = est_name,
          hyperparameters = estimator_hparam,
          loss = 1 / nrow(valid_data) * (elwise_sq - 2 * had_crossprod) +
            est_square,
          fold = origami::fold_index(fold = fold)
        )
      } else {

        # fit the covariance matrix estimator on the full dataset
        # NOTE: this is located here out of convenience... not computationally
        # efficient, but not all that important since only for simulations,
        # and this will not be released in the main package.
        est_mat_full <- est_fun(dat)

        out <- tibble::tibble(
          estimator = est_name,
          hyperparameters = estimator_hparam,
          loss = 1 / nrow(valid_data) * (elwise_sq - 2 * had_crossprod) +
            est_square,
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
          if (length(estimator_hparam) > 1) {
            estimator_hparam <- paste(estimator_hparam, collapse = ", ")
          }
          est_mat <- rlang::exec(
            est_fun,
            train_data,
            !!!as.list(unlist(hparam_grid[idx, ]))
          )

          # get the hadamard product and sum
          had_crossprod <- matrixStats::sum2(cross_prod * est_mat)

          # get the elementwise square and sum it
          est_square <- matrixStats::sum2(est_mat^2)

          # return the results from the fold
          if (is.null(true_cov_mat)) {
            out <- tibble::tibble(
              estimator = est_name,
              hyperparameters = estimator_hparam,
              loss = 1 / nrow(valid_data) * (elwise_sq - 2 * had_crossprod) +
                est_square,
              fold = origami::fold_index(fold = fold)
            )
          } else {

            # fit the covariance matrix estimator on the full dataset
            # NOTE: this is located here out of convenience... not computationally
            # efficient, but not all that important since only for simulations,
            # and this will not be released in the main package.
            est_mat_full <- rlang::exec(
              est_fun,
              dat,
              !!!as.list(unlist(hparam_grid[idx, ]))
            )

            out <- tibble::tibble(
              estimator = est_name,
              hyperparameters = estimator_hparam,
              loss = 1 / nrow(valid_data) * (elwise_sq - 2 * had_crossprod) +
                est_square,
              true_loss = trueFrobeniusLoss(est_mat, true_cov_mat),
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
#' @importFrom matrixStats sum2
#'
#' @keywords internal
trueFrobeniusLoss <- function(estimate, true_covar) {

  # compute the matrix of the cross products of variances
  diag_true_covar <- diag(true_covar)
  cross_prod_mat <- matrixStats::sum2(base::tcrossprod(
    diag_true_covar,
    diag_true_covar
  ))

  # compute the element-wise square of the true covariance matrix
  elem_square_true <- matrixStats::sum2(true_covar^2)

  # compute the element wise multiplication of the estimate and true covariance
  elem_mult <- matrixStats::sum2(estimate * true_covar)

  # compute the element-wise square of the estimate
  elem_square_est <- matrixStats::sum2(estimate^2)

  # combine all loss entries
  loss <- cross_prod_mat + 2 * elem_square_true -
    2 * elem_mult + elem_square_est

  return(loss)
}
################################################################################

#' Cross-Validation Function for Matrix Frobenius Loss
#'
#' @description \code{cvMatrixFrobeniusLoss} evaluates the Matrix Frobenius loss
#'   over a \code{fold} object (from \pkg{origami}
#'   \insertCite{Coyle2018}{cvCovEst}). This loss function is equivalent to that
#'   presented in \code{\link[cvCovEst]{cvFrobeniusLoss}} in terms of estimator
#'   selections, but is more computationally efficient.
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
#' @importFrom dplyr bind_rows
#' @importFrom origami training validation fold_index
#' @importFrom tibble tibble
#' @importFrom rlang eval_tidy !!! exec quo_get_expr
#' @importFrom matrixStats sum2
#'
#' @return A \code{\link[tibble]{tibble}} providing information on estimators,
#'  their hyperparameters (if any), and their scaled Matrix Frobenius loss
#'  evaluated on a given \code{fold}.
#'
#' @references
#'   \insertAllCited{}
#'
#' @export
cvMatrixFrobeniusLoss <- function(fold, dat, estimator_funs,
                                  estimator_params = NULL,
                                  true_cov_mat = NULL) {

  # split the data into training and validation
  train_data <- origami::training(dat)
  valid_data <- origami::validation(dat)

  # compute the sample covariance of the validation set
  sample_cov_valid <- coop::covar(valid_data)

  # get number of estimators
  estimator_funs <- rlang::quo_get_expr(estimator_funs)
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
      out <- tibble::tibble(
        estimator = est_name,
        hyperparameters = estimator_hparam,
        loss = matrixStats::sum2((est_mat - sample_cov_valid)^2),
        fold = origami::fold_index(fold = fold)
      )
    } else {

      # Compute the grid of hyperparameters
      hparam_grid <- expand.grid(estimator_params[[est_name]])

      # loop through the estimator hyperparameters
      param_out <- lapply(
        seq_len(nrow(hparam_grid)),
        function(idx) {

          # fit the covariance matrix estimator on the training set
          estimator_hparam <- paste(hyp_name, "=", hparam_grid[idx, ])
          if (length(estimator_hparam) > 1) {
            estimator_hparam <- paste(estimator_hparam, collapse = ", ")
          }
          est_mat <- rlang::exec(
            est_fun,
            train_data,
            !!!as.list(unlist(hparam_grid[idx, ]))
          )

          # return the results from the fold
          out <- tibble::tibble(
            estimator = est_name,
            hyperparameters = estimator_hparam,
            loss = matrixStats::sum2((est_mat - sample_cov_valid)^2),
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
