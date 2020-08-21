#' Penalized Scaled Frobenius Loss
#'
#' @description \code{penScaledFrobeniusLoss} computes the the penalized scaled
#'  Frobenius matrix loss of a covariance matrix estimate.
#'
#' @param est_mat A \code{matrix} output from an estimator function.
#' @param sample_cov_mat A \code{matrix}, the sample covariance matrix.
#' @param cov_sum A \code{numeric} representing the sum of covariance estimates
#'  between the identically positioned elements of \code{est_mat} and
#'  \code{sample_cov_mat}.
#'
#' @return The \code{numeric} penalized scaled Frobenius loss of \code{est_mat}.
#'
#' @importFrom matrixStats sum2
#'
#' @keywords internal
penScaledFrobeniusLoss <- function(est_mat, sample_cov_mat, cov_sum) {
  # TODO: implement sparse option
  # NOTE: don't check format here; will be computationally taxing

  # compute the scaling factor
  scaling_factor <- 1 / nrow(est_mat)

  # compute the Frobenius norm
  frobenius_norm <- matrixStats::sum2((est_mat - sample_cov_mat)^2)

  # return the loss
  loss <- scaling_factor * frobenius_norm + 2 * scaling_factor * cov_sum
  return(loss)
}


#' Aggregated Frobenius Loss
#'
#' @description \code{frobeniusLoss} computes the aggregated, cross-validated
#'   loss of \code{est_mat}, an estimate output by a covariance matrix
#'   estimator, over a validation set.
#'
#' @param est_mat A \code{matrix} output from an estimator function.
#' @param valid_data A \code{numeric} \code{matrix} corresponding to the
#'   validation set.
#'
#' @return The \code{numeric} aggregated Frobenius loss of \code{est_mat} for
#'   over \code{valid_data}.
#'
#' @importFrom matrixStats sum2
#' @importFrom Matrix tcrossproc
#'
#' @keywords internal
frobeniusloss <- function(est_mat, valid_data) {
  # TODO: implement sparse option
  # NOTE: don't check format here; will be computationally taxing

  # compute the Frobenius loss for each observation in the validation set
  valid_losses <- apply(
    valid_data, 1,
    function(obs) {

      # compute the rank one estimate of the validation set observation
      rank_one_cov_mat <- Matrix::tcrossprod(obs)

      # compute the Frobenius loss
      matrixStats::sum2((est_mat - rank_one_cov_mat)^2)

    }
  )

  # compute the aggregated Frobenius loss
  return(sum(valid_losses))
}
