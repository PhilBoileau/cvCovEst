#' Scaled Frobenius Loss
#'
#' @description \code{scaledFrobeniusLoss} computes the tractable version of the
#'   scaled Frobenius loss for an estimate of the covariance matrix.
#'
#' @param est_mat A \code{matrix} output from an estimator function.
#' @param sample_cov_mat A \code{matrix} consisting of the sample covariance
#'   matrix.
#' @param cov_sum A \code{numeric} representing the sum of covariance estimates
#'   between the identically positioned elements of \code{est_mat} and
#'   \code{sample_cov_mat}.
#'
#' @return The \code{numeric} scaled Frobenius loss of \code{est_mat}.
#'
#' @importFrom matrixStats sum2
#'
#' @keywords internal
scaledFrobeniusLoss <- function(est_mat, sample_cov_mat, cov_sum) {

  # TODO: implement sparse option
  # ASIDE: don't check format here; will be computationally taxing

  # compute the scaling factor
  scaling_factor <- 1/nrow(est_mat)

  # compute the Frobenius norm
  frob_norm <- matrixStats::sum2((est_mat - sample_cov_mat)^2)

  # return the loss
  return(scaling_factor*frob_norm + 2*scaling_factor*cov_sum)

}
