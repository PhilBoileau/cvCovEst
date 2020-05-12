#' Linear Shrinkage Estimator of Covariance Matrix
#'
#' @description \code{linearShrinkEst} computes the linear shrinkage estimate
#'   of the covariance matrix for a given value of \code{alpha}. The linear
#'   shrinkage estimator is defined as the convex combination of the sample
#'   covariance matrix and the identity matrix. The choice of \code{alpha}
#'   determines the bias-variance tradeoff of the estimators in this class:
#'   values near 1 are more likely to exhibit high variance but low bias, and
#'   values near 0 are more likely to be be very biased but have low variance.
#'
#' @param dat A numeric \code{data.frame} or matrix.
#' @param alpha A \code{numeric} between 0 and 1 defining the convex combination
#'   of the sample covariance matrix and the identity. \code{alpha = 1} returns
#'   the sample covariance matrix, and \code{alpha = 0} returns the identity.
#'
#' @return A \code{matrix} corresponding to the estimate of the covariance
#'   matrix.
#'
#' @importFrom coop covar
#'
#' @export
linearShrinkEst <- function(dat, alpha) {

  # compute the sample covariance matrix
  sample_cov_mat <- coop::covar(dat)

  # define the identity
  I_pn <- diag(ncol(dat))

  # shrink the sample covariance matrix
  return(alpha*sample_cov_mat + (1-alpha)*I_pn)
}


#' Simple Thresholding Estimator of Covariance Matrix
#'
#' @description \code{thresholdingEst} computes the thresholding estimate of the
#'   covariance matrix for a given value of \code{gamma}. The simple threshold
#'   estimator of convariance matrix applies a hard thresholding operator to
#'   each element of the sample covariance matrix.
#'
#' @param dat A numeric \code{data.frame} or matrix.
#' @param gamma A \code{numeric} larger than or equal to 0 defining the hard
#'   threshold applied to each element of the \code{dat}'s sample covariance
#'   matrix.
#'
#' @return A \code{matrix} corresponding to the estimate of the covariance
#'   matrix.
#'
#' @importFrom coop covar
#'
#' @export
thresholdingEst <- function(dat, gamma) {

  # compute the sample covariance matrix
  sample_cov_mat <- coop::covar(dat)

  # remove all elements smaller than gamma
  return(replace(sample_cov_mat, abs(sample_cov_mat) < gamma, 0))
}


#' Sample Covariance Matrix
#'
#' @description \code{sampleCovEst} computes the sample covariance matrix. This
#'   function is a wrapper around \code{\link[coop]{covar}}.
#'
#' @param dat A numeric \code{data.frame} or matrix.
#'
#' @return A \code{matrix} corresponding to the estimate of the covariance
#'   matrix.
#'
#' @importFrom coop covar
#'
#' @export
sampleCovEst <- function(dat) {

  # compute the same covariance matrix
  return(coop::covar(dat))

}
