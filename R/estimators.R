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
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#' @param alpha A \code{numeric} between 0 and 1 defining convex combinations
#'  of the sample covariance matrix and the identity. \code{alpha = 1} returns
#'  the sample covariance matrix, and \code{alpha = 0} returns the identity.
#'
#' @importFrom coop covar
#'
#' @return A \code{matrix} corresponding to the estimate of the covariance
#'  matrix.
#'
#' @export
linearShrinkEst <- function(dat, alpha) {
  # compute the sample covariance matrix
  sample_cov_mat <- coop::covar(dat)

  # define the identity
  idn_pn <- diag(ncol(dat))

  # shrink the sample covariance matrix
  return(alpha * sample_cov_mat + (1 - alpha) * idn_pn)
}

################################################################################

#' Ledoit-Wolf Linear Shrinkage Estimator
#'
#' @description \code{linearShrinkLWEst} computes the asymptotically optimal
#'  convex combination of the sample covariance matrix and the identity. This
#'  convex combination effectively shrinks the eigenvalues of the sample
#'  covariance matrix towards the identity. This estimator is more accurate
#'  than the sample covariance matrix in high-dimensional settings under loose
#'  assumptions. For more information, review the manuscript by
#'  \insertRef{Ledoit2004}{cvCovEst}).
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#'
#' @return A \code{matrix} corresponding to the Ledoit-Wolf linear shrinkgage
#'  estimate of the covariance matrix.
#'
#' @importFrom matrixStats sum2
#'
#' @export
linearShrinkLWEst <- function(dat) {

  # get the number of variables and observations
  p_n <- ncol(dat)
  n <- nrow(dat)

  # compute the sample covariance matrix
  sample_cov_mat <- coop::covar(dat)

  # define the identity
  idn_pn <- diag(p_n)

  # estimate the scalers
  dat <- as.matrix(dat)
  m_n <- matrixStats::sum2(sample_cov_mat * idn_pn) / p_n
  d_n_2 <- matrixStats::sum2((sample_cov_mat - m_n*idn_pn)^2) / p_n
  b_bar_n_2 <- apply(dat, 1,
    function(x) {

      matrixStats::sum2((tcrossprod(x)  - sample_cov_mat)^2)

    }
  )
  b_bar_n_2 <- 1/n^2 * 1/p_n * sum(b_bar_n_2)
  b_n_2 <- min(b_bar_n_2, d_n_2)

  # compute the estimator
  return(b_n_2/d_n_2*m_n*idn_pn + (d_n_2 - b_n_2)/d_n_2*sample_cov_mat)

}

################################################################################

#' Simple Thresholding Estimator of Covariance Matrix
#'
#' @description \code{thresholdingEst} computes the thresholding estimate of
#'  the covariance matrix for a given value of \code{gamma}. The threshold
#'  estimator of the covariance matrix applies a hard thresholding operator to
#'  each element of the sample covariance matrix. For more information on this
#'  estimator, review #'  \insertRef{Bickel2008}{cvCovEst}).
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#' @param gamma A \code{numeric} larger than or equal to 0 defining the hard
#'  threshold applied to each element of \code{dat}'s sample covariance matrix.
#'
#' @importFrom coop covar
#'
#' @return A \code{matrix} corresponding to the estimate of the covariance
#'   matrix.
#'
#' @export
thresholdingEst <- function(dat, gamma) {
  # compute the sample covariance matrix
  sample_cov_mat <- coop::covar(dat)

  # apply threshold by removing all elements smaller than gamma
  return(replace(sample_cov_mat, abs(sample_cov_mat) < gamma, 0))
}

################################################################################

#' Sample Covariance Matrix
#'
#' @description \code{sampleCovEst} computes the sample covariance matrix. This
#'   function is a simple wrapper around \code{\link[coop]{covar}}.
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#'
#' @importFrom coop covar
#'
#' @return A \code{matrix} corresponding to the estimate of the covariance
#'   matrix.
#'
#' @export
sampleCovEst <- function(dat) {
  # compute the sample covariance matrix
  return(coop::covar(dat))
}


################################################################################

#' Banding Estimator
#'
#' @description \code{bandingEst} estimates the covariance matrix of a data frame
#'   with ordered variables by forcing off-diagonal entries to be zero for
#'   indicies that are far removed from one another.  The i, j - entry of the
#'   estimated covariance matrix will be zero if the absolute value of i - j is
#'   greater than some non-negative constant, k.  \emph{Note: argument checks for
#'   this function were removed for computational efficiency.}
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#'
#' @param k A non-negative, numeric integer
#'
#' @importFrom coop covar
#'
#' @return A \code{matrix} corresponding to the estimate of the covariance
#'   matrix.
#'
#' @export
bandingEst <- function(dat, k) {

  # compute the sample covariance matrix
  sam_cov <- coop::covar(dat)

  n <- ncol(sam_cov)

  # loop over different indicies to create an indicator matrix
  indicator_list <- lapply(1:n, function(i) {
    # only consider the lower triangular matrix entries
    j <- i:n

    # calculate/indicate any differences greater than k
    di <- ifelse(abs(i - j) > k, 0, 1)

    # create a new vector corresponding to lower triangular matrix
    di <- c(rep(0, i-1), di)

    di

  })

  # combine vectors
  indicator_matrix <- dplyr::bind_cols(indicator_list)

  # flip the matrix
  indicator_matrix <- indicator_matrix + t(indicator_matrix) - diag(1, n)

  # replace the sample covariance matrix
  sam_cov <- replace(sam_cov, which(indicator_matrix == 0), 0)

  return(sam_cov)
}

