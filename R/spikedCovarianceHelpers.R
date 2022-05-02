#' Estimate Noise in Spiked Covariance Matrix Model
#'
#' @description \code{estimateNoise()} estimates the unknown noise term in a
#'   Gaussian spiked covariance matrix model, where the covariance matrix is
#'   assumed to be the identity matrix multiplied by the unknown noise, save for
#'   a few "spiked" entries. This procedures is described in
#'   \insertCite{donoho2018;textual}{cvCovEst}.
#'
#' @param eig_vals A \code{numeric} vector of estimated eigenvalues.
#' @param p_n_ratio A \code{numeric} indicating the asymptotic ratio of the
#'   number of features, p, and the number of observations, n. This ratio is
#'   assumed to be between 0 and 1.
#'
#' @return A \code{numeric} estimate of the noise term in a spiked covariance
#'   matrix model.
#'
#' @importFrom RMTstat qmp
#'
#' @references
#'   \insertAllCited{}
#'
#' @keywords internal
#'
estimateNoise <- function(eig_vals, p_n_ratio) {

  # extract the median of the Marceko-Pastur distribution
  median_mp_dist <- RMTstat::qmp(p = 0.5, svr = 1/p_n_ratio)

  # estimate the noise
  noise <- median(eig_vals) / median_mp_dist

  return(noise)
}

################################################################################

#' Extract Estimated Scaled Eigenvalues in Spiked Covariance Matrix Model
#'
#' @description \code{scaleEigVals()} computes the scaled eigenvalues, and
#'   filters out all eigenvalues that do not need to be shrunk.
#'
#' @param eig_vals A \code{numeric} vector of estimated eigenvalues.
#' @param noise  \code{numeric} representing the known scalar multiple of the
#'   identity matrix giving the approximate population covariance matrix.
#' @param p_n_ratio A \code{numeric} indicating the asymptotic ratio of the
#'   number of features, p, and the number of observations, n. This ratio is
#'   assumed to be between 0 and 1.
#' @param num_spikes \code{numeric} integer equal to or larger than one which
#'   providing the known number of spikes in the population covariance matrix.
#'   If set to \code{NULL}, the number of spikes is estimated.
#'
#' @return A \code{numeric} vector of the scaled eigenvalues to be shrunk.
#'
#' @keywords internal
#'
scaleEigVals <- function(eig_vals, noise, p_n_ratio, num_spikes) {

  # scale the eigenvalues
  scaled_eig_vals <- eig_vals / noise

  # filter the scaled eigenvalues
  # when the number of spikes is unknown
  if (is.null(num_spikes)) {

    eig_val_cutoff <- (1 + sqrt(p_n_ratio))^2
    scaled_eig_vals <- scaled_eig_vals[which(scaled_eig_vals > eig_val_cutoff)]

  # when the number of spikes is known
  } else {

    scaled_eig_vals <- scaled_eig_vals[1:num_spikes]

  }

  return(scaled_eig_vals)
}

################################################################################

#' Estimate Ell of Spiked Covariance Matrix Estimator
#'
#' @description \code{computeEll()} computes the ell value described in
#'   \insertCite{donoho2018;textual}{cvCovEst}.
#'
#' @param scaled_eig_vals A \code{numeric} vector of scaled estimated
#'   eigenvalues.
#' @param p A \code{numeric} integer indicating the number of features in the
#'   data.
#' @param p_n_ratio A \code{numeric} indicating the asymptotic ratio of the
#'   number of features, p, and the number of observations, n. This ratio is
#'   assumed to be between 0 and 1.
#'
#' @return A \code{numeric} vector.
#'
#' @references
#'   \insertAllCited{}
#'
#' @keywords internal
#'
computeEll <- function(scaled_eig_vals, p, p_n_ratio) {

  # get the number of spikes
  num_spikes <- length(scaled_eig_vals)
  if (num_spikes > 0) {

    # compute the shrinkage factor for the spiked eigenvalues
    shrink_factor <- (scaled_eig_vals + 1 - p_n_ratio +
                        sqrt((scaled_eig_vals + 1 - p_n_ratio)^2 -
                             4  * scaled_eig_vals)) / 2

    # define ell
    ell <- c(shrink_factor, rep(1, p - num_spikes))

  } else {
    ell <- rep(1, p)
  }

  return(ell)
}

################################################################################

#' Estimate C of Spiked Covariance Matrix Estimator
#'
#' @description \code{computeC()} computes the c(ell) value described in
#'   \insertCite{donoho2018;textual}{cvCovEst}.
#'
#' @param scaled_eig_vals A \code{numeric} vector of scaled estimated
#'   eigenvalues.
#' @param p A \code{numeric} integer indicating the number of features in the
#'   data.
#' @param p_n_ratio A \code{numeric} indicating the asymptotic ratio of the
#'   number of features, p, and the number of observations, n. This ratio is
#'   assumed to be between 0 and 1.
#'
#' @return A \code{numeric} vector.
#'
#' @references
#'   \insertAllCited{}
#'
#' @keywords internal
#'
computeC <- function(ell, p_n_ratio) {

  # get the dimension of the covariance matrix
  p <- length(ell)

  # filter ell
  ell <- ell[which(ell > (1 + sqrt(p_n_ratio)))]

  # compute c for entrie in ell that are large enough
  c_donoho <- sqrt((1 - p_n_ratio / (ell - 1)^2) / (1 + p_n_ratio / (ell - 1)))

  # all other entries are set to 0
  c_donoho <- c(c_donoho, rep(0, p - length(c_donoho)))

  return(c_donoho)
}
