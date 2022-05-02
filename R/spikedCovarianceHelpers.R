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
#' @examples
#' eig_vals <- rep(1, 50)
#' p_n_ratio <- 0.8
#' estimateNoise(eig_vals, p_n_ratio)
#'
estimateNoise <- function(eig_vals, p_n_ratio) {

  # extract the median of the Marceko-Pastur distribution
  median_mp_dist <- RMTstat::qmp(p = 0.5, svr = 1/p_n_ratio)

  # estimate the noise
  noise <- median(eig_vals) / median_mp_dist

  return(noise)
}


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
