#' Smoothly Clipped Absolute Deviation Thresholding Function
#'
#' @description \code{scadThreshold()} applies the smoothly clipped absolute
#'   deviation thresholding function to the entries of a \code{matrix}.
#'   In particular, it is meant to be applied to the sample covariance matrix.
#'
#' @param entry A \code{numeric} entry in a covariance matrix estimate.
#' @param lambda A non-negative \code{numeric} defining the degree of
#'  thresholding applied to each element of \code{dat}'s sample covariance
#'  matrix.
#' @param a A \code{numeric} larger than or equal to \code{2} defining the
#'  point at which the SCAD thresholding functions becomes equal to the hard
#'  thresholding function.
#'
#' @importFrom assertthat assert_that
#'
#' @return A regularized \code{numeric}.
#'
#' @keywords internal
scadThreshold <- function(entry, lambda, a) {
  # size safety of equivalence between SCAD and hard thresholding
  assertthat::assert_that(a >= 2)

  # Vectorized Version
  e1 <- abs(entry) <= 2 * lambda
  e2 <- abs(entry) > 2 * lambda & abs(entry) <= a * lambda

  entry[e1] <- ifelse(
    abs(entry[e1]) - lambda > 0, sign(entry[e1]) * (abs(entry[e1]) - lambda), 0
  )

  entry[e2] <- ((a - 1) * entry[e2] - sign(entry[e2]) * a * lambda) / (a - 2)

  return(entry)
}

###############################################################################

#' Adaptive LASSO Thresholding Function
#'
#' @description \code{adaptiveLassoThreshold()} applies the adaptive LASSO
#'  thresholding function to the entries of a \code{matrix}. In particular, it
#'  is meant to be applied to sample covariance matrix
#'
#' @param entry A \code{numeric} entry in a covariance matrix estimate.
#' @param lambda A non-negative \code{numeric} defining the amount of
#'  thresholding applied to each element of \code{dat}'s sample covariance
#'  matrix.
#' @param n A non-negative \code{numeric} defining the adaptive weight
#'  applied to each element of \code{dat}'s sample covariance matrix.
#'
#' @return A regularized \code{numeric}.
#'
#' @keywords internal
adaptiveLassoThreshold <- function(entry, lambda, n) {
  # define thresholding multiplier
  s <- abs(entry) - (lambda^(n + 1)) * abs(entry)^(-n)

  # apply regularization
  s[s < 0] <- 0
  s[s > 0] <- sign(entry[s > 0]) * s[s > 0]

  # output regularized entry
  return(s)
}
