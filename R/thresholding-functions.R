#' Smoothly Clipped Absolute Deviation Thresholding Function
#'
#' @param entry A \code{numeric} entry in a covariance matrix estimate.
#' @param lambda A non-negative \code{numeric} defining the amount of
#'   thresholding applied to each element of \code{dat}'s sample covariance
#'   matrix.
#' @param a A \code{numeric} larger than or equal to \code{2} defining the point
#'   at which the SCAD thresholding functions becomes equal to the hard
#'   thresholding function.
#'
#' @return A regularized \code{numeric}.
#'
#' @keywords internal
scadThreshold <- function(entry, lambda, a) {
  if(abs(entry) <= 2*lambda) {
    reg_entry <- abs(entry) - lambda
    if (reg_entry > 0)
      reg_entry <- sign(entry) * reg_entry
    else
      reg_entry <- 0
  } else if (abs(entry) <= a*lambda) {
    reg_entry <- ((a - 1)*entry - sign(entry)*a*lambda)/(a - 2)
  } else {
    reg_entry <- entry
  }

  return(reg_entry)
}
