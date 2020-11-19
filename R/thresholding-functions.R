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
  if (abs(entry) <= 2 * lambda) {
    reg_entry <- abs(entry) - lambda
    if (reg_entry > 0) {
      reg_entry <- sign(entry) * reg_entry
    } else {
      reg_entry <- 0
    }
  } else if (abs(entry) <= a * lambda) {
    reg_entry <- ((a - 1) * entry - sign(entry) * a * lambda) / (a - 2)
  } else {
    reg_entry <- entry
  }

  return(reg_entry)
}

################################################################################

#' Adaptive LASSO Thresholding Function
#'
#' @param entry A \code{numeric} entry in a covariance matrix estimate.
#' @param lambda A non-negative \code{numeric} defining the amount of
#'   thresholding applied to each element of \code{dat}'s sample covariance
#'   matrix.
#' @param n A non-negative \code{numeric} defining the adaptive weight
#'   applied to each element of \code{dat}'s sample covariance matrix.
#'
#' @return A regularized \code{numeric}.
#'
#' @keywords internal
adaptiveLassoThreshold <- function(entry, lambda, n) {
  s <- abs(entry) - (lambda^(n + 1)) * abs(entry)^(-n)

  if (s > 0) {
    reg_entry <- sign(entry) * s
  } else {
    reg_entry <- 0
  }

  return(reg_entry)
}

################################################################################

#' Symmetric Apply Function for Covariance Matrices
#'
#' @description \code{symmetricApply} implements a version of \code{apply}
#'   expressly made for looping over symmetric matrices.
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#' @param sym_fun a function to apply to each element of the covariance matrix.
#' @param sym_args a named vector or \code{list} of arguments to be passed to
#'   \code{sym_fun}.
#'
#' @return A \code{matrix}.
#'
#' @importFrom rlang exec
#'
#' @keywords internal
symmetricApply <- function(dat, sym_fun, sym_args) {

  # get number of columns
  n <- ncol(dat)

  # loop over different columns in dat
  lower_matrix <- sapply(1:n, function(i) {
    # extract upper triangular entries of dat
    lt_vec <- dat[i, i:n]

    # apply function to each element
    app_vec <- sapply(lt_vec, function(k) {
      f_args <- as.list(c(k, sym_args))
      return(rlang::exec(sym_fun, !!!f_args))
    })

    # return a new vector corresponding to lower triangular matrix column
    new_vec <- c(rep(0, i - 1), app_vec)

    return(new_vec)
  })

  # combine vectors
  sym_matrix <- suppressMessages(dplyr::bind_cols(lower_matrix))
  sym_matrix <- as.matrix(sym_matrix)

  # flip the matrix
  sym_matrix <- sym_matrix + t(sym_matrix) - diag(diag(sym_matrix))

  # rename the columns and rows
  colnames(sym_matrix) <- colnames(dat)
  rownames(sym_matrix) <- colnames(dat)

  # return the new symmetric matrix
  return(sym_matrix)
}
