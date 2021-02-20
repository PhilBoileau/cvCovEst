#' Linear Shrinkage Estimator
#'
#' @description \code{linearShrinkEst()} computes the linear shrinkage estimate
#'  of the covariance matrix for a given value of \code{alpha}. The linear
#'  shrinkage estimator is defined as the convex combination of the sample
#'  covariance matrix and the identity matrix. The choice of \code{alpha}
#'  determines the bias-variance tradeoff of the estimators in this class:
#'  values near 1 are more likely to exhibit high variance but low bias, and
#'  values near 0 are more likely to be be very biased but have low variance.
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#' @param alpha A \code{numeric} between 0 and 1 defining convex combinations
#'  of the sample covariance matrix and the identity. \code{alpha = 1} produces
#'  the sample covariance matrix, and \code{alpha = 0} returns the identity.
#'
#' @importFrom coop covar
#'
#' @return A \code{matrix} corresponding to the estimate of the covariance
#'  matrix.
#'
#' @examples
#' linearShrinkEst(dat = mtcars, alpha = 0.1)
#' @export
linearShrinkEst <- function(dat, alpha) {
  # compute the sample covariance matrix
  sample_cov_mat <- coop::covar(dat)

  # define the identity
  idn_pn <- diag(ncol(dat))

  # shrink the sample covariance matrix
  estimate <- alpha * sample_cov_mat + (1 - alpha) * idn_pn
  return(estimate)
}

###############################################################################

#' Ledoit-Wolf Linear Shrinkage Estimator
#'
#' @description \code{linearShrinkLWEst()} computes an asymptotically optimal
#'  convex combination of the sample covariance matrix and the identity matrix.
#'  This convex combination effectively shrinks the eigenvalues of the sample
#'  covariance matrix towards the identity. This estimator is more accurate
#'  than the sample covariance matrix in high-dimensional settings under fairly
#'  loose assumptions. For more information, consider reviewing the manuscript
#'  by \insertCite{Ledoit2004;textual}{cvCovEst}.
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#'
#' @importFrom matrixStats sum2
#'
#' @return A \code{matrix} corresponding to the Ledoit-Wolf linear shrinkage
#'  estimate of the covariance matrix.
#'
#' @references
#'  \insertAllCited{}
#'
#' @examples
#' linearShrinkLWEst(dat = mtcars)
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
  d_n_2 <- matrixStats::sum2((sample_cov_mat - m_n * idn_pn)^2) / p_n
  b_bar_n_2 <- apply(
    dat, 1,
    function(x) {
      matrixStats::sum2((tcrossprod(x) - sample_cov_mat)^2)
    }
  )
  b_bar_n_2 <- 1 / n^2 * 1 / p_n * sum(b_bar_n_2)
  b_n_2 <- min(b_bar_n_2, d_n_2)

  # compute the estimator
  estimate <- (b_n_2 / d_n_2) * m_n * idn_pn +
    (d_n_2 - b_n_2) / d_n_2 * sample_cov_mat
  return(estimate)
}

###############################################################################

#' Hard Thresholding Estimator
#'
#' @description \code{thresholdingEst()} computes the hard thresholding estimate
#'  of the covariance matrix for a given value of \code{gamma}. The threshold
#'  estimator of the covariance matrix applies a hard thresholding operator to
#'  each element of the sample covariance matrix. For more information on this
#'  estimator, review \insertCite{Bickel2008_thresh;textual}{cvCovEst}.
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#' @param gamma A non-negative \code{numeric} defining the degree of hard
#'  thresholding applied to each element of \code{dat}'s sample covariance
#'  matrix.
#'
#' @importFrom coop covar
#'
#' @return A \code{matrix} corresponding to the estimate of the covariance
#'  matrix.
#'
#' @references
#'   \insertAllCited{}
#'
#' @examples
#' thresholdingEst(dat = mtcars, gamma = 0.2)
#' @export
thresholdingEst <- function(dat, gamma) {
  # compute the sample covariance matrix
  sample_cov_mat <- coop::covar(dat)

  # apply threshold by removing all elements smaller than gamma
  estimate <- replace(sample_cov_mat, abs(sample_cov_mat) < gamma, 0)
  return(estimate)
}

###############################################################################

#' Sample Covariance Matrix
#'
#' @description \code{sampleCovEst()} computes the sample covariance matrix. This
#'  function is a simple wrapper around \code{\link[coop]{covar}()}.
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#'
#' @importFrom coop covar
#'
#' @return A \code{matrix} corresponding to the estimate of the covariance
#'  matrix.
#'
#' @examples
#' sampleCovEst(dat = mtcars)
#' @export
sampleCovEst <- function(dat) {
  # compute the sample covariance matrix
  estimate <- coop::covar(dat)
  return(estimate)
}

###############################################################################

#' Banding Estimator
#'
#' @description \code{bandingEst()} estimates the covariance matrix of data with
#'  ordered variables by forcing off-diagonal entries to be zero for indices
#'  that are far removed from one another. The {i, j} entry of the estimated
#'  covariance matrix will be zero if the absolute value of {i - j} is greater
#'  than some non-negative constant \code{k}. This estimator was proposed by
#'  \insertCite{bickel2008_banding;textual}{cvCovEst}.
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#' @param k A non-negative, \code{numeric} integer.
#'
#' @importFrom coop covar
#'
#' @return A \code{matrix} corresponding to the estimate of the covariance
#'  matrix.
#'
#' @references
#'  \insertAllCited{}
#'
#' @examples
#' bandingEst(dat = mtcars, k = 2L)
#' @export
bandingEst <- function(dat, k) {
  # compute the sample covariance matrix
  sam_cov <- coop::covar(dat)
  ncol_sampcov <- ncol(sam_cov)

  # loop over different indices to create an indicator matrix
  indicator_list <- lapply(seq_len(ncol_sampcov), function(i) {
    # only consider the lower triangular matrix entries
    j <- seq(i, ncol_sampcov)

    # calculate/indicate any differences greater than k
    di <- ifelse(abs(i - j) > k, 0, 1)

    # create a new vector corresponding to lower triangular matrix
    di <- c(rep(0, i - 1), di)
    return(di)
  })

  # combine vectors
  indicator_matrix <- do.call(cbind, indicator_list)

  # flip the matrix
  indicator_matrix <- indicator_matrix + t(indicator_matrix) -
    diag(1, ncol_sampcov)

  # replace the sample covariance matrix
  sam_cov <- replace(sam_cov, which(indicator_matrix == 0), 0)
  sam_cov <- unname(sam_cov)

  return(sam_cov)
}

###############################################################################

#' Tapering Estimator
#'
#' @description \code{taperingEst()} estimates the covariance matrix of a
#'  \code{data.frame}-like object with ordered variables by gradually shrinking
#'  the bands of the sample covariance matrix towards zero. The estimator is
#'  defined as the Hadamard product of the sample covariance matrix and a
#'  weight matrix. The amount of shrinkage is dictated by the weight matrix
#'  and is specified by a hyperparameter \code{k}. This estimator is attributed
#'  to \insertCite{cai2010;textual}{cvCovEst}.
#'
#'  The weight matrix is a Toeplitz matrix with entries defined as follows. Let
#'  i and j index the rows and columns of the weight matrix, respectively. If
#'  \code{abs(i - j) <= k / 2}, then entry {i, j} in the weight matrix is equal
#'  to 1. If \code{k / 2 < abs(i - j) < k}, then entry {i, j} is equal to
#'  \code{2 - 2 * abs(i - j) / k}. Otherwise, entry {i, j} is equal to 0.
#'
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#' @param k A non-negative, even \code{numeric} integer.
#'
#' @importFrom coop covar
#'
#' @return A \code{matrix} corresponding to the estimate of the covariance
#'  matrix.
#'
#' @references
#'  \insertAllCited{}
#'
#' @examples
#' taperingEst(dat = mtcars, k = 0.1)
#' @export
taperingEst <- function(dat, k) {
  # compute the sample covariance matrix
  sam_cov <- coop::covar(dat)
  ncol_sampcov <- ncol(sam_cov)
  k_h <- k / 2

  # loop over different indices to create weight vectors
  weight_list <- lapply(seq_len(ncol_sampcov), function(i) {
    # only consider the lower triangular matrix entries
    j <- seq(i, ncol_sampcov)

    # calculate the difference in indices
    di <- abs(i - j)

    # loop over elements in the difference vector and assign weights
    w <- sapply(di, function(d) {
      if (d <= k_h) {
        wi <- 1
      } else if (d > k_h & d < k) {
        wi <- 2 - (d / k_h)
      } else {
        wi <- 0
      }
      return(wi)
    })

    # multiply by corresponding entries in sample covariance matrix
    sam_vec <- sam_cov[j, i]
    sam_vec <- sam_vec * w

    # create a new vector corresponding to lower triangular matrix column
    sam_vec <- c(rep(0, i - 1), sam_vec)
    return(sam_vec)
  })

  # combine vectors
  weight_matrix <- do.call(cbind, weight_list)

  # flip the matrix
  weight_matrix <- weight_matrix + t(weight_matrix) - diag(diag(sam_cov))

  # return the new weight matrix
  return(unname(weight_matrix))
}

###############################################################################

#' Analytical Non-Linear Shrinkage Estimator
#'
#' @description \code{nlShrinkLWEst()} invokes the analytical estimator
#'  presented by \insertCite{Ledoit2020;textual}{cvCovEst} for applying a
#'  nonlinear shrinkage function to the sample eigenvalues of the covariance
#'  matrix. The shrinkage function relies on an application of the Hilbert
#'  Transform to an estimate of the sample eigenvalues' limiting spectral
#'  density. This estimated density is computed with the Epanechnikov kernel
#'  using a global bandwidth parameter of \code{n^(-1/3)}. The resulting
#'  shrinkage function pulls eigenvalues towards the nearest mode of their
#'  empirical distribution, thus creating a localized shrinkage effect rather
#'  than a global one.
#'
#'  We do not recommend that this estimator be employed when
#'  the estimand is the correlation matrix. The diagonal entries of the
#'  resulting estimate are not guaranteed to be equal to one.
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#'
#' @importFrom coop covar
#' @importFrom matrixStats rowMeans2
#'
#' @return A \code{matrix} corresponding to the estimate of the covariance
#'  matrix.
#'
#' @references
#'   \insertAllCited{}
#'
#' @examples
#' nlShrinkLWEst(dat = mtcars)
#' @export
nlShrinkLWEst <- function(dat) {
  # get the dimensions of the data
  n <- nrow(dat)
  p <- ncol(dat)

  # Compute the Sample Covariance Matrix
  sam_cov <- coop::covar(dat)

  # Get the sorted eigenvalues and eigenvectors
  sam_eig <- eigen(sam_cov, symmetric = TRUE)
  lambda <- sam_eig$values
  u <- sam_eig$vectors

  # Analytical Nonlinear Shrinkage Kernal Formula
  # accept a tolerance of 1e-2
  eig_nonzero_tol <- sum(lambda > 0.01)
  i <- min(p, eig_nonzero_tol)
  lambda <- lambda[seq_len(i)]
  r <- length(lambda)
  L <- matrix(lambda, nrow = r, ncol = i)

  # LW Equation 4.9
  h <- n^(-1 / 3)
  H <- h * t(L)
  x <- (L - t(L)) / H

  # LW Equation 4.7
  s1 <- (3 / 4) / sqrt(5)
  s2 <- -(3 / 10) / pi
  pos_x <- (1 - (x^2) / 5)
  pos_x <- replace(pos_x, pos_x < 0, 0)
  f_tilde <- s1 * matrixStats::rowMeans2(pos_x / H)

  # LW Equation 4.8
  log_term <- log(abs((sqrt(5) - x) / (sqrt(5) + x)))
  Hftemp <- (s2 * x) + (s1 / pi) * (1 - (x^2) / 5) * log_term
  sq5 <- which(abs(x) == sqrt(5))
  Hftemp[sq5] <- s2 * x[sq5]
  H_tilde <- matrixStats::rowMeans2(Hftemp / H)

  # LW Equation 4.3
  s3 <- pi * (p / n)
  s4 <- 1 / (h^2)
  if (p <= eig_nonzero_tol) {
    d_tilde <- lambda / ((s3 * lambda * f_tilde)^2 +
      (1 - (p / n) - s3 * lambda * H_tilde)^2)
  } else {
    ones <- rep(1, p - eig_nonzero_tol)
    log_term <- log((1 + sqrt(5) * h) / (1 - sqrt(5) * h))
    m <- mean(1 / lambda)

    # LW Equation C.8
    Hf_tilde0 <- (1 / pi) * ((3 / 10) * s4 +
      (s1 / h) * (1 - (1 / 5) * s4) * log_term) * m

    # LW Equation C.5
    d_tilde0 <- 1 / (pi * (p - eig_nonzero_tol) / eig_nonzero_tol * Hf_tilde0)

    # LW Equation C.4
    d_tilde1 <- lambda / ((pi^2 * lambda^2) * (f_tilde^2 + H_tilde^2))
    d_tilde <- c(d_tilde1, d_tilde0 * ones)
  }

  # LW Equation 4.4
  sigma_tilde <- u %*% tcrossprod(diag(d_tilde), u)
  return(sigma_tilde)
}

###############################################################################

#' Linear Shrinkage Estimator, Dense Target
#'
#' @description \code{denseLinearShrinkEst()} computes the asymptotically
#'  optimal convex combination of the sample covariance matrix and a dense
#'  target matrix. This target matrix's diagonal elements are equal to the
#'  average of the sample covariance matrix estimate's diagonal elements, and
#'  its off-diagonal elements are equal to the average of the sample covariance
#'  matrix estimate's off-diagonal elements. For information on this
#'  estimator's derivation, see \insertCite{Ledoit2020b;textual}{cvCovEst} and
#'  \insertCite{shafer2005;textual}{cvCovEst}.
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#'
#' @importFrom coop covar
#' @importFrom stats var
#' @importFrom Matrix triu
#' @importFrom matrixStats sum2
#'
#' @return A \code{matrix} corresponding to the estimate of the covariance
#'  matrix.
#'
#' @references
#'   \insertAllCited{}
#'
#' @examples
#' denseLinearShrinkEst(dat = mtcars)
#' @export
denseLinearShrinkEst <- function(dat) {
  # get the number of variables and observations
  p_n <- ncol(dat)

  # compute the sample covariance matrix
  sample_cov_mat <- coop::covar(dat)

  # compute elements of the dense target
  samp_diag <- diag(sample_cov_mat)
  samp_off_diag <- rep(
    Matrix::triu(sample_cov_mat)[upper.tri(sample_cov_mat)],
    each = 2
  )
  mean_var <- mean(samp_diag)
  mean_cov <- mean(samp_off_diag)
  f_mat <- matrix(data = mean_cov, nrow = p_n, ncol = p_n)
  diag(f_mat) <- mean_var

  # compute shrinkage factor
  numerator <- sum(stats::var(diag(samp_diag))) +
    sum(stats::var(samp_off_diag))
  denominator <- matrixStats::sum2((f_mat - sample_cov_mat)^2)
  shrink_factor <- min(1, max(0, numerator / denominator))

  # compute the estimate
  return(shrink_factor * f_mat + (1 - shrink_factor) * sample_cov_mat)
}

###############################################################################

#' Smoothly Clipped Absolute Deviation Estimator
#'
#' @description \code{scadEst()} applies the SCAD thresholding function of
#'  \insertCite{fan2001;textual}{cvCovEst} to each entry of the sample
#'  covariance matrix. This penalized estimator constitutes a compromise
#'  between hard and soft thresholding of the sample covariance matrix: it is
#'  a linear interpolation between soft thresholding up to \code{2 * lambda}
#'  and hard thresholding after \code{3.7 * lambda}
#'  \insertCite{rothman2009}{cvCovEst}.
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#' @param lambda A non-negative \code{numeric} defining the degree of
#'  thresholding applied to each element of \code{dat}'s sample covariance
#'  matrix.
#'
#' @importFrom coop covar
#'
#' @return A \code{matrix} corresponding to the estimate of the covariance
#'  matrix.
#'
#' @references
#'   \insertAllCited{}
#'
#' @examples
#' scadEst(dat = mtcars, lambda = 0.2)
#' @export
scadEst <- function(dat, lambda) {
  # compute the sample covariance matrix
  sample_cov_mat <- coop::covar(dat)

  # apply threshold by removing all elements smaller than gamma
  # TODO: Create a symmertric apply for covariance matrices
  estimate <- apply(sample_cov_mat, c(1, 2), scadThreshold,
    lambda = lambda, a = 3.7
  )
  return(estimate)
}

###############################################################################

#' POET Estimator
#'
#' @description \code{poetEst()} implements the Principal Orthogonal complEment
#'  Thresholding (POET) estimator, a nonparametric, unobserved-factor-based
#'  estimator of the covariance matrix \insertCite{fan2013}{cvCovEst}. The
#'  estimator is defined as the sum of the sample covariance matrix'
#'  rank-\code{k} approximation and its post-thresholding principal orthogonal
#'  complement. The hard thresholding function is used here, though others
#'  could be used instead.
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#' @param k An \code{integer} indicating the number of unobserved latent
#'  factors. Empirical evidence suggests that the POET estimator is robust to
#'  overestimation of this hyperparameter \insertCite{fan2013}{cvCovEst}. In
#'  practice, it is therefore preferable to use larger values.
#' @param lambda A non-negative \code{numeric} defining the amount of
#'  thresholding applied to each element of sample covariance matrix's
#'  orthogonal complement.
#'
#' @importFrom coop covar
#' @importFrom RSpectra eigs_sym
#'
#' @return A \code{matrix} corresponding to the estimate of the covariance
#'  matrix.
#'
#' @references
#'  \insertAllCited{}
#'
#' @examples
#' poetEst(dat = mtcars, k = 2L, lambda = 0.1)
#' @export
poetEst <- function(dat, k, lambda) {
  # compute the sample covariance matrix
  sample_cov_mat <- coop::covar(dat)

  # perform the eigenvalue decomposition
  eig_decomp <- RSpectra::eigs_sym(sample_cov_mat, k)

  # compute spectral decomposition
  spectral_decomp <- lapply(
    seq_len(k),
    function(i) {
      eig_decomp$values[i] * tcrossprod(eig_decomp$vectors[, i])
    }
  )
  spectral_decomp <- Reduce(`+`, spectral_decomp)

  # compute principal orthogonal complement
  poc <- sample_cov_mat - spectral_decomp
  poc_diag <- diag(diag(poc))

  # regularize principal orthogonal complement
  # TODO: Create a symmetric apply for covariance matrices
  poc <- replace(poc, abs(poc) < lambda, 0)
  poc <- poc - diag(diag(poc)) + poc_diag

  # return the estimate
  estimate <- spectral_decomp + poc
  return(estimate)
}

###############################################################################

#' Robust POET Estimator for Elliptical Distributions
#'
#' @description \code{robustPoetEst()} implements the robust version of
#'  Principal Orthogonal complEment Thresholding (POET) estimator, a
#'  nonparametric, unobserved-factor-based estimator of the covariance matrix
#'  when the underlying distribution is ellipitcal
#'  \insertCite{fan2018}{cvCovEst}. The estimator is defined as the sum of the
#'  sample covariance matrix's rank-\code{k} approximation and its
#'  post-thresholding principal orthogonal complement. The rank-\code{k}
#'  approximation is constructed from the sample covariance matrix, its leading
#'  eigenvalues, and its leading eigenvectors.  The sample covariance matrix and
#'  leading eigenvalues are initially estimated via an M-estimation procedure
#'  and the marginal Kendall's tau estimator. The leading eigenvectors are
#'  estimated using spatial Kendall's tau estimator. The hard thresholding
#'  function is used to regularize the idiosyncratic errors' estimated
#'  covariance matrix, though other regularization schemes could be used.
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#' @param k An \code{integer} indicating the number of unobserved latent
#'  factors. Empirical evidence suggests that the POET estimator is robust to
#'  overestimation of this hyperparameter \insertCite{fan2013}{cvCovEst}. In
#'  practice, it is therefore preferable to use larger values.
#' @param lambda A non-negative \code{numeric} defining the amount of
#'  thresholding applied to each element of sample covariance matrix's
#'  orthogonal complement.
#' @param var_est A \code{character} dictating which variance estimator to
#'  use. This must be one of the strings \code{"sample"}, \code{"mad"}, or
#'  \code{"huber"}. \code{"sample"} uses sample variances; \code{"mad"}
#'  estimates variances via median absolute deviation; \code{"huber"} uses an
#'  M-estimator for variance under the Huber loss.
#'
#' @importFrom RSpectra eigs_sym
#' @importFrom matrixStats colSds colMads colVars
#' @importFrom stats optimize
#'
#' @return A \code{matrix} corresponding to the estimate of the covariance
#'   matrix.
#'
#' @references
#'   \insertAllCited{}
#'
#' @examples
#' robustPoetEst(dat = mtcars, k = 2L, lambda = 0.1, var_est = "sample")
#' @export
robustPoetEst <- function(dat, k, lambda,
                          var_est = c("sample", "mad", "huber")) {
  # if data frame, coerce to matrix
  if (!is.matrix(dat)) {
    dat <- as.matrix(dat)
  }

  # set default base covariance estimator
  var_est <- match.arg(var_est)

  # get the dimensions of the data
  n <- nrow(dat)

  # use M-estimator and Huber loss to robustly estimate variance
  if (var_est == "sample") {
    D_est <- diag(matrixStats::colSds(dat))
  } else if (var_est == "mad") {
    D_est <- diag(matrixStats::colMads(dat))
  } else if (var_est == "huber") {
    # This method is proposed by Fan et al but most computationally expensive
    alpha <- sqrt(1 / (8 * max(matrixStats::colVars(dat))))
    huber <- function(x, alpha) {
      if (abs(x) > 1 / alpha) {
        return(2 / alpha * abs(x) - 1 / alpha^2)
      } else {
        return(x^2)
      }
    }
    mest <- function(y) {
      stats::optimize(f = function(x, alpha) {
        sum(sapply(x - y, FUN = huber, alpha = alpha))
      }, alpha = alpha, lower = min(y), upper = max(y))$minimum
    }
    D_est <- diag(sqrt(pmax(apply(dat^2, 2, mest) -
      apply(dat, 2, mest)^2, 1e-6)))
  }

  # Marginal Kendall's tau estimator can be vectorized as the multiplication of
  # the matrix of signs of elementwise differences for each variable
  # with its transpose.
  diff_mat <- function(x) {
    outer(x, x, FUN = "-")[lower.tri(sign(outer(x, x, FUN = "-")))]
  }
  Diff <- apply(dat, 2, diff_mat)

  # calculate the estimator of R
  R_est <- sin(crossprod(sign(Diff)) * (pi / (n * (n - 1))))

  # calculate the first estimator for covariance matrix
  # NOTE: D_est %*% R_est %*% D_est is equivalent but (slightly) slower than
  Sigma_est1 <- tcrossprod(crossprod(D_est, R_est), D_est)

  # calculate the estimator of leading eigenvalues
  if (k == 1) {
    Lambda_est <- as.matrix(RSpectra::eigs_sym(Sigma_est1, k)$values)
  } else {
    Lambda_est <- diag(RSpectra::eigs_sym(Sigma_est1, k)$values)
  }

  # calculate spatial Kendall's tau estimator
  Sigma_est2 <- (2 / (n * (n - 1))) *
    crossprod(Diff / apply(Diff, 1, function(x) sqrt(sum(x^2))))

  # calculate the estimator for leading eigenvectors
  Gamma_est <- RSpectra::eigs_sym(Sigma_est2, k)$vectors

  # calculate the low rank structure
  # NOTE: tcrossprod is faster than manually computing the transpose, as in
  #       Gamma_est %*% Lambda_est %*% t(Gamma_est)
  Low_rank_est <- Gamma_est %*% tcrossprod(Lambda_est, Gamma_est)

  # regularize the principal orthogonal component
  Sigma_estu <- Sigma_est1 - Low_rank_est
  Sigma_estu <- replace(Sigma_estu, abs(Sigma_estu) < lambda, 0)
  estimate <- Sigma_estu + Low_rank_est
  return(estimate)
}

###############################################################################

#' Adaptive LASSO Estimator
#'
#' @description \code{adaptiveLassoEst()} applied the adaptive LASSO to the
#'  entries of the sample covariance matrix. The thresholding function is
#'  inspired by the penalized regression introduced by
#'  \insertCite{zou2006;textual}{cvCovEst}. The thresholding function assigns
#'  a weight to each entry of the sample covariance matrix based on its
#'  initial value. This weight then determines the relative size of the penalty
#'  resulting in larger values being penalized less and reducing bias
#'  \insertCite{rothman2009}{cvCovEst}.
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#' @param lambda A non-negative \code{numeric} defining the amount of
#'  thresholding applied to each element of \code{dat}'s sample covariance
#'  matrix.
#' @param n A non-negative \code{numeric} defining the exponent of the adaptive
#'  weight applied to each element of \code{dat}'s sample covariance matrix.
#'
#' @importFrom coop covar
#'
#' @return A \code{matrix} corresponding to the estimate of the covariance
#'  matrix.
#'
#' @references
#'   \insertAllCited{}
#'
#' @examples
#' adaptiveLassoEst(dat = mtcars, lambda = 0.9, n = 0.9)
#' @export
adaptiveLassoEst <- function(dat, lambda, n) {
  # compute the sample covariance matrix
  sample_cov_mat <- coop::covar(dat)

  # apply adaptive thresholding to the sample covariance matrix
  adaptive_cov_mat <- apply(
    sample_cov_mat, c(1, 2),
    adaptiveLassoThreshold,
    lambda = lambda, n = n
  )

  # output the post-thresholding estimate
  return(adaptive_cov_mat)
}
