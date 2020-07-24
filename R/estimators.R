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
#'  estimator, review \insertRef{Bickel2008_thresh}{cvCovEst}.
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
#'   greater than some non-negative constant, \code{k}. This estimator was
#'   put forth by \insertRef{bickel2008_banding}{cvCovEst}.
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#'
#' @param k A non-negative, numeric integer
#'
#' @importFrom coop covar
#' @importFrom dplyr bind_cols
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

    return(di)

  })

  # combine vectors
  indicator_matrix <- suppressMessages(dplyr::bind_cols(indicator_list))

  # flip the matrix
  indicator_matrix <- indicator_matrix + t(indicator_matrix) - diag(1, n)

  # replace the sample covariance matrix
  sam_cov <- replace(sam_cov, which(indicator_matrix == 0), 0)

  return(sam_cov)
}


################################################################################

#' Tapering Estimator
#'
#' @description \code{taperingEst} estimates the covariance matrix of a
#'  \code{data.frame}-like object with ordered variables by gradually shrinking
#'  the bands of the sample covariance matrix towards zero. The estimator is
#'  defined as the hadamard product of the sample covariance matrix and a weight
#'  matrix. The amount of shrinkage is dictated by the weight matrix, and is
#'  controlled by a hyperparameter, \code{k}. This estimator is attributed to
#'  \insertRef{cai2010}{cvCovEst}.
#'
#'  The weight matrix is a Toeplitz matrix whose entries are defined as follows:
#'  Let i and j index the rows and columns of the weight matrix, respectively.
#'  If \code{abs(i-j) <= k/2}, then entry i,j in the weight matrix is equal to
#'  1. If \code{k/2 < abs(i-j) < k}, then entry i,j is equal to
#'  \code{2 - 2*abs(i-j)/k}. Otherwise, entry i,j is equal to 0.
#'
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#'
#' @param k A non-negative, even \code{numeric} integer.
#'
#' @importFrom coop covar
#' @importFrom dplyr bind_cols
#'
#' @return A \code{matrix} corresponding to the estimate of the covariance
#'  matrix.
#'
#' @export
taperingEst <- function(dat, k) {

  # compute the sample covariance matrix
  sam_cov <- coop::covar(dat)

  n <- ncol(sam_cov)

  k_h <- k/2

  # loop over different indicies to create weight vectors
  weight_list <- lapply(1:n, function(i) {
    # only consider the lower triangular matrix entries
    j <- i:n

    # calculate the difference in indicies
    di <- abs(i - j)

    # loop over elements in the difference vector and assign weights
    w <- sapply(di, function(d) {

      if (d <= k_h) {
        wi <- 1
      } else if (d > k_h & d < k) {
        wi <- 2 - (d/k_h)
      } else {
        wi <- 0
      }
      return(wi)
    })

    # multiply by corresponding entries in sample covariance matrix
    sam_vec <- sam_cov[j, i]
    sam_vec <- sam_vec * w

    # create a new vector corresponding to lower triangular matrix column
    sam_vec <- c(rep(0, i-1), sam_vec)

    return(sam_vec)
    })

  # combine vectors
  weight_matrix <- suppressMessages(dplyr::bind_cols(weight_list))

  # flip the matrix
  weight_matrix <- weight_matrix + t(weight_matrix) - diag(diag(sam_cov))
  weight_matrix <- as.matrix(weight_matrix)

  # return the new weight matrix
  return(weight_matrix)

}


################################################################################

#' Non-Linear Shrinkage Estimator (Description Not Updated)
#'
#' @description \code{nonlinearShrinkLWEst} invokes the analytical formula
#'  presented by Ledoit and Wolfe for applying a nonlinear shrinkage function
#'  to the sample eigenvalues of the covariance matrix.  The shrinkage function
#'  utilizes the limiting spectral density of the sample eigenvalues and its
#'  Hilbert Transform, both of which are representated using nonparametric
#'  kernel estimation.  The result of the Hilbert transform is that smaller
#'  clusters of eigenvalues are pushed up towards larger clusters, thus creating
#'  a localized shrinkage effect rather than a global one.  For more
#'  information, see \insertRef{Ledoit2018}{cvCovEst}.
#'
#'
#' @param dat A numeric \code{data.frame}, \code{matrix}, or similar object.
#'
#' @importFrom coop covar
#'
#' @return A \code{matrix} corresponding to the estimate of the covariance
#'  matrix.
#'
#' @export
nonlinearShrinkLWEst <- function(dat) {

  n <- nrow(dat)
  p <- ncol(dat)

  # Compute the Sample Covariance Matrix
  sam_cov <- coop::covar(dat)

  # Get the sorted eigenvalues and eigenvectors
  sam_eig <- eigen(sam_cov)

  lambda <- sort(sam_eig$values, decreasing = FALSE, index.return = TRUE)

  u <- sam_eig$vectors[,lambda$ix]

  # Analytical Nonlinear Shrinkage Kernal Formula
  i <- max(1, p - n + 1)
  lambda <- lambda$x[i:p]
  r <- length(lambda)
  c <- min(n, p)
  L <- matrix(lambda, nrow = r, ncol = c)

  # LW Equation 4.9
  h <- n^(-1/3)
  H <- h * t(L)
  x <- (L - t(L))/H

  # LW Equation 4.7
  s1 <- (3/4)/sqrt(5)
  s2 <- -(3/10)/pi
  pos_x <- (1 - (x^2)/5)
  pos_x <- replace(pos_x, list = which(pos_x < 0), 0)
  f_tilde = s1 * rowMeans(pos_x/H)

  # LW Equation 4.8
  log_term <- log(abs((sqrt(5) - x)/(sqrt(5) + x)))
  Hftemp <- (s2 * x) + (s1/pi) * (1 - (x^2)/5) * log_term
  sq5 <- which(abs(x) == sqrt(5))
  Hftemp[sq5] <- s2*x[sq5]
  H_tilde <- rowMeans(Hftemp/H)

  # LW Equation 4.3
  s3 <- pi*(p/n)
  s4 <- 1/(h^2)
  if (p <= n) {
    d_tilde <- lambda/((s3 * lambda * f_tilde)^2 +
                         (1 - (p/n) - s3 *lambda*H_tilde)^2)

  }
  else {
    ones <- rep(1, p-n)
    log_term <- log((1 + sqrt(5)*h)/(1 - sqrt(5)*h))
    m <- mean(1/lambda)

    # LW Equation C.8
    Hf_tilde0 <- (1/pi) * ((3/10)*s4 + (s1/h)*(1 - (1/5)*s4) * log_term) * m

    # LW Equation C.5
    d_tilde0 <- 1/(pi*(p - n))/(n*Hf_tilde0)

    # LW Equation C.4
    d_tilde1 <- lambda/((pi^2 * lambda^2)*(f_tilde^2 + H_tilde^2))
    d_tilde <- c(dtilde0*ones, dtilde1)
  }

  # LW Equation 4.4
  sigma_tilde <- u %*% diag(d_tilde) %*% t(u)

  return(sigma_tilde)
}




