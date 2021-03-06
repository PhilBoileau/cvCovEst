% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{cvEigenPlot}
\alias{cvEigenPlot}
\title{Eigenvalue Plot}
\usage{
cvEigenPlot(
  dat,
  estimator,
  stat = "min",
  dat_orig,
  k,
  leading = TRUE,
  plot_type = "eigen",
  cv_details,
  has_hypers
)
}
\arguments{
\item{dat}{A named \code{list}.  Specifically, this is the standard output of
\code{\link{cvCovEst}()}.}

\item{estimator}{A \code{character} vector specifying one or more classes of
estimators to compare.}

\item{stat}{A \code{character} vector containing the names of various
cross-validated risk summary statistics.  Within each class of estimator,
eigenvalues will be plot for the estimators corresponding to each stat.}

\item{dat_orig}{The \code{numeric data.frame}, \code{matrix}, or similar
object originally passed to \code{\link{cvCovEst}()}.}

\item{k}{A \code{numeric} indicating the number of eigenvalues to plot. Must
be less than or equal to the number of columns of the original data matrix.}

\item{leading}{A \code{logical} indicating if the leading eigenvalues should
be used.  Default is \code{TRUE}.  If \code{FALSE}, the trailing
eigenvalues will be used instead.}

\item{plot_type}{A \code{character} detailing the type of plot. Passed to
\code{theme_cvCovEst}, defaults to \code{"risk"}}

\item{cv_details}{A \code{character} vector summarizing key arguments passed
to \code{\link{cvCovEst}()}.}

\item{has_hypers}{A \code{character} vector containing the names of current
estimators with hyperparameters.}
}
\value{
A plot, or grid of plots, showing the \code{k} leading or trailing
eigenvalues of the specified estimators and associated summary statistics of
the cross-validated risk.
}
\description{
\code{cvEigenPlot()} plots the eigenvalues of one or more
estimators produced by \code{\link{cvCovEst}()}.
}
\keyword{internal}
