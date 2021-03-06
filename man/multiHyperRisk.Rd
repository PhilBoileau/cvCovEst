% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{multiHyperRisk}
\alias{multiHyperRisk}
\title{Multi-Hyperparameter Risk Plots}
\usage{
multiHyperRisk(dat, estimator, switch_vars = FALSE, min_max = FALSE)
}
\arguments{
\item{dat}{A \code{data.frame} of cross-validated risks. Specifically, this
is the \code{risk_df} table output by \code{\link{cvCovEst}()}.}

\item{estimator}{A \code{character} vector specifying one or more classes of
estimators to compare.}

\item{switch_vars}{A \code{logical} indicating if the x-axis and factor
variables should be switched.  Default is \code{FALSE}.}

\item{min_max}{A \code{logical}. If \code{TRUE}, only the minimum and
maximum values of the factor hyperparameter will be used. Defaults to
\code{FALSE}.}
}
\value{
A named \code{list} of plots.
}
\description{
\code{multiHyperRisk()} produces plots of the cross-validated
risk for estimators with more than one hyperparameter.  The function
transforms one of the hyperparameters into a factor and uses it to
distinguish between the risk of various estimators.  If one of the
hyperparameters has only one unique value, that hyperparameter is used as
the factor variable.  If all hyperparameters have only one unique value, a
plot is not generated for that estimator class.
}
\keyword{internal}
