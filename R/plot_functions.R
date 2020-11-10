################################################################################
# Plotting Functions for cvCovEst


################################################################################
#' Summary Statistics of Empirical Risk by Estimator Class
#'
#' @description The \code{empRiskByClass} function calculates the following
#' summary statistics for the empirical risk within each class of estimator
#' passed to \code{cvCovEst}: minimum, Q1, median, mean, Q3, and maximum.  The
#' results are output as a \code{tibble}.
#'
#' @param dat The table of empirical risk calculations which is output by
#' \code{cvCovEst}.
#'
#' @return A \code{data.frame} with rows corresponding to estimator classes and
#' columns corresponding to each summary statistic.
#'
#' @importFrom dplyr group_by summarize arrange %>%
#' @importFrom stats median quantile
#'
#' @keywords internal
empRiskByClass <- function(dat) {
  empRisk <- dat %>%
    dplyr::group_by(estimator) %>%
    dplyr::summarise(min_risk = min(empirical_risk),
                     Q1_risk = quantile(empirical_risk, probs = 0.25),
                     median_risk = quantile(empirical_risk, probs = 0.5),
                     mean_risk = mean(empirical_risk),
                     Q3_risk = quantile(empirical_risk, probs = 0.75),
                     max_risk = max(empirical_risk),
                     .groups = "keep") %>%
    dplyr::arrange(mean_risk)

  return(empRisk)
}

################################################################################
#' Showing Best Estimator Within Each Class of Estimators
#'
#' @description The \code{bestInClass} function finds the best performing
#'  estimator within each class of estimator passed to \code{cvCovEst} and
#'  finds the associated hyperparameters if applicable.
#'
#' @param dat The table of empirical risk calculations which is output by
#' \code{cvCovEst}.
#'
#' @return A \code{data.frame} with rows corresponding to estimator classes and
#'  columns for hyperparameter values and empirical risk for the best estimator
#'  in that class.
#'
#' @importFrom dplyr group_by summarize arrange first %>%
#'
#' @keywords internal
bestInClass <- function(dat) {
  bestEst <- dat %>%
    dplyr::group_by(estimator) %>%
    dplyr::summarise(hyperparameter = first(hyperparameters),
                     empirical_risk = first(empirical_risk),
                     .groups = "keep") %>%
    dplyr::arrange(empirical_risk)

  return(bestEst)
}

