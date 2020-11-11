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

################################################################################
#' Summarize Risk by Class with Hyperparameter
#'
#' @description The \code{hyperRisk} function groups together estimators of the
#'  same class and selects the hyperparameter values over quantiles of the risk.
#'
#' @param dat The table of empirical risk calculations which is output by
#' \code{cvCovEst}.
#'
#' @return A named \code{list} of data frames.  Each list element corresponds to
#'  a \code{data.frame} of summary statistics for a specific estimator class.
#'  If no estimators have hyper-parameters, a message is returned.
#'
#' @importFrom dplyr filter mutate first %>%
#' @importFrom stats quantile
#'
#' @keywords internal
hyperRisk <- function(dat) {
  # These are the estimators with hyperparameters
  has_hypers <- c("linearShrinkEst", "thresholdingEst",
                  "bandingEst", "taperingEst",
                  "poetEst", "adaptiveLassoEst")

  estimators <- unique(dat$risk_df$estimators)

  if (any(has_hypers %in% estimators)) {

    hyper_est <- estimators[which(estimators %in% has_hypers)]

    hyperSumm <- lapply(hyper_est, function(est) {

      h <- dat$risk_df %>%
        dplyr::filter(estimator == est) %>%
        mutate(empirical_risk = round(empirical_risk))

      risk_stats <- quantile(h$empirical_risk,
                             probs = c(0, 0.25, 0.50, 0.75, 1),
                             type = 3)

      hyper_risk <- sapply(unname(risk_stats), function(r) {
        # Filter by the quantiles of the empirical risk
        hr <- h %>%
          dplyr::filter(empirical_risk == r)

        vec <- c(dplyr::first(hr$hyperparameters),
                 dplyr::first(hr$empirical_risk))

        return(vec)
      })

      df <- data.frame(t(hyper_risk), row.names = names(risk_stats))
      colnames(df) <- c("hyperparameters", "empirical_risk")

      return(df)
    })
    # Named list of data.frames corresponding to each estimator class
    names(hyperSumm) <- hyper_est
  }
  else{

    hyperSumm <- NULL
    message("No estimators have hyperparameters. hyperRisk = NULL")

  }

  return(hyperSumm)

}

