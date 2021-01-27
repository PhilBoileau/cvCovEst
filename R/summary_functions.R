# Main Summary Method and Summary Functions for cvCovEst

################################################################################
#' Summary Statistics of Empirical Risk by Estimator Class
#'
#' @description The \code{empRiskByClass} function calculates the following
#'  summary statistics for the empirical risk within each class of estimator
#'  passed to \code{cvCovEst}: minimum, Q1, median, mean, Q3, and maximum.  The
#'  results are output as a \code{tibble}.
#'
#' @param dat The table of empirical risk calculations which is output by
#'  \code{cvCovEst}.
#'
#' @return \code{\link[tibble]{tibble}} with rows corresponding to estimator
#'  classes and columns corresponding to each summary statistic.
#'
#' @importFrom dplyr group_by summarize arrange %>%
#' @importFrom stats median quantile
#' @importFrom rlang .data
#'
#' @keywords internal
empRiskByClass <- function(dat) {
  emp_risk <- dat %>%
    dplyr::group_by(.data$estimator) %>%
    dplyr::summarise(
      min = min(.data$empirical_risk),
      Q1 = quantile(.data$empirical_risk, probs = 0.25, type = 3),
      median = quantile(.data$empirical_risk, probs = 0.5, type = 3),
      Q3 = quantile(.data$empirical_risk, probs = 0.75, type = 3),
      max = max(.data$empirical_risk),
      mean_risk = mean(.data$empirical_risk),
      .groups = "keep") %>%
    dplyr::arrange(.data$mean_risk)

  return(emp_risk)
}

################################################################################
#' Showing Best Estimator Within Each Class of Estimators
#'
#' @description The \code{bestInClass} function finds the best performing
#'  estimator within each class of estimator passed to \code{cvCovEst} and
#'  finds the associated hyper-parameters if applicable.
#'
#' @param dat The table of empirical risks which is output by \code{cvCovEst}.
#' @param worst This facilitates the option to choose the worst performing
#'  estimator in each class.  Default is \code{FALSE}.
#'
#' @return \code{\link[tibble]{tibble}} with rows corresponding to estimator
#' classes and columns for hyperparameter values and empirical risk for the best
#' estimator in that class.
#'
#' @importFrom dplyr group_by summarize arrange first %>%
#' @importFrom rlang .data
#'
#' @keywords internal
bestInClass <- function(dat, worst = FALSE) {

  if (worst) {
    worst_est <- dat %>%
      dplyr::group_by(.data$estimator) %>%
      dplyr::summarise(
        hyperparameter = dplyr::last(.data$hyperparameters),
        empirical_risk = dplyr::last(.data$empirical_risk),
        .groups = "keep") %>%
      dplyr::arrange(.data$empirical_risk)

    return(worst_est)

  }
  else{
    best_est <- dat %>%
      dplyr::group_by(.data$estimator) %>%
      dplyr::summarise(
        hyperparameter = dplyr::first(.data$hyperparameters),
        empirical_risk = dplyr::first(.data$empirical_risk),
        .groups = "keep") %>%
      dplyr::arrange(.data$empirical_risk)

    return(best_est)

  }
}

################################################################################
#' Summarize Empirical Risks by Class with Hyperparameter
#'
#' @description The \code{hyperRisk} function groups together estimators of the
#'  same class and parses the hyperparameter values over quantiles of the risk.
#'
#' @param dat The table of empirical risk calculations which is output by
#'  \code{cvCovEst}.
#'
#' @return A named \code{list} of data frames. Each list element corresponds to
#'  a \code{\link[tibble]{tibble}} of summary statistics for a specific
#'  estimator class. If no estimators have hyper-parameters, a message is
#'  returned.
#'
#' @importFrom dplyr filter mutate first %>%
#' @importFrom stats quantile
#' @importFrom tibble tibble
#'
#' @keywords internal
hyperRisk <- function(dat) {
  # These are the estimators with hyperparameters
  has_hypers <- c(
    "linearShrinkEst", "thresholdingEst",
    "bandingEst", "taperingEst",
    "scadEst", "poetEst", "robustPoetEst",
    "adaptiveLassoEst"
  )

  estimators <- unique(dat$estimator)

  if (any(has_hypers %in% estimators)) {

    hyper_est <- estimators[which(estimators %in% has_hypers)]

    hyper_summ <- lapply(hyper_est, function(est) {

      h <- dat %>%
        dplyr::filter(.data$estimator == est) %>%
        dplyr::mutate(empirical_risk = round(.data$empirical_risk))

      risk_stats <- quantile(
        h$empirical_risk,
        probs = c(0, 0.25, 0.50, 0.75, 1),
        type = 3
      )

      hyper_risk <- sapply(unname(risk_stats), function(r) {
        # Filter by the quantiles of the empirical risk
        hr <- h %>%
          dplyr::filter(.data$empirical_risk == r)

        vec <- c(
          dplyr::first(hr$hyperparameters),
          dplyr::first(hr$empirical_risk)
        )

        return(vec)
      })

      df <- tibble::tibble(
        t(hyper_risk),
        row.names = c("min", "Q1", "median", "Q3", "max"))

      colnames(df) <- c("hyperparameters", "empirical_risk")

      return(df)
    })

    # Named list of data.frames corresponding to each estimator class
    names(hyper_summ) <- hyper_est
  }
  else{
    hyper_summ <- NULL
    message("No candidate estimators have hyperparameters. hyperRisk = NULL")
  }
  return(hyper_summ)
}

################################################################################
#' Generic Summary Method for cvCovEst
#'
#' @description The \code{summary} method provides summary statistics regarding
#'  the performance of \code{cvCovEst} and can be used for diagnostic plotting.
#'
#' @param object A named \code{list} of class \code{"cvCovEst"}.
#' @param summ_fun A character vector specifying which summaries to output.
#' @param ... Additional arguments passed to summary method.
#'
#' @return A named \code{list} where each element corresponds to the output of
#' of the requested summaries.
#'
#' @importFrom rlang exec
#'
#' @export
summary.cvCovEst <- function(
  object,
  summ_fun = c(
    "empRiskByClass",
    "bestInClass",
    "worstInClass",
    "hyperRisk"),
  ...
) {

  summary_functions <- c(
    "empRiskByClass", "bestInClass", "worstInClass", "hyperRisk")

  # Check cvCovEst credentials
  checkPlotSumArgs(
    object,
    which_fun = "summary",
    summ_fun = summ_fun)

  risk_dat <- object$risk_df

  sums_to_exec <- summary_functions[which(
    summary_functions %in% summ_fun)]

  out <- lapply(sums_to_exec, function(sum_fun) {
    if (sum_fun == "worstInClass") {
      f <- rlang::exec("bestInClass", risk_dat, worst = TRUE)
    }
    else{
      f <- rlang::exec(sum_fun, risk_dat)
    }
    return(f)
  })

  names(out) <- sums_to_exec
  return(out)
}