# Main Summary Method and Summary Functions for cvCovEst

################################################################################
#' Summary Statistics of Cross-Validated Risk by Estimator Class
#'
#' @description \code{cvRiskByClass()} calculates the following
#'  summary statistics for the cross-validated risk within each class of
#'  estimator passed to \code{\link{cvCovEst}()}: minimum, Q1, median, mean, Q3,
#'  and maximum.  The results are output as a \code{tibble}.
#'
#' @param dat The \code{\link[tibble]{tibble}} of cross-validated risk
#'  calculations which is output by \code{\link{cvCovEst}()}.
#'
#' @return \code{\link[tibble]{tibble}} with rows corresponding to estimator
#'  classes and columns corresponding to each summary statistic.
#'
#' @importFrom dplyr group_by summarize arrange %>% rename
#' @importFrom stats median quantile
#' @importFrom rlang .data
#'
#' @keywords internal
cvRiskByClass <- function(dat) {
  cv_risk <- dat %>%
    dplyr::group_by(.data$estimator) %>%
    dplyr::summarise(
      Min = min(.data$empirical_risk),
      Q1 = quantile(.data$empirical_risk, probs = 0.25, type = 3),
      Median = quantile(.data$empirical_risk, probs = 0.5, type = 3),
      Q3 = quantile(.data$empirical_risk, probs = 0.75, type = 3),
      Max = max(.data$empirical_risk),
      Mean = mean(.data$empirical_risk),
      .groups = "keep"
    ) %>%
    dplyr::arrange(.data$Mean)

  cv_risk <- dplyr::rename(cv_risk, Estimator = estimator)

  return(cv_risk)
}

################################################################################
#' Showing Best Estimator Within Each Class of Estimators
#'
#' @description \code{bestInClass()} finds the best performing estimator within
#'  each class of estimator passed to \code{\link{cvCovEst}()} and
#'  finds the associated hyper-parameters if applicable.
#'
#' @param dat The \code{\link[tibble]{tibble}} of cross-validated risks which is
#'  output by \code{\link{cvCovEst}()}.
#' @param worst This facilitates the option to choose the worst performing
#'  estimator in each class.  Default is \code{FALSE}.
#'
#' @return \code{\link[tibble]{tibble}} with rows corresponding to estimator
#'  classes and columns for hyperparameter values and cross-validated risk for
#'  the best estimator in that class.
#'
#' @importFrom dplyr group_by summarize arrange first %>%
#' @importFrom rlang .data
#'
#' @keywords internal
bestInClass <- function(dat, worst = FALSE) {
  if (worst) {
    inClass <- dat %>%
      dplyr::group_by(.data$estimator) %>%
      dplyr::summarise(
        hyperparameter = dplyr::last(.data$hyperparameters),
        empirical_risk = dplyr::last(.data$empirical_risk),
        .groups = "keep"
      ) %>%
      dplyr::arrange(.data$empirical_risk)
  }
  else {
    inClass <- dat %>%
      dplyr::group_by(.data$estimator) %>%
      dplyr::summarise(
        hyperparameter = dplyr::first(.data$hyperparameters),
        empirical_risk = dplyr::first(.data$empirical_risk),
        .groups = "keep"
      ) %>%
      dplyr::arrange(.data$empirical_risk)
  }
  return(inClass)
}

################################################################################
#' Summarize Cross-Validated Risks by Class with Hyperparameter
#'
#' @description \code{hyperRisk()} groups together estimators of the
#'  same class and parses the hyperparameter values over quantiles of the risk.
#'
#' @param dat The \code{\link[tibble]{tibble}} of cross-validated risk calculations
#'  which is output by \code{\link{cvCovEst}()}.
#'
#' @return A named \code{list} of data frames. Each list element corresponds to
#'  a \code{\link[tibble]{tibble}} of summary statistics for a specific
#'  estimator class. If no estimators have hyper-parameters, a message is
#'  returned.
#'
#' @importFrom dplyr filter mutate first %>%
#' @importFrom stats quantile
#' @importFrom tibble as_tibble
#'
#' @keywords internal
hyperRisk <- function(dat) {
  estimators <- unique(dat$estimator)

  # Get Attributes
  attr_df <- estAttributes(estimator = estimators)
  attr_df <- dplyr::bind_rows(attr_df)
  attr_df$estimator <- estimators
  has_hypers <- attr_df$estimator[which(attr_df$has_hypers)]

  if (any(estimators %in% has_hypers)) {
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

      df <- data.frame(
        t(hyper_risk),
        stat = c("min", "Q1", "median", "Q3", "max")
      )

      colnames(df) <- c("hyperparameters", "empirical_risk", "stat")
      df <- tibble::as_tibble(df)

      return(df)
    })

    # Named list of data.frames corresponding to each estimator class
    names(hyper_summ) <- hyper_est
  }
  else {
    hyper_summ <- NULL
    message("No candidate estimators have hyperparameters. hyperRisk = NULL")
  }
  return(hyper_summ)
}

################################################################################
#' Generic Summary Method for cvCovEst
#'
#' @description \code{summary()} provides summary statistics regarding
#'  the performance of \code{\link{cvCovEst}()} and can be used for diagnostic
#'  plotting.
#'
#' @param object A named \code{list} of class \code{"cvCovEst"}.
#' @param summ_fun A \code{character} vector specifying which summaries to
#'  output.  See Details for function descriptions.
#' @param ... Additional arguments passed to \code{summary()}These are not
#'  explicitly used and should be ignored by the user.
#'
#' @details \code{summary()} accepts four different choices for the
#'  \code{summ_fun} argument.  The choices are:
#'  \itemize{
#'     \item \code{"cvRiskByClass"} - Returns the minimum, first quartile,
#'       median, third quartile, and maximum of the cross-validated risk
#'       associated with each class of estimator passed to
#'       \code{\link{cvCovEst}()}.
#'     \item \code{"bestInClass"} - Returns the specific hyperparameters, if
#'       applicable, of the best performing estimator within each class.
#'     \item \code{"worstInClass"} - Returns the specific hyperparameters, if
#'       applicable, of the worst performing estimator within each class.
#'     \item \code{"hyperRisk"} - For estimators that take hyperparameters as
#'       arguments, this function returns the hyperparameters associated with
#'       the minimum, first quartile, median, third quartile, and maximum of the
#'       cross-validated risk within each class of estimator. Each class has
#'       its own \code{\link[tibble]{tibble}}, which are returned as a
#'       \code{list}.
#'  }
#'
#' @return A named \code{list} where each element corresponds to the output of
#' of the requested summaries.
#'
#' @importFrom rlang exec
#'
#' @examples
#' cv_dat <- cvCovEst(
#'   dat = mtcars,
#'   estimators = c(
#'     linearShrinkEst, thresholdingEst, sampleCovEst
#'   ),
#'   estimator_params = list(
#'     linearShrinkEst = list(alpha = seq(0.1, 0.9, 0.1)),
#'     thresholdingEst = list(gamma = seq(0.1, 0.9, 0.1))
#'   ),
#'   center = TRUE,
#'   scale = TRUE
#' )
#'
#' summary(cv_dat)
#' @export
summary.cvCovEst <- function(
                             object,
                             summ_fun = c(
                               "cvRiskByClass",
                               "bestInClass",
                               "worstInClass",
                               "hyperRisk"
                             ),
                             ...) {
  summary_functions <- c(
    "cvRiskByClass", "bestInClass", "worstInClass", "hyperRisk"
  )

  # Check cvCovEst credentials
  checkPlotSumArgs(
    object,
    which_fun = "summary",
    summ_fun = summ_fun
  )

  risk_dat <- object$risk_df

  sums_to_exec <- summary_functions[which(
    summary_functions %in% summ_fun
  )]

  out <- lapply(sums_to_exec, function(sum_fun) {
    if (sum_fun == "worstInClass") {
      f <- rlang::exec("bestInClass", risk_dat, worst = TRUE)
    }
    else {
      f <- rlang::exec(sum_fun, risk_dat)
    }
    return(f)
  })

  names(out) <- sums_to_exec
  return(out)
}
