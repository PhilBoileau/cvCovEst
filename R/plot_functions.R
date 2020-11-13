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
    dplyr::summarise(
      min_risk = min(empirical_risk),
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
#' @param worst This facilitates the option to choose the worst performing
#' estimator in each class.  The default is FALSE.
#'
#' @return A \code{data.frame} with rows corresponding to estimator classes and
#'  columns for hyperparameter values and empirical risk for the best estimator
#'  in that class.
#'
#' @importFrom dplyr group_by summarize arrange first %>%
#'
#' @keywords internal
bestInClass <- function(dat, worst = FALSE) {

  if (worst) {
    worstEst <- dat %>%
      dplyr::group_by(estimator) %>%
      dplyr::summarise(
        hyperparameter = dplyr::last(hyperparameters),
        empirical_risk = dplyr::last(empirical_risk),
        .groups = "keep") %>%
      dplyr::arrange(empirical_risk)

    return(worstEst)

  }
  else{
    bestEst <- dat %>%
      dplyr::group_by(estimator) %>%
      dplyr::summarise(
        hyperparameter = dplyr::first(hyperparameters),
        empirical_risk = dplyr::first(empirical_risk),
        .groups = "keep") %>%
      dplyr::arrange(empirical_risk)

    return(bestEst)

  }
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
  has_hypers <- c(
    "linearShrinkEst", "thresholdingEst",
    "bandingEst", "taperingEst",
    "poetEst", "adaptiveLassoEst"
    )

  estimators <- unique(dat$estimator)

  if (any(has_hypers %in% estimators)) {

    hyper_est <- estimators[which(estimators %in% has_hypers)]

    hyperSumm <- lapply(hyper_est, function(est) {

      h <- dat %>%
        dplyr::filter(
          estimator == est
          ) %>%
        dplyr::mutate(
          empirical_risk = round(empirical_risk)
          )

      risk_stats <- quantile(
        h$empirical_risk,
        probs = c(0, 0.25, 0.50, 0.75, 1),
        type = 3
        )

      hyper_risk <- sapply(unname(risk_stats), function(r) {
        # Filter by the quantiles of the empirical risk
        hr <- h %>%
          dplyr::filter(
            empirical_risk == r
            )

        vec <- c(
          dplyr::first(hr$hyperparameters),
          dplyr::first(hr$empirical_risk)
          )

        return(vec)
      })

      df <- data.frame(
        t(hyper_risk),
        row.names = names(risk_stats)
        )

      colnames(df) <- c(
        "hyperparameters",
        "empirical_risk"
        )

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

################################################################################
#' Summary Function for cvCovEst
#'
#' @description The \code{cvSummary} provides summary statistics regarding the
#' performance of \code{cvCovEst} and can be used for diagnostic plotting.
#'
#' @param dat A named \code{list}.  Specifically, this is the standard output of
#' \code{cvCovEst}.
#'
#' @param summary A character vector specifying which summaries to output.  The
#' default is \code{'all'}.
#'
#' @return A named \code{list} where each element corresponds to the output of
#' of the requested summaries.
#'
#' @importFrom rlang exec
#'
#' @keywords external
cvSummary <- function(dat, summary = 'all') {

  risk_dat <- dat$risk_df

  summary_functions <- c(
    "empRiskByClass",
    "bestInClass",
    "worstInClass",
    "hyperRisk"
  )

  if (summary == 'all') {
    sums_to_exec <- summary_functions
  }
  else{
    which_sum <- which(
      summary_functions %in% summary
    )

    sums_to_exec <- summary_functions[which_sum]
  }

  out = lapply(
    sums_to_exec,
    function(sum_fun) {
      if (sum_fun == 'worstInClass') {
        f <- rlang::exec(
          'bestInClass',
          risk_dat,
          worst = TRUE
          )
      }
      else{
        f <- rlang::exec(
          sum_fun,
          risk_dat
          )
      }
      return(f)
    }
  )
  names(out) <- sums_to_exec
  return(out)

  }

################################################################################
#' Single Heat Map Plot
#'
#' @description The \code{cvSingleMelt} plots a heat map visualization of a
#' single estimator for the covariance matrix.
#'
#' @param dat A named \code{list}.  Specifically, this is the standard output of
#' \code{cvCovEst}.
#'
#' @param estimator A single string specifying which class of estimator to
#' visualize.
#'
#' @param stat A single specifying the performance of the chosen estimator.  The
#'  two options are \code{"best"} for the best performing estimator and
#'  \code{"worst"} for the worst performing estimator.
#'
#' @param dat_orig The numeric \code{data.frame}, \code{matrix}, or similar
#'  object originally passed to \code{cvCovEst}.
#'
#' @return A heat map plot of the desired covariance matrix estimator.
#'
#' @importFrom rlang exec
#' @importFrom ggplot2 ggplot geom_raster
#'
#' @keywords external
cvSingleMelt <- function(dat, estimator, stat, dat_orig) {
  # add checks here for assuring cvCovEst output

  has_hypers <- c(
    "linearShrinkEst", "thresholdingEst",
    "bandingEst", "taperingEst",
    "poetEst", "adaptiveLassoEst"
  )

  best_vs_worst <- cvSummary(
    dat,
    summary = c(
      'bestInClass',
      'worstInClass'
      )
    )

  if (stat == 'best') {
    est <- best_vs_worst$bestInClass %>%
      dplyr::filter(
        estimator == estimator
      )
  }
  else{
    est <- best_vs_worst$worstInClass %>%
      dplyr::filter(
        estimator == estimator
      )
  }

  if (estimator %in% has_hypers) {
    estHypers <- lapply(
      stringr::str_split(
        est[1, 2], ", "
        ),
      function(s) {
        h <- stringr::str_split(
          s, "= "
          ) %>% unlist()

        return(
          as.numeric(h[2])
          )
        }
      )

    estArgs <- list(
      dat = dat_orig,
      estHypers %>% unlist()
    )

    estimate <- rlang::exec(
      estimator,
      !!!estArgs
      )
  }
  else{
    estimate <- rlang::exec(
      estimator
    )
  }

  meltEst <- reshape2::melt(estimate)
  return(meltEst)
}




