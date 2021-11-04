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
      Min = min(.data$cv_risk),
      Q1 = quantile(.data$cv_risk, probs = 0.25, type = 3),
      Median = quantile(.data$cv_risk, probs = 0.5, type = 3),
      Q3 = quantile(.data$cv_risk, probs = 0.75, type = 3),
      Max = max(.data$cv_risk),
      Mean = mean(.data$cv_risk),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$Mean)

  cv_risk <- dplyr::rename(cv_risk, "Estimator" = "estimator")

  return(cv_risk)
}

################################################################################
#' Showing Best Estimator Within Each Class of Estimators
#'
#' @description \code{bestInClass()} finds the best performing estimator within
#'  each class of estimator passed to \code{\link{cvCovEst}()} and
#'  finds the associated hyperparameters if applicable.
#'
#' @param dat The \code{\link[tibble]{tibble}} of cross-validated risks which is
#'  output by \code{\link{cvCovEst}()}.
#' @param worst This facilitates the option to choose the worst performing
#'  estimator in each class.  Default is \code{FALSE}.
#'
#' @return \code{\link[tibble]{tibble}} with rows corresponding to estimator
#'  classes and columns for hyperparameter values, cross-validated risk, and
#'  other summary metrics for the best (or worst) estimator in that class.
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
        cv_risk = dplyr::last(.data$cv_risk),
        cond_num = dplyr::last(.data$cond_num),
        sign = dplyr::last(.data$sign),
        sparsity = dplyr::last(.data$sparsity),
        .groups = "drop"
      ) %>%
      dplyr::arrange(.data$cv_risk)
  }
  else {
    inClass <- dat %>%
      dplyr::group_by(.data$estimator) %>%
      dplyr::summarise(
        hyperparameter = dplyr::first(.data$hyperparameters),
        cv_risk = dplyr::first(.data$cv_risk),
        cond_num = dplyr::first(.data$cond_num),
        sign = dplyr::first(.data$sign),
        sparsity = dplyr::first(.data$sparsity),
        .groups = "drop"
      ) %>%
      dplyr::arrange(.data$cv_risk)
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
        dplyr::mutate(cv_risk = round(.data$cv_risk))

      risk_stats <- quantile(
        h$cv_risk,
        probs = c(0, 0.25, 0.50, 0.75, 1),
        type = 3
      )

      hyper_risk <- sapply(unname(risk_stats), function(r) {
        # Filter by the quantiles of the empirical risk
        hr <- h %>%
          dplyr::filter(.data$cv_risk == r)

        vec <- c(
          dplyr::first(hr$hyperparameters),
          dplyr::first(hr$cv_risk)
        )

        return(vec)
      })

      df <- data.frame(
        t(hyper_risk),
        stat = c("min", "Q1", "median", "Q3", "max")
      )

      colnames(df) <- c("hyperparameters", "cv_risk", "stat")
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
#' General Matrix Metrics
#'
#' @description \code{matrixMetrics} computes the condition number, sparsity,
#'  and sign of a covariance matrix estimate.
#'
#' @param estimator A \code{matrix} corresponding to a single covariance matrix
#'  estimator.
#'
#' @return A named \code{list} containing the three values.
#'
#' @keywords internal
matrixMetrics <- function(estimate) {
  # Compute the Eigenvalues
  e_vals <- eigen(estimate, symmetric = TRUE, only.values = TRUE)$values
  n <- length(e_vals)

  # Calculate Condition Number
  if (e_vals[n] != 0){
    cn <- round(e_vals[1]/e_vals[n], digits = 3)
  }
  else{
    cn <- 0
  }

  # Determine Matrix "Sign" (positive-definite, positive-semi-definite, etc)
  if (all(e_vals > 0)) {
    sign <- "PD"
  }
  if (any(e_vals == 0) & any(e_vals > 0) & !any(e_vals < 0)) {
    sign <- "PSD"
  }
  if (all(e_vals < 0)) {
    sign <- "ND"
  }
  if (any(e_vals == 0) & any(e_vals < 0) & !any(e_vals > 0)) {
    sign <- "NSD"
  }
  if (any(e_vals > 0) & any(e_vals < 0)) {
    sign <- "IND"
  }

  if (all(estimate == 0)){
    sign <- "NA"
  }

  # Calculate Sparsity Measure
  sp <- round(sum(estimate == 0)/length(estimate), digits = 3)

  metrics <- list(cond_num = cn, sign = sign, sparsity = sp)

  return(metrics)
}


################################################################################
#' Matrix Metrics for cvCovEst Object
#'
#' @description \code{cvMatrixMetrics} computes various metrics and properties
#'  for each covariance matrix estimator candidate's estimate.
#'
#' @param object A named list of class \code{"cvCovEst"} containing the
#'  cross-validated risk assessment.
#'
#' @param dat_orig The \code{numeric data.frame}, \code{matrix}, or similar
#'  object originally passed to \code{\link{cvCovEst}()}.
#'
#' @return A named list of class \code{"cvCovEst"} whose cross-validated risk
#'  assessment is now a \code{\link[tibble]{tibble}} containing the
#'  corresponding metrics for each estimate.  The \code{\link[tibble]{tibble}}
#'  is grouped by estimator and ordered by the primary hyperparameter if
#'  applicable.
#'
#' @importFrom rlang exec .data
#' @importFrom dplyr bind_rows bind_cols group_by arrange %>% select
#'
#' @keywords internal
cvMatrixMetrics <- function(object, dat_orig) {

  mat_mets <- lapply(1:nrow(object$risk_df), function(e){
    # Subset by individual estimator
    est_dat <- object$risk_df[e, ]
    estimator <- as.character(est_dat[1, "estimator"])
    est_args <- list(dat = dat_orig)
    est_attr <- estAttributes(estimator)
    # Get Hypers if Applicable
    if (est_attr[[estimator]][["has_hypers"]]) {
      est_hypers <- getHypers(est_dat, summ_stat = NULL)
      est_args <- append(est_args, est_hypers$hyper_values)
      names(est_args) <- c("dat", est_hypers$hyper_names)
      h1 <- as.numeric(est_hypers$hyper_values[1])
      if (est_attr[[estimator]][["n_hypers"]] == 2) {
        h2 <- as.numeric(est_hypers$hyper_values[2])
      }
      else{
        h2 <- NA
      }
    }
    else{
      h1 <- NA
      h2 <- NA
    }
    # Compute estimate
    est <- rlang::exec(estimator, !!!est_args)
    # Compute metrics
    est_metrics <- matrixMetrics(est)
    met_names <- names(est_metrics)
    est_metrics <- append(est_metrics, c(h1, h2))
    names(est_metrics) <- c(met_names, "hyper1", "hyper2")

    return(est_metrics)
  })

  mat_mets <- dplyr::bind_rows(mat_mets) %>%
    dplyr::bind_cols(object$risk_df) %>%
    dplyr::group_by(.data$estimator) %>%
    dplyr::arrange(.data$hyper1, by_group = TRUE) %>%
    dplyr::select(.data$estimator, .data$hyperparameters, .data$cv_risk,
                  .data$cond_num, .data$sign, .data$sparsity, .data$hyper1,
                  .data$hyper2)
  object$risk_df <- mat_mets

  return(object)
}

################################################################################
#' Generic Summary Method for cvCovEst
#'
#' @description \code{summary()} provides summary statistics regarding
#'  the performance of \code{\link{cvCovEst}()} and can be used for diagnostic
#'  plotting.
#'
#' @param object A named \code{list} of class \code{"cvCovEst"}.
#' @param dat_orig The \code{numeric data.frame}, \code{matrix}, or similar
#'  object originally passed to \code{\link{cvCovEst}()}.
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
#'       applicable, of the best performing estimator within each class along
#'       with other metrics.
#'     \item \code{"worstInClass"} - Returns the specific hyperparameters, if
#'       applicable, of the worst performing estimator within each class along
#'       with other metrics.
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
#' summary(cv_dat, mtcars)
#' @export
summary.cvCovEst <- function(
                             object,
                             dat_orig,
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

  object <- cvMatrixMetrics(object, dat_orig)
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

  if (length(out) == 1)
    out <- out[[1]]
  else
    names(out) <- sums_to_exec

  return(out)
}
