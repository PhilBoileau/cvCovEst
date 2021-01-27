#' Check Arguments Passed to plot.cvCovEst and summary.cvCovEst
#'
#' @description The \code{checkPlotSumArgs} function verifies that all arguments
#' passed to the \code{plot.cvCovEst} and \code{summary.cvCovEst} functions meet
#' their specifications.  Some additional arguments may be checked at the
#' individual function level.
#'
#' @param dat An object of class, \code{"cvCovEst"}.  Specifically, this is the
#'  standard output of the function \code{cvCovEst}.
#' @param dat_orig The numeric \code{data.frame}, \code{matrix}, or similar
#'  object originally passed to \code{cvCovEst}.
#' @param which_fun A choice of \code{"plot"} or \code{"summary"} depending on.
#'  which function is being checked.
#' @param estimator A character vector specifying one or more classes of
#'  estimators to compare.
#' @param plot_type A character vector specifying one of four choices of
#'  diagnostic plots.
#' @param summ_fun A character vector specifying which summaries to output.
#' @param stat A character vector of one or more summary statistics to use when
#'  comparing estimators.
#' @param k The number of leading/trailing eigenvalues to plot.
#' @param leading A \code{logical} indicating if the leading eigenvalues should
#'  be used.
#' @param abs_v A \code{logical} determining if the absolute value of the matrix
#'  entries should be used for plotting the matrix heatmap.
#'
#' @importFrom assertthat assert_that
#'
#' @return Whether all argument conditions are satisfied.
#'
#' @keywords internal
checkPlotSumArgs <- function(
  dat,
  dat_orig,
  which_fun,
  estimator,
  plot_type,
  summ_fun,
  stat,
  k,
  leading,
  abs_v) {

  # Define possible valid arguments for both functions
  cv_names <- c("estimate", "estimator", "risk_df", "cv_df", "args")

  # Check that names of dat match what is expected of cvCovEst output
  if (is.cvCovEst(dat)) {
    assertthat::assert_that(
      all(cv_names %in% names(dat)) == TRUE,
      msg = paste(
        "cvCovEst object is missing data.",
        "Make sure that the cvCovEst object possess the following elements:",
        "'estimate', 'estimator', 'risk_df', 'cv_df', and 'args'."
      )
    )
  }
  # For plot functions only:
  if (which_fun == "plot") {
    # Define valid plot arguments
    plot_types <- c("summary", "risk", "eigen", "heatmap")
    stat_choices <- c("min", "Q1", "median", "Q3", "max")
    cv_estimators <- unique(dat$risk_df$estimator)

    # Check valid summary stat choices
    assertthat::assert_that(
      all(stat %in% stat_choices) == TRUE,
      msg = "Non-supported summary statistic provided.")
    # Check valid estimators - estimators must have been passed through cvCovEst
    assertthat::assert_that(
      all(estimator %in% cv_estimators) == TRUE,
      msg = "Can only use estimators passed to the cvCovEst function.")
    # Check valid k values
    if (!is.null(k)) {
      assertthat::assert_that(
        k <= ncol(dat_orig),
        msg = "k cannot exceed the number of columns in dat_orig."
      )
    }
    # Check valid plot types
    assertthat::assert_that(
      all(
        plot_type %in% plot_types,
        length(plot_type) == 1) == TRUE,
      msg = "Must provide a single valid plot type.")
    # Check valid logicals
    assertthat::assert_that(
      all(
        is.logical(leading),
        is.logical(abs_v)
      ) == TRUE,
      msg = "Invalid value for logical argument."
    )
    # Check valid summary plot conditions
    if (plot_type == "summary" & !is.null(estimator)) {
      message("Summary plot defaults to the optimal selected estimator.")
    }
  }
  # For summary functions only (which_fun == 'summary')
  else{
    # Define valid summary functions
    summary_functions <- c(
      "empRiskByClass", "bestInClass", "worstInClass", "hyperRisk")

    # Check valid summary functions
    assertthat::assert_that(
      all(summ_fun %in% summary_functions) == TRUE,
      msg = "Must provide a valid summary function."
    )
  }
}
