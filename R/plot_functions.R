# Plotting Functions for cvCovEst

################################################################################
#' Multiple Heat Map Plot
#'
#' @description \code{cvMultiMelt()} visualizes the structure of one or more
#'  covariance matrix estimators through a grid of heat maps, where each heat
#'  map corresponds to a different estimator.
#'
#' @param dat A named \code{list}.  Specifically, this is the standard output of
#'  \code{cvCovEst}.
#' @param estimator A \code{character} vector specifying one or more classes of
#'  estimators to visualize and compare.
#' @param stat A \code{character} vector containing the names of various
#'  cross-validated risk summary statistics. Estimators corresponding to each
#'  statistics will be visualized with a different heatmap. Default is
#'  \code{'min'} for minimum cross-validated risk.
#' @param dat_orig The \code{numeric data.frame}, \code{matrix}, or similar
#'  object originally passed to \code{\link{cvCovEst}()}.
#' @param plot_type A \code{character} detailing the type of plot. Passed to
#'  \code{theme_cvCovEst}, defaults to \code{"risk"}
#' @param cv_details A \code{character} vector summarizing key arguments passed
#'  to \code{\link{cvCovEst}()}.
#' @param has_hypers A \code{character} vector containing the names of current
#'  estimators with hyperparameters.
#' @param abs_v A \code{logical} determining if the absolute value of the matrix
#'  entries should be displayed versus the signed value. Default is \code{TRUE}.
#'
#' @return A grid of heat map plots comparing the desired covariance matrix
#' estimators.
#'
#' @importFrom rlang exec .data
#' @importFrom dplyr filter %>% bind_rows mutate_if
#' @importFrom ggplot2 ggplot aes geom_raster facet_wrap labs facet_grid vars scale_fill_gradient scale_fill_gradient2
#' @importFrom assertthat assert_that
#' @importFrom RColorBrewer brewer.pal
#'
#' @keywords internal
cvMultiMelt <- function(
                        dat,
                        estimator,
                        stat = "min",
                        dat_orig,
                        plot_type = "heatmap",
                        cv_details,
                        has_hypers,
                        abs_v = TRUE) {
  stat_choices <- c("min", "Q1", "median", "Q3", "max")

  single_stat <- ifelse(length(stat) == 1, TRUE, FALSE)

  single_est <- ifelse(length(estimator) == 1, TRUE, FALSE)

  # Center and Scale Original Data to Match Call to cvCovEst
  dat_orig <- safeColScale(
    dat_orig,
    center = dat$args$center,
    scale = dat$args$scale
  ) %>%
    unname()

  # Call cvSummary
  cv_sum <- summary.cvCovEst(dat)

  # Setup Values
  blues <- RColorBrewer::brewer.pal(n = 9, name = "Blues")

  legend_title <- ifelse(abs_v, "Absolute Value", "Value")

  plot_title <- ifelse(
    dat$args$scale, "Correlation Matrix Heatmap", "Covariance Matrix Heatmap"
  )

  # Single Stat Option
  if (single_stat) {
    stat_melts <- lapply(estimator, function(est) {
      # For Estimators With Hyper-parameters
      if (est %in% has_hypers) {
        estimator_stats <- cv_sum$hyperRisk[[est]]
        # Get The Associated Hyper-parameters
        est_hypers <- getHypers(dat = estimator_stats, summ_stat = stat)

        # Collect Args & Run The Associated Estimator
        dat <- list(dat_orig)
        arg_names <- c("dat", est_hypers$hyper_names)

        est_args <- append(dat, est_hypers$hyper_values)
        names(est_args) <- arg_names

        estimate <- rlang::exec(est, !!!est_args)
      }
      # If No Hyper-parameters
      else {
        estimate <- rlang::exec(est, dat_orig)
      }

      # Create Melted Data Frame
      if (abs_v) {
        estimate <- abs(estimate)
      }

      melt_est <- as.data.frame.table(
        estimate,
        responseName = "value"
      ) %>%
        dplyr::mutate_if(is.factor, as.integer)

      # Label by Estimator
      est_name <- rep(est, nrow(melt_est))

      melt_est$Var1 <- rev(melt_est$Var1)
      melt_est$estimator <- est_name

      return(melt_est)
    })

    # Combine and Re-factor
    stat_melts <- dplyr::bind_rows(stat_melts)
    stat_melts$estimator <- factor(
      stat_melts$estimator,
      levels = cv_sum$bestInClass$estimator
    )

    plot <- ggplot(stat_melts, aes(x = .data$Var1, y = .data$Var2)) +
      geom_raster(aes(fill = .data$value)) +
      facet_wrap(facets = vars(.data$estimator)) +
      labs(
        title = plot_title, fill = legend_title,
        caption = cv_details, y = "Covariates"
      )

    if (abs_v) {
      plot <- plot +
        scale_fill_gradient(low = "white", high = blues[9]) +
        theme_cvCovEst(plot_type = plot_type)
    }
    else {
      plot <- plot +
        scale_fill_gradient2(low = "gold", mid = "white", high = blues[9]) +
        theme_cvCovEst(plot_type = plot_type)
    }
    return(plot)
  }

  # Multi-Stat Option
  else {
    # Multi-Stat for Single Estimator
    if (single_est) {
      assertthat::assert_that(
        estimator %in% has_hypers,
        msg = "The chosen estimator has no hyper-parameters.
        All plots will be the same."
      )

      estimator_stats <- cv_sum$hyperRisk[[estimator]]
      stat_melts <- lapply(stat, function(stat_est) {
        # Get The Associated Hyper-parameters
        est_hypers <- getHypers(dat = estimator_stats, summ_stat = stat_est)

        # Collect Args & Run The Associated Estimator
        dat <- list(dat_orig)
        arg_names <- c("dat", est_hypers$hyper_names)

        est_args <- append(dat, est_hypers$hyper_values)
        names(est_args) <- arg_names

        estimate <- rlang::exec(estimator, !!!est_args)

        # Create Melted Data Frame
        if (abs_v) {
          estimate <- abs(estimate)
        }

        melt_est <- as.data.frame.table(
          estimate,
          responseName = "value"
        ) %>%
          dplyr::mutate_if(is.factor, as.integer)

        # Label by Stat
        stat_name <- rep(stat_est, nrow(melt_est))

        melt_est$Var1 <- rev(melt_est$Var1)
        melt_est$summary_stat <- stat_name

        return(melt_est)
      })

      # Combine and Re-factor
      stat_melts <- dplyr::bind_rows(stat_melts)
      stat_melts$summary_stat <- factor(
        stat_melts$summary_stat,
        levels = stat_choices
      )

      plot <- ggplot(stat_melts, aes(x = .data$Var1, y = .data$Var2)) +
        geom_raster(aes(fill = .data$value)) +
        facet_wrap(facets = vars(.data$summary_stat), nrow = 1) +
        labs(title = plot_title, fill = legend_title, caption = cv_details)

      if (abs_v) {
        plot <- plot +
          scale_fill_gradient(low = "white", high = blues[9]) +
          theme_cvCovEst(plot_type = plot_type)
      }
      else {
        plot <- plot +
          scale_fill_gradient2(low = "gold", mid = "white", high = blues[9]) +
          theme_cvCovEst(plot_type = plot_type)
      }
      return(plot)
    }
    else {
      # Multi-Stat for Several Estimators
      multi_melts <- lapply(estimator, function(est) {
        estimator_stats <- cv_sum$hyperRisk[[est]]

        stat_melts <- lapply(stat, function(stat_est) {

          # Get The Associated Hyper-parameters
          est_hypers <- getHypers(
            dat = estimator_stats,
            summ_stat = stat_est
          )

          # Run The Associated Estimator
          dat <- list(dat_orig)
          arg_names <- c("dat", est_hypers$hyper_names)

          est_args <- append(dat, est_hypers$hyper_values)
          names(est_args) <- arg_names

          estimate <- rlang::exec(est, !!!est_args)

          # Create Melted Data Frame
          if (abs_v) {
            estimate <- abs(estimate)
          }

          melt_est <- as.data.frame.table(
            estimate,
            responseName = "value"
          ) %>%
            dplyr::mutate_if(is.factor, as.integer)

          # Label by Stat
          stat_name <- rep(stat_est, nrow(melt_est))

          # Label by Estimator
          est_name <- rep(est, nrow(melt_est))

          melt_est$Var1 <- rev(melt_est$Var1)
          melt_est$summary_stat <- stat_name
          melt_est$estimator <- est_name

          return(melt_est)
        })

        stat_melts <- dplyr::bind_rows(stat_melts)

        return(stat_melts)
      })

      # Combine and Re-factor
      multi_melts <- dplyr::bind_rows(multi_melts)

      multi_melts$summary_stat <- factor(
        multi_melts$summary_stat,
        levels = stat_choices
      )
      multi_melts$estimator <- factor(
        multi_melts$estimator,
        levels = cv_sum$bestInClass$estimator
      )

      plot <- ggplot(multi_melts, aes(x = .data$Var1, y = .data$Var2)) +
        geom_raster(aes(fill = .data$value)) +
        facet_grid(
          rows = vars(.data$estimator),
          cols = vars(.data$summary_stat)
        ) +
        labs(title = plot_title, fill = legend_title, caption = cv_details)

      if (abs_v) {
        plot <- plot +
          scale_fill_gradient(low = "white", high = blues[9]) +
          theme_cvCovEst(plot_type = plot_type)
      }
      else {
        plot <- plot +
          scale_fill_gradient2(low = "gold", mid = "white", high = blues[9]) +
          theme_cvCovEst(plot_type = plot_type)
      }
      return(plot)
    }
  }
}

################################################################################
#' Eigenvalue Plot
#'
#' @description \code{cvEigenPlot()} plots the eigenvalues of one or more
#' estimators produced by \code{\link{cvCovEst}()}.
#'
#' @param dat A named \code{list}.  Specifically, this is the standard output of
#'  \code{\link{cvCovEst}()}.
#' @param estimator A \code{character} vector specifying one or more classes of
#'  estimators to compare.
#' @param stat A \code{character} vector containing the names of various
#'  cross-validated risk summary statistics.  Within each class of estimator,
#'  eigenvalues will be plot for the estimators corresponding to each stat.
#' @param dat_orig The \code{numeric data.frame}, \code{matrix}, or similar
#'  object originally passed to \code{\link{cvCovEst}()}.
#' @param k A \code{numeric} indicating the number of eigenvalues to plot. Must
#'  be less than or equal to the number of columns of the original data matrix.
#' @param leading A \code{logical} indicating if the leading eigenvalues should
#'  be used.  Default is \code{TRUE}.  If \code{FALSE}, the trailing
#'  eigenvalues will be used instead.
#' @param plot_type A \code{character} detailing the type of plot. Passed to
#'  \code{theme_cvCovEst}, defaults to \code{"risk"}
#' @param cv_details A \code{character} vector summarizing key arguments passed
#'  to \code{\link{cvCovEst}()}.
#' @param has_hypers A \code{character} vector containing the names of current
#'  estimators with hyperparameters.
#'
#' @return A plot, or grid of plots, showing the \code{k} leading or trailing
#'  eigenvalues of the specified estimators and associated summary statistics of
#'  the cross-validated risk.
#'
#' @importFrom ggplot2 ggplot aes vars geom_point geom_path scale_x_continuous facet_wrap labs scale_color_viridis_d
#' @importFrom RColorBrewer brewer.pal
#' @importFrom RSpectra eigs_sym
#' @importFrom rlang exec .data
#' @importFrom dplyr %>%
#' @importFrom tibble tibble
#'
#' @keywords internal
cvEigenPlot <- function(
                        dat,
                        estimator,
                        stat = "min",
                        dat_orig,
                        k,
                        leading = TRUE,
                        plot_type = "eigen",
                        cv_details,
                        has_hypers) {
  stat_choices <- c("min", "Q1", "median", "Q3", "max")

  single_stat <- ifelse(length(stat) == 1, TRUE, FALSE)

  single_est <- ifelse(length(estimator) == 1, TRUE, FALSE)

  # Center and Scale Original Data to Match Call to cvCovEst
  dat_orig <- safeColScale(
    dat_orig,
    center = dat$args$center,
    scale = dat$args$scale
  ) %>% unname()

  # Determine leading/trailing and set index accordingly
  eig_type <- ifelse(leading, "LA", "SA")
  p <- ncol(dat_orig)
  if (leading) {
    index <- seq(k)
  }
  else {
    index <- seq(p - k + 1, p)
  }

  # Get Summary Output
  cv_sum <- summary.cvCovEst(dat)
  blues <- RColorBrewer::brewer.pal(9, "Blues")
  b <- ifelse(leading, "1", as.character(p))

  # Single Estimator Option
  if (single_est) {
    # Has Hypers
    if (estimator %in% has_hypers) {
      estimator_stats <- cv_sum$hyperRisk[[estimator]]

      stat_eigs <- lapply(stat, function(stat_est) {
        # Get The Associated Hyper-parameters
        est_hypers <- getHypers(dat = estimator_stats, summ_stat = stat_est)

        # Run The Associated Estimator
        dat <- list(dat_orig)
        arg_names <- c("dat", est_hypers$hyper_names)

        est_args <- append(dat, est_hypers$hyper_values)
        names(est_args) <- arg_names

        estimate <- rlang::exec(estimator, !!!est_args)

        # Get the eigenvalues
        est_eigs <- tibble::tibble(
          index = index,
          eigenvalues = suppressWarnings(RSpectra::eigs_sym(
            estimate,
            k = k, which = eig_type
          )$values),
          stat = rep(stat_est, k),
          estimator = rep(estimator, k)
        )
        return(est_eigs)
      })

      # Combine and Re-factor
      stat_eigs <- dplyr::bind_rows(stat_eigs)
    }
    else {
      # No Hypers
      if (!single_stat) {
        message(
          "Estimator has no hyperparameters.  All stats yield same estimate."
        )
      }
      estimate <- rlang::exec(estimator, dat_orig)

      # Create Data Frame
      stat_eigs <- tibble::tibble(
        index = index,
        eigenvalues = suppressWarnings(RSpectra::eigs_sym(
          estimate,
          k = k, which = eig_type
        )$values),
        stat = rep(stat[1], k),
        estimator = rep(estimator, k)
      )
    }
    # Re-factor
    stat_eigs$stat <- factor(stat_eigs$stat, levels = stat_choices)

    # Generate Plot
    if (k == 1) {
      plot <- ggplot(
        stat_eigs,
        aes(x = .data$index, y = .data$eigenvalues, color = .data$stat)
      ) +
        geom_point() +
        scale_x_continuous(n.breaks = 3, labels = c("", b, "")) +
        facet_wrap(facets = vars(.data$estimator)) +
        scale_color_viridis_d(name = "CV Risk", begin = 0, end = 0.8) +
        labs(
          title = "Estimator Eigenvalues", caption = cv_details,
          x = "Eigenvalue Index", y = "Eigenvalue"
        ) +
        theme_cvCovEst(plot_type = plot_type)
    }
    else {
      if ("summary" %in% plot_type) {
        plot <- ggplot(
          stat_eigs,
          aes(x = .data$index, y = .data$eigenvalues)
        ) +
          geom_path(color = blues[9]) +
          scale_x_continuous(n.breaks = min(10, k)) +
          facet_wrap(facets = vars(.data$estimator)) +
          labs(
            title = "Estimator Eigenvalues", caption = cv_details,
            x = "Eigenvalue Index", y = "Eigenvalue"
          ) +
          theme_cvCovEst(plot_type = plot_type)
      }
      else {
        plot <- ggplot(
          stat_eigs,
          aes(x = .data$index, y = .data$eigenvalues, color = .data$stat)
        ) +
          geom_path() +
          scale_x_continuous(n.breaks = min(10, k)) +
          facet_wrap(facets = vars(.data$estimator)) +
          scale_color_viridis_d(name = "CV Risk", begin = 0, end = 0.8) +
          labs(
            title = "Estimator Eigenvalues", caption = cv_details,
            x = "Eigenvalue Index", y = "Eigenvalue"
          ) +
          theme_cvCovEst(plot_type = plot_type)
      }
    }

    return(plot)
  }
  # Multiple Estimators
  else {
    stat_eigs <- lapply(stat, function(stat_est) {
      est_eigs <- lapply(estimator, function(est) {
        # Has Hypers
        if (est %in% has_hypers) {
          estimator_stats <- cv_sum$hyperRisk[[est]]

          est_hypers <- getHypers(dat = estimator_stats, summ_stat = stat_est)

          # Run The Associated Estimator
          dat <- list(dat_orig)
          arg_names <- c("dat", est_hypers$hyper_names)

          est_args <- append(dat, est_hypers$hyper_values)
          names(est_args) <- arg_names

          estimate <- rlang::exec(est, !!!est_args)
        }
        else {
          # No Hypers
          estimate <- rlang::exec(est, dat_orig)
        }
        # Create Data Frame
        eigs_df <- tibble::tibble(
          index = index,
          eigenvalues = suppressWarnings(RSpectra::eigs_sym(
            estimate,
            k = k, which = eig_type
          )$values),
          stat = rep(stat_est, k),
          estimator = rep(est, k)
        )

        return(eigs_df)
      })

      return(est_eigs)
    })

    # Combine and Re-factor
    stat_eigs <- dplyr::bind_rows(stat_eigs)

    stat_eigs$stat <- factor(
      stat_eigs$stat,
      levels = stat_choices
    )

    stat_eigs$estimator <- factor(
      stat_eigs$estimator,
      levels = cv_sum$bestInClass$estimator
    )

    # Generate Plot
    if (k == 1) {
      plot1 <- ggplot(
        stat_eigs,
        aes(
          x = .data$index,
          y = .data$eigenvalues,
          color = .data$stat
        )
      ) +
        geom_point() +
        scale_x_continuous(n.breaks = 3, labels = c("", b, ""))
    }
    else {
      plot1 <- ggplot(
        stat_eigs,
        aes(
          x = .data$index,
          y = .data$eigenvalues,
          color = .data$stat
        )
      ) +
        geom_path() +
        scale_x_continuous(n.breaks = min(10, k))
    }

    plot <- plot1 +
      facet_wrap(facets = vars(.data$estimator)) +
      scale_color_viridis_d(name = "CV Risk", begin = 0, end = 0.8) +
      labs(
        title = "Estimator Eigenvalues", caption = cv_details,
        x = "Eigenvalue Index", y = "Eigenvalue"
      ) +
      theme_cvCovEst(plot_type = plot_type)

    return(plot)
  }
}


################################################################################
#' Multi-Hyperparameter Risk Plots
#'
#' @description \code{multiHyperRisk()} produces plots of the cross-validated
#'  risk for estimators with more than one hyperparameter.  The function
#'  transforms one of the hyperparameters into a factor and uses it to
#'  distinguish between the risk of various estimators.  If one of the
#'  hyperparameters has only one unique value, that hyperparameter is used as
#'  the factor variable.  If all hyperparameters have only one unique value, a
#'  plot is not generated for that estimator class.
#'
#' @param dat A \code{data.frame} of cross-validated risks. Specifically, this
#'  is the \code{risk_df} table output by \code{\link{cvCovEst}()}.
#' @param estimator A \code{character} vector specifying one or more classes of
#'  estimators to compare.
#' @param switch_vars A \code{logical} indicating if the x-axis and factor
#'  variables should be switched.  Default is \code{FALSE}.
#' @param min_max A \code{logical}. If \code{TRUE}, only the minimum and
#'  maximum values of the factor hyperparameter will be used. Defaults to
#'  \code{FALSE}.
#'
#' @importFrom dplyr filter %>%
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot geom_path aes vars facet_wrap scale_color_viridis_d labs
#'
#' @return A named \code{list} of plots.
#'
#' @keywords internal
multiHyperRisk <- function(
                           dat,
                           estimator,
                           switch_vars = FALSE,
                           min_max = FALSE) {
  plot_list <- list()

  for (e in estimator) {
    p <- switch(e,
      robustPoetEst = plotRobustPoetEst(
        dat = dat,
        switch_vars = switch_vars,
        min_max = min_max
      ),
      poetEst = plotPoetEst(
        dat = dat,
        switch_vars = switch_vars,
        min_max = min_max
      ),
      adaptiveLassoEst = plotAdaptiveLassoEst(
        dat = dat,
        switch_vars = switch_vars,
        min_max = min_max
      )
    )

    plot_list <- append(plot_list, p)
  }

  return(plot_list)
}

################################################################################
#' Cross-Validated Risk Plot
#'
#' @description \code{cvRiskPlot()} plots the cross-validated risk for a given
#'  estimator, or set of estimators, as a function of the hyperparameters.
#'
#' @param dat A named \code{list}.  Specifically, this is the standard output of
#'  \code{\link{cvCovEst}()}.
#' @param estimator A \code{character} vector specifying one or more classes of
#'  estimators to compare.
#' @param plot_type A \code{character} detailing the type of plot. Passed to
#'  \code{\link{theme_cvCovEst}()}, defaults to \code{"risk"}
#' @param cv_details A \code{character} vector summarizing key arguments passed
#'  to \code{\link{cvCovEst}()}.
#' @param switch_vars A \code{logical}. If \code{TRUE},
#'  the hyperparameters used for the x-axis and factor variables are switched.
#'  Only applies to estimators with more than one hyperparameter. Defaults to
#'  \code{FALSE}.
#' @param min_max A \code{logical}. If \code{TRUE}, only the minimum and
#'  maximum values of the factor hyperparameter will be used.  Only applies to
#'  estimators with more than one hyperparameter. Defaults to \code{FALSE}.
#'
#' @return A single plot or grid of plots for each estimator specified.
#'
#' @importFrom stringr str_split
#' @importFrom dplyr group_by filter %>% summarise arrange count
#' @importFrom stats t.test
#' @importFrom assertthat assert_that
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggpubr ggarrange annotate_figure
#' @importFrom ggplot2 ggplot aes vars geom_path geom_ribbon alpha facet_wrap scale_x_continuous labs
#' @importFrom rlang .data
#'
#' @keywords internal
cvRiskPlot <- function(
                       dat,
                       est,
                       plot_type = "risk",
                       cv_details,
                       switch_vars = FALSE,
                       min_max = FALSE) {

  # Get Attributes
  attr_df <- estAttributes(estimator = est)
  attr_df <- dplyr::bind_rows(attr_df)
  attr_df$estimator <- est

  # Check for estimators with no hypers
  assertthat::assert_that(
    any(attr_df$has_hypers),
    msg = "Risk plot unavailable for estimators without hyperparameters."
  )

  if (!all(attr_df$has_hypers)) {
    message("Only plots for estimators with hyperparameters will be displayed")
  }

  # Filter only estimators with hypers
  has_hypers <- attr_df$estimator[which(attr_df$has_hypers)]
  hyper_dat <- dat$risk_df %>%
    dplyr::filter(.data$estimator %in% has_hypers)

  # Estimators with single or  2+ hypers
  single_hyper <- attr_df$estimator[which(attr_df$n_hypers == 1)]
  multi_hypers <- attr_df$estimator[which(attr_df$n_hypers > 1)]

  # Send multi hyper estimators to multiHyperRisk
  if (any(est %in% multi_hypers)) {
    est_names <- est[which(est %in% multi_hypers)]

    multi_plots <- multiHyperRisk(
      dat = hyper_dat,
      estimator = est_names,
      switch_vars = switch_vars,
      min_max = min_max
    )

    # Parse list of plots into a single output
    num_plots <- length(multi_plots)

    if (num_plots == 1) {
      final_plot <- multi_plots[[1]] +
        labs(caption = cv_details)
    }
    if (num_plots > 1 & num_plots <= 3) {
      final_plot <- ggpubr::ggarrange(
        plotlist = multi_plots,
        ncol = num_plots,
        nrow = 1,
        align = "h"
      )

      final_plot <- ggpubr::annotate_figure(
        final_plot,
        bottom = cv_details
      )
    }
    if (num_plots == 4) {
      final_plot <- ggpubr::ggarrange(
        plotlist = multi_plots,
        ncol = 2,
        nrow = 2,
        align = "hv"
      )

      final_plot <- ggpubr::annotate_figure(
        final_plot,
        bottom = cv_details
      )
    }
    if (num_plots > 4) {
      final_plot <- ggpubr::ggarrange(
        plotlist = multi_plots,
        ncol = 3,
        nrow = 2,
        align = "hv"
      )

      final_plot <- ggpubr::annotate_figure(
        final_plot,
        bottom = cv_details
      )
    }
    final_plot1 <- final_plot
  }
  else {
    final_plot1 <- NULL
  }

  # For estimators with only one hyper
  if (any(est %in% single_hyper)) {
    est <- est[which(est %in% single_hyper)]

    risk <- hyper_dat %>%
      dplyr::filter(.data$estimator %in% est)

    # Check for occurrences of only 1 hyperparameter
    hyper_count <- risk %>%
      dplyr::group_by(.data$estimator) %>%
      dplyr::count()

    if (any(hyper_count$n == 1)) {
      removed <- hyper_count$estimator[which(hyper_count$n == 1)]
      keep <- hyper_count$estimator[which(hyper_count$n != 1)]

      assertthat::assert_that(
        length(keep) > 0,
        msg = "Cannot plot estimators with only one hyperparameter instance."
      )

      risk <- risk %>%
        dplyr::filter(.data$estimator %in% keep)

      rem_message <- paste(
        "The following estimators were omitted due to insufficient",
        "observations",
        removed,
        sep = " "
      )

      message(rem_message)
    }

    hyper_vals <- lapply(risk$hyperparameters, function(h) {
      h_split <- stringr::str_split(
        h, "= "
      ) %>% unlist()

      return(as.numeric(h_split[2]))
    }) %>%
      unlist()

    risk$hyperparameters <- hyper_vals
    risk <- risk %>%
      dplyr::group_by(.data$estimator) %>%
      dplyr::arrange(.data$hyperparameters, .by_group = TRUE)

    final_plot2 <- ggplot(risk) +
      geom_path(aes(x = .data$hyperparameters, y = .data$empirical_risk)) +
      facet_wrap(facets = vars(.data$estimator), scales = "free_x") +
      labs(
        title = "Estimator Cross-Validated Risk",
        caption = cv_details,
        x = "Hyperparameter",
        y = "Risk"
      ) +
      scale_x_continuous(n.breaks = 10) +
      theme_cvCovEst(plot_type = plot_type)
  }
  else {
    final_plot2 <- NULL
  }

  plot_list <- list(
    multi_plots = final_plot1,
    single_plots = final_plot2
  )

  return(plot_list)
}

################################################################################
#' Summary Plot
#'
#' @description \code{cvSummaryPlot()} combines plots of the eigenvalues and the
#'  covariance heatmap for the optimal estimator selected by
#'  \code{\link{cvCovEst}()}, and also provides a table showing the best
#'  estimator within each class.  A plot the risk of the optimal estimator's
#'  class is also provided if applicable.
#'
#' @param dat A named \code{list}.  Specifically, this is the standard output of
#'  \code{cvCovEst}.
#' @param estimator A \code{character} vector specifying which class of
#'  estimator to plot.
#' @param dat_orig The \code{numeric data.frame}, \code{matrix}, or similar
#'  object originally passed to \code{cvCovEst}.
#' @param plot_type A \code{character} detailing the type of plot. Passed to
#'  \code{\link{theme_cvCovEst}()}, defaults to \code{"risk"}
#' @param cv_details Character vector summarizing key arguments passed to
#'  \code{\link{cvCovEst}()}.
#' @param has_hypers A \code{character} vector containing the names of current
#'  estimators with hyperparameters.
#' @param multi_hypers A \code{character} vector containing the names of current
#'  estimators with multiple hyperparameters.
#' @param abs_v A \code{logical} determining if the absolute value of the matrix
#'  entries should be used for plotting the matrix heatmap.  Default is
#'  \code{TRUE}.
#' @param switch_vars A \code{logical}. If \code{TRUE},
#'  the hyperparameters used for the x-axis and factor variables are switched.
#'  Only applies to estimators with more than one hyperparameter. Defaults to
#'  \code{FALSE}.
#' @param min_max A \code{logical}. If \code{TRUE}, only the minimum and
#'  maximum values of the factor hyperparameter will be used.  Only applies to
#'  estimators with more than one hyperparameter. Defaults to \code{FALSE}.
#'
#' @return A collection of plots and summary statistics for the optimal
#'  estimator selected by \code{cvCovEst}.
#'
#' @importFrom ggpubr annotate_figure ggtexttable ttheme rownames_style colnames_style tbody_style tbody_add_border ggparagraph tab_add_title ggarrange
#' @importFrom ggplot2 alpha
#'
#' @keywords internal
cvSummaryPlot <- function(
                          dat,
                          estimator,
                          dat_orig,
                          stat,
                          k,
                          leading,
                          plot_type = "summary",
                          cv_details,
                          has_hypers,
                          multi_hypers,
                          abs_v,
                          switch_vars,
                          min_max) {
  if (is.null(k)) {
    k <- ncol(dat_orig)
  }

  # Upper Left
  # If estimator has hypers, use Risk Plot, otherwise use risk by class table
  if (estimator %in% has_hypers) {
    p1 <- cvRiskPlot(
      dat = dat,
      est = estimator,
      plot_type = plot_type,
      cv_details = cv_details,
      switch_vars = switch_vars,
      min_max = min_max
    )

    if (estimator %in% multi_hypers) {
      p1 <- p1$multi_plots
    }
    else {
      p1 <- p1$single_plots
    }
  }
  else {
    p1 <- summary.cvCovEst(object = dat, summ_fun = "cvRiskByClass")

    p1 <- ggpubr::ggtexttable(
      p1$empRiskByClass,
      rows = NULL,
      theme = ggpubr::ttheme(
        base_style = "lBlueWhite",
        base_size = 8
      )
    )

    p1 <- ggpubr::tab_add_title(
      p1,
      text = "Cross-Validated Risk Summary By Class",
      face = "bold",
      size = 10,
      hjust = -0.5,
      padding = unit(2, "line")
    )
  }

  # Upper Right - Eigenvalue Plot
  p2 <- cvEigenPlot(
    dat = dat,
    estimator = estimator,
    dat_orig = dat_orig,
    stat = stat,
    k = k,
    leading = leading,
    plot_type = plot_type,
    cv_details = cv_details,
    has_hypers = has_hypers
  )

  # Lower Left - Heatmap
  p3 <- cvMultiMelt(
    dat = dat,
    estimator = estimator,
    stat = stat,
    dat_orig = dat_orig,
    plot_type = c("summary", "heatmap"),
    cv_details = cv_details,
    has_hypers = has_hypers,
    abs_v = abs_v
  )

  # Lower Right - Summary Table - bestInClass
  p4_a <- summary.cvCovEst(object = dat, summ_fun = "bestInClass")

  best_est <- p4_a$bestInClass$estimator[1]

  colnames(p4_a$bestInClass) <- c("Estimator", "Hyperparameter(s)", "CV Risk")
  p4_a <- ggpubr::ggtexttable(
    p4_a$bestInClass,
    rows = NULL,
    theme = ttheme(
      base_style = "lBlueWhite",
      base_size = 8
    )
  )

  p4_a <- ggpubr::tab_add_title(
    p4_a,
    text = "Best Hyperparameter Values by Class",
    face = "bold",
    size = 10,
    hjust = -0.125,
    padding = unit(2, "line")
  )

  p4 <- ggpubr::annotate_figure(
    p4_a,
    fig.lab = cv_details,
    fig.lab.pos = "bottom.left",
    fig.lab.size = 12,
    fig.lab.face = "italic"
  )

  plot <- ggpubr::ggarrange(
    p1, p2, p3, p4,
    ncol = 2, nrow = 2,
    heights = c(1, 1.2)
  )

  plot <- ggpubr::annotate_figure(
    plot,
    top = paste("cvCovEst Selection: ", best_est, sep = "")
  )

  return(plot)
}



################################################################################
#' Generic Plot Method for cvCovest
#'
#' @description The \code{plot} method is a generic method for plotting objects
#'  of class, \code{"cvCovEst"}.  The method is designed as a tool for diagnostic
#'  and exploratory analysis purposes when selecting a covariance matrix
#'  estimator using \code{cvCovEst}.
#'
#' @param x An object of class, \code{"cvCovEst"}.  Specifically, this is the
#'  standard output of the function \code{cvCovEst}.
#' @param dat_orig The \code{numeric data.frame}, \code{matrix}, or similar
#'  object originally passed to \code{cvCovEst}.
#' @param plot_type A \code{character} vector specifying one of four choices of
#'  diagnostic plots.  Default is \code{"summary"}.  See Details for more about
#'  each plotting choice.
#' @param estimator A \code{character} vector specifying one or more classes of
#'  estimators to compare.  If \code{NULL}, the class of estimator associated
#'  with optimal \code{cvCovEst} selection is used.
#' @param stat A \code{character} vector of one or more summary statistics to
#'  use when comparing estimators.  Default is \code{"min"} for minimum
#'  cross-validated risk.  See Details for more options.
#' @param k A \code{integer} indicating the number of leading/trailing
#'  eigenvalues to plot. If \code{NULL}, will default to the number of columns
#'  in \code{dat_orig}.
#' @param leading A \code{logical} indicating if the leading eigenvalues should
#'  be used.  Default is \code{TRUE}.  If \code{FALSE}, the trailing eigenvalues
#'  are used instead.
#' @param abs_v A \code{logical} determining if the absolute value of the matrix
#'  entries should be used for plotting the matrix heat map.  Default is
#'  \code{TRUE}.
#' @param switch_vars A \code{logical}. If \code{TRUE}, the hyperparameters used
#'  for the x-axis and factor variables are switched in the plot of the
#'  cross-validated risk.  Only applies to estimators with more than one
#'  hyperparameter. Default is \code{FALSE}.
#' @param min_max A \code{logical}.  If \code{TRUE}, only the minimum and
#'  maximum values of the factor hyperparameter will be used.  Only applies to
#'  estimators with more than one hyperparameter. Default is \code{FALSE}.
#' @param ... Additional arguments passed to the plot method.  These are not
#'  explicitly used and should be ignored by the user.
#'
#' @details This plot method is designed to aide users in understanding the
#'  estimation procedure carried out in \code{\link{cvCovEst}()}. There are
#'  currently four different values for \code{plot_type} that can be called:
#'  \itemize{
#'     \item \code{"eigen"} - Plots the eigenvalues associated with the
#'       specified \code{estimator} and \code{stat} arguments in decreasing
#'       order.
#'     \item \code{"risk"} - Plots the cross-validated risk of the specified
#'       \code{estimator} as a function of the hyperparameter values passed to
#'       \code{\link{cvCovEst}()}.  This type of plot is only compatible with
#'       estimators which take hyperparameters as arguments.
#'     \item \code{"heatmap"} - Plots a covariance heat map associated with the
#'       specified \code{estimator} and \code{stat} arguments.  Multiple
#'       estimators and performance stats may be specified to produce grids of
#'       heat maps.
#'     \item \code{"summary"} - Specifying this plot type will run all of the
#'       above plots for the best performing estimator selected by
#'       \code{\link{cvCovEst}()}.  These plots are then combined into a single
#'       panel along with a table containing the best performing estimator
#'       within each class.  If the optimal estimator selected by
#'       \code{\link{cvCovEst}()} does not have hyperparameters, then the risk
#'       plot is replaced with a table displaying the minimum, first quartile,
#'       median, third quartile, and maximum of the cross-validated risk
#'       associated with each class of estimator.
#'  }
#'
#'   The \code{stat} argument accepts five values.  They each correspond to a
#'   summary statistic of the cross-validated risk distribution within a class
#'   of estimator.  Possible values are:
#'   \itemize{
#'     \item \code{"min"} - minimum
#'     \item \code{"Q1"} - first quartile
#'     \item \code{"median"} - median
#'     \item \code{"Q3"} - third quartile
#'     \item \code{"max"} - maximum
#'  }
#'
#' @importFrom assertthat assert_that
#' @importFrom rlang as_label
#'
#' @return A plot object
#'
#' @examples
#' cv_dat <- cvCovEst(
#'   dat = mtcars,
#'   estimators = c(
#'     thresholdingEst, sampleCovEst
#'   ),
#'   estimator_params = list(
#'     thresholdingEst = list(gamma = seq(0.1, 0.9, 0.1))
#'   ),
#'   center = TRUE,
#'   scale = TRUE
#' )
#'
#' plot(x = cv_dat, dat_orig = mtcars)
#' @export
plot.cvCovEst <- function(
                          x,
                          dat_orig,
                          estimator = NULL,
                          plot_type = c("summary"),
                          stat = c("min"),
                          k = NULL,
                          leading = TRUE,
                          abs_v = TRUE,
                          switch_vars = FALSE,
                          min_max = FALSE,
                          ...) {

  # Check cvCovEst credentials
  checkPlotSumArgs(
    dat = x,
    dat_orig = dat_orig,
    which_fun = "plot",
    estimator = estimator,
    plot_type = plot_type,
    stat = stat,
    k = k,
    leading = leading,
    abs_v = abs_v
  )

  # Define cv_details
  pretty_args <- list(
    mc = "Monte Carlo CV",
    v_fold = "V-fold CV",
    cvMatrixFrobeniusLoss = "Matrix Frobenius Loss",
    cvFrobeniusLoss = "Scaled Frobenius Loss",
    cvScaledMatrixFrobeniusLoss = "Scaled Matrix Frobenius Loss"
  )

  if (x$args$cv_scheme == "mc") {
    scheme <- paste(pretty_args$mc, "split", x$args$mc_split, sep = " ")
  }
  else {
    scheme <- pretty_args$v_fold
  }

  folds <- paste(x$args$v_fold, "folds", sep = " ")
  loss <- pretty_args[[rlang::as_label(x$args$cv_loss)]]
  cv_details <- paste(scheme, folds, loss, sep = "  ||  ")

  # Plot the winning estimator for summary plot or NULL estimator
  if (plot_type == "summary" | is.null(estimator)) {
    estimator <- unlist(stringr::str_split(x$estimator, ", "))[1]
  }

  if (is.null(k)) {
    k <- ncol(dat_orig)
  }

  # Get Attributes
  attr_df <- estAttributes(estimator = estimator)
  attr_df <- dplyr::bind_rows(attr_df)
  attr_df$estimator <- estimator
  has_hypers <- attr_df$estimator[which(attr_df$has_hypers)]
  multi_hypers <- attr_df$estimator[which(attr_df$n_hypers > 1)]

  plot <- switch(
    plot_type,
    summary = cvSummaryPlot(
      dat = x,
      estimator = estimator,
      dat_orig = dat_orig,
      stat = stat,
      k = k,
      leading = leading,
      plot_type = "summary",
      cv_details = cv_details,
      has_hypers = has_hypers,
      multi_hypers = multi_hypers,
      abs_v = abs_v,
      switch_vars = switch_vars,
      min_max = min_max
    ),
    risk = cvRiskPlot(
      dat = x,
      est = estimator,
      plot_type = "risk",
      cv_details = cv_details,
      switch_vars = switch_vars,
      min_max = min_max
    ),
    eigen = cvEigenPlot(
      dat = x,
      estimator = estimator,
      stat = stat,
      dat_orig = dat_orig,
      k = k,
      leading = leading,
      plot_type = "eigen",
      cv_details = cv_details,
      has_hypers = has_hypers
    ),
    heatmap = cvMultiMelt(
      dat = x,
      estimator = estimator,
      stat = stat,
      dat_orig = dat_orig,
      plot_type = "heatmap",
      cv_details = cv_details,
      has_hypers = has_hypers,
      abs_v = abs_v
    )
  )

  if (plot_type == "risk") {
    if (is.null(plot$multi_plots)) {
      return(plot$single_plots)
    }
    if (is.null(plot$single_plots)) {
      return(plot$multi_plots)
    }
    return(plot)
  }
  else {
    return(plot)
  }
}
