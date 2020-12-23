################################################################################
# Plotting Functions for cvCovEst

################################################################################
#' Check for cvCovEst Class
#'
#' @description \code{is.cvCovEst} provides a generic method for checking if
#' input is of class \code{cvCovEst}.
#'
#' @param x The specific object to test.
#'
#' @return A \code{logical} indicating \code{TRUE} if \code{x} inherits from
#' class \code{cvCovEst}.
#'
#' @export
is.cvCovEst <- function(x) {
  inherits(x, "cvCovEst")
}

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
      min = min(empirical_risk),
      Q1 = quantile(empirical_risk, probs = 0.25, type = 3),
      Q2 = quantile(empirical_risk, probs = 0.5, type = 3),
      Q3 = quantile(empirical_risk, probs = 0.75, type = 3),
      max = max(empirical_risk),
      mean_risk = mean(empirical_risk),
      .groups = "keep") %>%
    dplyr::arrange(mean_risk)

  return(empRisk)
}

################################################################################
#' Showing Best Estimator Within Each Class of Estimators
#'
#' @description The \code{bestInClass} function finds the best performing
#'  estimator within each class of estimator passed to \code{cvCovEst} and
#'  finds the associated hyper-parameters if applicable.
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
    "scadEst", "poetEst", "robustPoetEst",
    "adaptiveLassoEst"
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
        row.names = c(
          "min",
          "Q1", "Q2","Q3",
          "max"
          )
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
#' @description The \code{summary} method provides summary statistics regarding
#' the performance of \code{cvCovEst} and can be used for diagnostic plotting.
#'
#' @param dat A named \code{list} of class \code{"cvCovEst"}.
#'
#' @param summ_fun A character vector specifying which summaries to output.
#'
#' @return A named \code{list} where each element corresponds to the output of
#' of the requested summaries.
#'
#' @importFrom rlang exec
#'
#' @export
summary.cvCovEst <- function(
  dat,
  summ_fun = c("empRiskByClass",
               "bestInClass",
               "worstInClass",
               "hyperRisk")) {

  cv_names <- c(
    "estimate",
    "estimator",
    "risk_df",
    "cv_df",
    "args")

  summary_functions <- c(
    "empRiskByClass",
    "bestInClass",
    "worstInClass",
    "hyperRisk"
  )

  if (is.cvCovEst(dat)) {
    assertthat::assert_that(
      all(cv_names %in% names(dat)) == TRUE,
      msg = "cvCovEst object is missing data."
      )
    assertthat::assert_that(
      all(summ_fun %in% summary_functions) == TRUE,
      msg = "Must provide a valid summary function."
    )
  }

  risk_dat <- dat$risk_df

  sums_to_exec <- summary_functions[which(
    summary_functions %in% summ_fun)]

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
#' Hyperparameter Retrieval Function
#'
#' @description The \code{getHypers} retrieves the names and values of all
#'  hyperparameters associated with an estimator passed to \code{cvCovEst}.
#'
#' @param dat A \code{data.frame} of estimators and their hyperparameter values.
#'  Specifically, this is one of the outputs of \code{summary.cvCovEst} or
#'  \code{cvCovEst}.
#'
#' @param summ_stat A character vector specifying the summary statistic of
#'  interest.
#'
#' @param new_df A \code{logical} indicating whether a new \code{data.frame}
#'  should be returned with columns for individual hyperparameters.  Default is
#'  \code{FALSE}.
#'
#' @return A named \code{list} containing the names of all hyperparameters and
#'  their associated values, or a new wider \code{data.frame}.
#'
#' @importFrom rlang exec
#' @importFrom stringr str_split
#'
#' @keywords internal
getHypers <- function(dat, summ_stat, new_df = FALSE) {

  if (new_df){

    hypers <- lapply(
      dat$hyperparameters,
      function(h){
        h_split <- stringr::str_split(
          h, ", "
        ) %>% unlist()

        n <- length(h_split)

        hyperNames <- rep(" ", n)
        hyperVals <- rep(0, n)

        for (i in 1:n) {
          h_split2 <- stringr::str_split(
            h_split[i], "= "
            ) %>% unlist()

          hyperNames[i] <- stringr::str_squish(h_split2[1])
          if (h_split2[2] %in% c('mad', 'sample', 'huber')){
            hyperVals[i] <- stringr::str_squish(h_split2[2])
          }
          else{
            hyperVals[i] <- as.numeric(h_split2[2])
          }
        }

        hyperList <- as.list(hyperVals)
        names(hyperList) <- hyperNames

        return(hyperList)
      })

    hypers <- data.table::rbindlist(hypers)
    hypers <- cbind(dat[,1], hypers, dat[,3])
  }
  else{

  hyperList <- as.list(
    stringr::str_split(
      dat[summ_stat, 1], ", "
    ) %>% unlist()
  )

  hyperValues <- lapply(
    hyperList,
    function(s) {
      hyper <- stringr::str_split(
        s, "= "
      ) %>% unlist()

      if (hyper[2] %in% c('mad', 'sample', 'huber')){
        return(hyper[2])
      }
      else{
        return(
          as.numeric(hyper[2])
        )}
    })

  hyperNames <- lapply(
    hyperList,
    function(s) {
      hyper <- stringr::str_split(
        s, "= "
      ) %>% unlist()

      return(
        stringr::str_squish(hyper[1])
      )
    }) %>% unlist()

  hypers <- list(
    hyperNames = hyperNames,
    hyperValues = hyperValues
  )
  }

  return(hypers)
}


################################################################################
#' Multiple Heat Map Plot
#'
#' @description The \code{cvMultiMelt} compares the structure of two or more
#'  covariance matrix estimators through a grid of heat maps, where each heat map
#'  corresponds to a different estimator.
#'
#' @param dat A named \code{list}.  Specifically, this is the standard output of
#' \code{cvCovEst}.
#'
#' @param estimator A character vector specifying one or more classes of
#' estimators to compare.
#'
#' @param stat A string specifying which summary statistics to use when
#'  comparing two or more estimators.  Default is \code{'min'} for minimum
#'  empirical risk.
#'
#' @param dat_orig The numeric \code{data.frame}, \code{matrix}, or similar
#'  object originally passed to \code{cvCovEst}.
#'
#' @return A grid of heat map plots comparing the desired covariance matrix
#' estimators.
#'
#' @importFrom rlang exec
#' @importFrom reshape2 melt
#' @importFrom stringr str_split
#' @importFrom dplyr filter %>%
#' @import ggplot2
#' @import assertthat
#' @import viridis
#' @import RColorBrewer
#'
#' @keywords internal
cvMultiMelt <- function(dat,
                        estimator,
                        stat = 'min',
                        dat_orig) {


  stat_choices <- c(
    "min",
    "Q1",
    "Q2",
    "Q3",
    "max"
  )

  has_hypers <- c(
    "linearShrinkEst", "thresholdingEst",
    "bandingEst", "taperingEst",
    "scadEst", "poetEst", "robustPoetEst",
    "adaptiveLassoEst"
  )

  single_stat <- ifelse(
    length(stat) == 1,
    TRUE,
    FALSE
  )

  single_est <- ifelse(
    length(estimator) == 1,
    TRUE,
    FALSE
  )

  # Center and Scale Original Data to Match Call to cvCovEst
  dat_orig <- cvCovEst::safeColScale(
    dat_orig,
    center = dat$args$center,
    scale = dat$args$scale
    )

  # Only Certain Summary Stats Supported
  assertthat::assert_that(
    all(stat %in% stat_choices) == TRUE,
    msg = paste(
      "Only the following summary statistics are currently supported: ",
      stat_choices
    )
  )

  # Only estimators called to cvCovEst can be used
  cv_estimators <- unique(dat$risk_df$estimator)
  assertthat::assert_that(
    all(
      estimator %in% cv_estimators == TRUE
    ),
    msg = "Only estimators passed to cvCovEst can be plotted."
  )

  # Call cvSummary
  cv_sum <- summary.cvCovEst(dat)
  blues <- RColorBrewer::brewer.pal(n = 9, name = "Blues")

  # Single Stat Option
  if (single_stat){

    stat_melts <- lapply(
      estimator,
      function(est) {
        # For Estimators With Hyper-parameters
        if (est %in% has_hypers){
          estimatorStats <- cv_sum$hyperRisk[[est]]

          # Get The Associated Hyper-parameters
          estHypers <- getHypers(
            dat = estimatorStats,
            summ_stat = stat
          )

          # Run The Associated Estimator
          dat = list(dat_orig)
          arg_names <- c('dat', estHypers$hyperNames)

          estArgs <- append(
            dat,
            estHypers$hyperValues
            )
          names(estArgs) <- arg_names

          estimate <- rlang::exec(
            est,
            !!!estArgs
          )
        }
        # If No Hyper-parameters
        else{
          estimate <- rlang::exec(
            est,
            dat_orig
          )
        }

        # Create Melted Data Frame
        meltEst <- reshape2::melt(
          abs(estimate)
        )

        # Label by Estimator
        est_name <- rep(
          est,
          nrow(meltEst)
          )

        meltEst$Var1 <- rev(meltEst$Var1)
        meltEst$estimator <- est_name

        return(meltEst)
        })

    # Combine and Re-factor
    stat_melts <- dplyr::bind_rows(stat_melts)
    stat_melts$estimator <- factor(
      stat_melts$estimator,
      levels = cv_sum$bestInClass$estimator
    )

    plot <- ggplot2::ggplot(
      stat_melts,
      aes(x = Var1, y = Var2)) +
      ggplot2::geom_raster(
        aes(fill = value)) +
      ggplot2::facet_wrap(
        facets = vars(estimator)) +
      scale_fill_viridis(
        name = "Covariance (abs. value)",
        option = 'cividis'
      ) +
      theme(legend.position = 'bottom',
            legend.key.width = unit(10, 'mm'),
            legend.title = element_text(
              size = 10,
              face = 'bold',
              vjust = 0.75
            ),
            legend.text = element_text(
              size = 8,
              face = 'bold'),
            strip.background = element_rect(
              fill = alpha(blues[4], alpha = 0.5),
              color = blues[9],
              size = 0.5
            ),
            strip.text = element_text(
              size = 10,
              face = 'bold',
              colour = blues[9]),
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank())

    return(plot)
  }

  # Multi-Stat Option
  else{
    # Multi-Stat for Single Estimator
    if (single_est) {
      assertthat::assert_that(
        estimator %in% has_hypers,
        msg = "The chosen estimator has no hyper-parameters.
        All plots will be the same."
      )

      estimatorStats <- cv_sum$hyperRisk[[estimator]]
      stat_melts <- lapply(
        stat,
        function(stat_est) {
          # Get The Associated Hyper-parameters
          estHypers <- getHypers(
            dat = estimatorStats,
            summ_stat = stat_est
          )

          # Run The Associated Estimator
          dat = list(dat_orig)
          arg_names <- c('dat', estHypers$hyperNames)

          estArgs <- append(
            dat,
            estHypers$hyperValues
          )
          names(estArgs) <- arg_names

          estimate <- rlang::exec(
            estimator,
            !!!estArgs
          )

          # Create Melted Data Frame
          meltEst <- reshape2::melt(
            abs(estimate)
            )


          n <- nrow(meltEst)

          # Label by Stat
          stat_name <- rep(
            stat_est,
            n
          )

          meltEst$Var1 <- rev(meltEst$Var1)
          meltEst$summary_stat <- stat_name

          return(meltEst)
          })

      # Combine and Re-factor
      stat_melts <- dplyr::bind_rows(stat_melts)
      stat_melts$summary_stat <- factor(
        stat_melts$summary_stat,
        levels = stat_choices
      )

      plot <- ggplot2::ggplot(
        stat_melts,
        aes(x = Var1, y = Var2)) +
        ggplot2::geom_raster(
          aes(fill = value)) +
        ggplot2::facet_wrap(
          facets = vars(summary_stat),
          nrow = 1) +
        scale_fill_viridis(
          name = "Covariance (abs. value)",
          option = 'cividis'
        ) +
        theme(legend.position = 'bottom',
              legend.key.width = unit(10, 'mm'),
              legend.title = element_text(
                size = 10,
                face = 'bold',
                vjust = 0.75
              ),
              legend.text = element_text(
                size = 8,
                face = 'bold'),
              strip.background = element_rect(
                fill = alpha(blues[4], alpha = 0.5),
                color = blues[9],
                size = 0.5
              ),
              strip.text = element_text(
                size = 10,
                face = 'bold',
                colour = blues[9]),
              panel.background = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank())



      return(plot)
    }
    else{
      # Multi-Stat for Several Estimators
      multi_melts <- lapply(
        estimator,
        function(est) {
          estimatorStats <- cv_sum$hyperRisk[[est]]

          stat_melts <- lapply(
            stat,
            function(stat_est) {

              # Get The Associated Hyper-parameters
              estHypers <- getHypers(
                dat = estimatorStats,
                summ_stat = stat_est
              )

              # Run The Associated Estimator
              dat = list(dat_orig)
              arg_names <- c('dat', estHypers$hyperNames)

              estArgs <- append(
                dat,
                estHypers$hyperValues
              )
              names(estArgs) <- arg_names

              estimate <- rlang::exec(
                est,
                !!!estArgs
              )

              # Create Melted Data Frame
              meltEst <- reshape2::melt(
                abs(estimate)
                )

              # Label by Stat
              stat_name <- rep(
                stat_est,
                nrow(meltEst)
              )

              # Label by Estimator
              est_name <- rep(
                est,
                nrow(meltEst)
              )

              meltEst$Var1 <- rev(meltEst$Var1)
              meltEst$summary_stat <- stat_name
              meltEst$estimator <- est_name

              return(meltEst)
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

      plot <- ggplot2::ggplot(
        multi_melts,
        aes(x = Var1, y = Var2)) +
        ggplot2::geom_raster(
          aes(fill = value)) +
        ggplot2::facet_grid(
          rows = vars(estimator),
          cols = vars(summary_stat)) +
        scale_fill_viridis(
          name = "Covariance (abs. value)",
          option = 'cividis'
        ) +
        theme(legend.position = 'bottom',
              legend.key.width = unit(10, 'mm'),
              legend.title = element_text(
                size = 10,
                face = 'bold',
                vjust = 0.75
              ),
              legend.text = element_text(
                size = 8,
                face = 'bold'),
              strip.background = element_rect(
                fill = alpha(blues[4], alpha = 0.5),
                color = blues[9],
                size = 0.5
              ),
              strip.text = element_text(
                size = 10,
                face = 'bold',
                colour = blues[9]),
              panel.background = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank())


      return(plot)
    }
  }
}

################################################################################
#' Eigenvalue Plot
#'
#' @description \code{cvEigenPlot} plots the eigenvalues of one or more
#' estimators produced by \code{cvCovEst}.
#'
#' @param dat A named \code{list}.  Specifically, this is the standard output of
#' \code{cvCovEst}.
#'
#' @param estimator A character vector specifying one or more classes of
#' estimators to compare.
#'
#' @param stat A string specifying a single summary statistic to use when
#'  comparing two or more estimators.  Default is \code{'min'} for minimum
#'  empirical risk.  Several stats can be used when only plotting one class of
#'  estimator.
#'
#' @param dat_orig The numeric \code{data.frame}, \code{matrix}, or similar
#'  object originally passed to \code{cvCovEst}.
#'
#' @param k The number of eigenvalues to plot.  Must be less than or equal to
#'  the number of columns of the original data matrix.
#'
#' @param leading A \code{logical} indicating if the leading eigenvalues should
#'  be used.  Defautlt is \code{TRUE}.  If \code{FALSE}, the trailing eigenvalues
#'  will be used instead.
#'
#' @return A plot, or grid of plots, showing the k leading or trailing
#'  eigenvalues of the specified estimators and associated summary statistics of
#'  the empirical risk.
#'
#' @import ggplot2
#' @import viridis
#' @importFrom RColorBrewer brewer.pal
#' @importFrom RSpectra eigs_sym
#' @importFrom rlang exec
#'
#'
#' @keywords internal
cvEigenPlot <- function(
  dat, estimator, stat = 'min', dat_orig, k, leading = TRUE) {

  stat_choices <- c("min", "Q1", "Q2", "Q3", "max")

  has_hypers <- c(
    "linearShrinkEst", "thresholdingEst",
    "bandingEst", "taperingEst",
    "scadEst", "poetEst", "robustPoetEst",
    "adaptiveLassoEst"
  )

  single_stat <- ifelse(
    length(stat) == 1,
    TRUE,
    FALSE
  )

  single_est <- ifelse(
    length(estimator) == 1,
    TRUE,
    FALSE
  )

  # Center and Scale Original Data to Match Call to cvCovEst
  dat_orig <- cvCovEst::safeColScale(
    dat_orig,
    center = dat$args$center,
    scale = dat$args$scale
  )

  # Only Certain Summary Stats Supported
  assertthat::assert_that(
    all(stat %in% stat_choices) == TRUE,
    msg = paste(
      "Only the following summary statistics are currently supported: ",
      stat_choices
    )
  )

  # Only estimators called to cvCovEst can be used
  cv_estimators <- unique(dat$risk_df$estimator)
  assertthat::assert_that(
    all(
      estimator %in% cv_estimators == TRUE
    ),
    msg = "Only estimators passed to cvCovEst can be plotted."
  )

  # Check for appropriate k value
  p <- ncol(dat_orig)
  assertthat::assert_that(
    k <= p,
    msg = 'k must be less than or equal to the number of columns of dat_orig.'
  )

  # Determine leading/trailing and set index accordingly
  eig_type <- ifelse(leading, "LA", "SA")
  if (leading) {
    index <- seq(k)
  }
  else{
    index <- seq(p - k + 1, p)
  }

  # Get Summary Output
  cv_sum <- summary.cvCovEst(dat)
  blues <- RColorBrewer::brewer.pal(9, "Blues")

  # Single Estimator Option
  if (single_est){
    # Has Hypers
    if (estimator %in% has_hypers) {
      estimatorStats <- cv_sum$hyperRisk[[estimator]]

      stat_eigs <- lapply(
        stat,
        function(stat_est) {
          # Get The Associated Hyper-parameters
          estHypers <- getHypers(
            dat = estimatorStats,
            summ_stat = stat_est
          )

          # Run The Associated Estimator
          dat = list(dat_orig)
          arg_names <- c('dat', estHypers$hyperNames)

          estArgs <- append(
            dat,
            estHypers$hyperValues
          )
          names(estArgs) <- arg_names

          estimate <- rlang::exec(
            estimator,
            !!!estArgs
          )

          # Get the eigenvalues
          estEigs <- data.frame(
            index = index,
            eigenvalues = suppressWarnings(
              RSpectra::eigs_sym(
                estimate,
                k = k,
                which = eig_type)$values),
            stat = rep(
              stat_est,
              k
            )
          )

          return(estEigs)
        })

      # Combine and Re-factor
      stat_eigs <- dplyr::bind_rows(stat_eigs)
    }
    # No Hypers
    else{
      if (!single_stat) {
        message(
          "Estimator has no hyperparameters.  All stats yield same estimate."
        )
      }
      estimate <- rlang::exec(
        estimator,
        dat_orig
      )
      # Create Data Frame
      stat_eigs <- data.frame(
        index = index,
        eigenvalues = suppressWarnings(
          RSpectra::eigs_sym(
            estimate,
            k = k,
            which = eig_type)$values),
        stat = rep(
          stat_est,
          k
        )
      )
    }
    # Re-factor
    stat_eigs$stat <- factor(
      stat_eigs$stat,
      levels = stat_choices
    )
    # Generate Plot
    plot <- ggplot2::ggplot(
      stat_eigs,
      aes(x = index,
          y = eigenvalues,
          color = stat)) +
      geom_path() +
      scale_color_viridis_d() +
      xlab("Eigenvalue Index") +
      ylab("Eigenvalue") +
      theme(
        legend.key = element_blank(),
        legend.box.background = element_rect(
          fill = NA,
          size = 0.75
        ),
        legend.title = element_text(
          size = 10,
          face = 'bold',
          vjust = 0.75
        ),
        legend.text = element_text(
          size = 10),
        strip.background = element_rect(
          fill = alpha(blues[4], alpha = 0.5),
          color = blues[9],
          size = 0.5),
        strip.text = element_text(
          size = 10,
          face = 'bold',
          colour = blues[9]),
        axis.title = element_text(
          size = 12),
        axis.text = element_text(
          size = 10),
        plot.title = element_text(
          hjust = 0.5,
          size = 14),
        plot.caption = element_text(
          hjust = 0,
          size = 10,
          face = 'italic'),
        panel.background = element_blank(),
        panel.border = element_rect(
          fill = NA),
        panel.grid = element_line(
          color = alpha(blues[3], 0.75)))

    return(plot)
  }
  # Multiple Estimators
  else{
    stat_eigs <- lapply(
      stat,
      function(stat_est) {

      est_eigs <- lapply(
        estimator,
        function(est) {
          # Has Hypers
          if (est %in% has_hypers){
            estimatorStats <- cv_sum$hyperRisk[[est]]

            estHypers <- getHypers(
              dat = estimatorStats,
              summ_stat = stat_est
            )

            # Run The Associated Estimator
            dat = list(dat_orig)
            arg_names <- c('dat', estHypers$hyperNames)

            estArgs <- append(
              dat,
              estHypers$hyperValues
            )
            names(estArgs) <- arg_names

            estimate <- rlang::exec(
              est,
              !!!estArgs
            )
          }
          # No Hypers
          else{
            estimate <- rlang::exec(
              est,
              dat_orig
            )
          }
          # Create Data Frame
          eigs_df <- data.frame(
            index = index,
            eigenvalues = suppressWarnings(
              RSpectra::eigs_sym(estimate, k = k, which = eig_type)$values
              ),
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
    plot <- ggplot2::ggplot(
      stat_eigs,
      aes(x = index,
          y = eigenvalues,
          color = stat)) +
      geom_path() +
      facet_wrap(facets = vars(estimator)) +
      scale_color_viridis_d(
        name = 'Risk Stats'
      ) +
      xlab("Eigenvalue Index") +
      ylab("Eigenvalue") +
      theme(
        legend.key = element_blank(),
        legend.box.background = element_rect(
          fill = NA,
          size = 0.75
        ),
        legend.title = element_text(
          size = 10,
          face = 'bold',
          vjust = 0.75
        ),
        legend.text = element_text(
          size = 10),
        strip.background = element_rect(
          fill = alpha(blues[4], alpha = 0.5),
          color = blues[9],
          size = 0.5),
        strip.text = element_text(
          size = 10,
          face = 'bold',
          colour = blues[9]),
        axis.title = element_text(
          size = 12),
        axis.text = element_text(
          size = 10),
        plot.title = element_text(
          hjust = 0.5,
          size = 14),
        plot.caption = element_text(
          hjust = 0,
          size = 10,
          face = 'italic'),
        panel.background = element_blank(),
        panel.border = element_rect(
          fill = NA),
        panel.grid = element_line(
          color = alpha(blues[3], 0.75)))
    return(plot)
  }
}

################################################################################
#' Empirical Risk Plot
#'
#' @description \code{cvRiskPlot} plots the empirical risk for a given estimator,
#'  or set of estimators, as a function of the hyperparameters.
#'
#' @param dat A named \code{list}.  Specifically, this is the standard output of
#' \code{cvCovEst}.
#'
#' @param estimator A character vector specifying one or more classes of
#' estimators to compare.
#'
#' @param conf A \code{logical} indicating whether confidence intervals should
#'  be displayed for line plots.  Default is \code{FALSE}.
#'
#' @return A single plot or grid of plots for each estimator specified.
#'
#' @importFrom stringr str_split
#' @importFrom dplyr group_by filter %>%
#' @importFrom stats t.test
#' @import ggplot2
#' @import assertthat
#' @import viridis
#' @importFrom RColorBrewer brewer.pal
#'
#' @keywords internal
cvRiskPlot <- function(dat, est, conf = FALSE) {

  # Exclude estimators with 2+ hyperparameters for now
  invalid_est <- c('poetEst',
                   'adaptiveLassoEst',
                   'robustPoetEst')

  est <- est[which(est %in% invalid_est == FALSE)]

  blues <- RColorBrewer::brewer.pal(9, 'Blues')
  pretty_args <- list(
    mc = 'Monte Carlo CV',
    v_fold = "V-fold CV",
    cvMatrixFrobeniusLoss = "Matrix Frobenius Loss",
    cvFrobeniusLoss = "Scaled Frobenius Loss"

  )

  if (dat$args$cv_scheme == 'mc'){
    scheme <- paste(
      pretty_args$mc, "split", dat$args$mc_split, sep = " "
      )
  }
  else{
    scheme <- pretty_args$v_fold
  }

  folds <- paste(dat$args$v_fold, "folds", sep = " ")
  loss <- pretty_args[[rlang::as_label(dat$args$cv_loss)]]

  cv_details <- paste(scheme, folds, loss, sep = "  ||  ")


  # For Confidence Intervals
  if (conf) {
    # Check for enough fold observations -- Can change to warning/change min folds
    assert_that(
      dat$args$v_folds >= 20,
      msg = 'v_folds must be >= 20 for reliable confidence intervals.')

    risk <- dat$cv_df %>%
      filter(estimator %in% est) %>%
      group_by(estimator, hyperparameters) %>%
      summarise(empirical_risk = mean(loss),
                lower = t.test(loss)$conf.int[1],
                upper = t.test(loss)$conf.int[2],
                .groups = 'keep')

    hyperVals <- lapply(
      risk$hyperparameters,
      function(h) {
        h_split <- stringr::str_split(
          h, "= ") %>% unlist()

        return(as.numeric(h_split[2]))

    }) %>% unlist()

    risk$hyperparameters <- hyperVals
    risk <- risk %>%
      group_by(estimator) %>%
      arrange(hyperparameters, .by_group = TRUE)

    plot <- ggplot(risk) +
      geom_path(aes(x = hyperparameters, y = empirical_risk)) +
      geom_ribbon(aes(x = hyperparameters, ymin = lower, ymax = upper),
                  fill = alpha(blues[8], 0.25)) +
      facet_wrap(facets = vars(estimator),
                 scales = "free_x") +
      labs(title = 'cvCovEst Empirical Risk',
           caption = cv_details) +
      xlab("Hyperparameter") +
      ylab("Risk") +
      theme(strip.background = element_rect(
              fill = alpha(blues[4], alpha = 0.5),
              color = blues[9],
              size = 0.5
            ),
            strip.text = element_text(
              size = 10,
              face = 'bold',
              colour = blues[9]),
            axis.title = element_text(
              size = 12),
            axis.text = element_text(
              size = 10),
            plot.title = element_text(
              hjust = 0.5,
              size = 14
            ),
            plot.caption = element_text(
              hjust = 0,
              size = 10,
              face = 'italic'
            ),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA),
            panel.grid = element_line(color = alpha(blues[3], 0.75)))

    return(plot)

  }
  # No Confidence Intervals
  else{
    risk <- dat$risk_df %>%
      filter(estimator %in% est)

    hyperVals <- lapply(
      risk$hyperparameters,
      function(h) {
        h_split <- stringr::str_split(
          h, "= ") %>% unlist()

        return(as.numeric(h_split[2]))

      }) %>% unlist()

    risk$hyperparameters <- hyperVals
    risk <- risk %>%
      group_by(estimator) %>%
      arrange(hyperparameters, .by_group = TRUE)

    plot <- ggplot(risk) +
      geom_path(aes(x = hyperparameters, y = empirical_risk)) +
      facet_wrap(facets = vars(estimator),
                 scales = "free_x") +
      labs(title = 'cvCovEst Empirical Risk',
           caption = cv_details) +
      xlab("Hyperparameter") +
      ylab("Risk") +
      theme(strip.background = element_rect(
        fill = alpha(blues[4], alpha = 0.5),
        color = blues[9],
        size = 0.5
      ),
      strip.text = element_text(
        size = 10,
        face = 'bold',
        colour = blues[9]),
      axis.title = element_text(
        size = 12),
      axis.text = element_text(
        size = 10),
      plot.title = element_text(
        hjust = 0.5,
        size = 14
      ),
      plot.caption = element_text(
        hjust = 0,
        size = 10,
        face = 'italic'
      ),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA),
      panel.grid = element_line(color = alpha(blues[3], 0.75)))

    return(plot)
  }
}








