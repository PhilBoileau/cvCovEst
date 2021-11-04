################################################################################
# Helper Functions for Plotting & Summarizing cvCovEst Results

################################################################################
#' Check for cvCovEst Class
#'
#' @description \code{is.cvCovEst()} provides a generic method for checking if
#' input is of class \code{cvCovEst}.
#'
#' @param x The specific object to test.
#'
#' @return A \code{logical} indicating \code{TRUE} if \code{x} inherits from
#' class \code{cvCovEst}.
#'
#' @examples
#' cv_dat <- cvCovEst(
#'   dat = mtcars,
#'   estimators = c(
#'     thresholdingEst, sampleCovEst
#'   ),
#'   estimator_params = list(
#'     thresholdingEst = list(gamma = seq(0.1, 0.3, 0.1))
#'   ),
#'   center = TRUE,
#'   scale = TRUE
#' )
#'
#' is.cvCovEst(cv_dat)
#' @export
is.cvCovEst <- function(x) {
  inherits(x, "cvCovEst")
}

################################################################################
#' cvCovEst Plot Theme
#'
#' @description \code{theme_cvCovEst()} defines the overall theme of the
#'  \code{cvCovEst} package plotting functions and makes changes depending on
#'  which plot function is being called.
#'
#' @param plot_type A \code{character} vector specifying which plot is to be
#'  displayed. Can contain more than one value in the case of
#'  \code{plot_type = c("heatmap", "summary")}.
#'
#' @return A \code{ggplot} theme.
#'
#' @importFrom ggplot2 theme_gray theme unit element_blank element_text element_rect element_line %+replace%
#' @importFrom RColorBrewer brewer.pal
#'
#' @keywords internal
theme_cvCovEst <- function(plot_type, ...) {
  blues <- RColorBrewer::brewer.pal(9, "Blues")

  # Base Theme
  cv_theme <- theme_gray(...) %+replace%
    theme(
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA),
      panel.grid.major = element_line(color = alpha(blues[3], 0.75)),
      panel.grid.minor = element_line(color = alpha(blues[3], 0.5)),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      plot.title = element_blank(),
      plot.caption = element_text(hjust = 0, size = 10, face = "italic"),
      legend.key = element_blank(),
      legend.title = element_text(vjust = 0.75, size = 10, face = "bold"),
      legend.text = element_text(size = 8, face = "bold"),
      strip.background = element_rect(
        fill = alpha(blues[4], alpha = 0.5), color = blues[9], size = 0.5
      ),
      strip.text = element_text(size = 12, face = "bold", colour = blues[9])
    )

  # Changes for cvMultiMelt
  if ("heatmap" %in% plot_type) {
    cv_theme <- cv_theme +
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(10, "mm")
      )
    # If cvMultiMelt & cvSummaryPlot
    if ("summary" %in% plot_type) {
      cv_theme <- cv_theme +
        theme(
          axis.title.y = element_text(size = 12),
          legend.key.width = unit(10, "mm"),
          legend.title = element_text(
            vjust = 0.75, size = 8, face = "bold"
          )
        )
    }
  }

  # Changes for cvSummaryPlot
  if ("summary" %in% plot_type) {
    cv_theme <- cv_theme +
      theme(
        plot.title = element_text(
          hjust = 0.5, vjust = 1.5, size = 10, face = "bold"
        ),
        plot.caption = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()
      )
  }

  return(cv_theme)
}

################################################################################
#' Hyperparameter Retrieval Function
#'
#' @description \code{getHypers()} retrieves the names and values of all
#'  hyperparameters associated with an estimator passed to \code{cvCovEst()}.
#'
#' @param dat A \code{data.frame} of estimators and their hyperparameter values.
#'  Specifically, this is one of the outputs of \code{summary.cvCovEst()} or
#'  \code{cvCovEst()}.
#' @param summ_stat A character vector specifying the summary statistic of
#'  interest.
#' @param new_df A \code{logical} indicating whether a new \code{data.frame}
#'  should be returned with columns for individual hyperparameters.  Default is
#'  \code{FALSE}.
#'
#' @return A named \code{list} containing the names of all hyperparameters and
#'  their associated values, or a new wider \code{data.frame}.
#'
#' @importFrom rlang exec
#' @importFrom stringr str_split str_squish
#' @importFrom dplyr bind_rows
#'
#' @keywords internal
getHypers <- function(dat, summ_stat, new_df = FALSE) {
  if (new_df) {
    n <- as.integer(ncol(dat))

    hypers <- lapply(dat$hyperparameters, function(h) {
      h_split <- stringr::str_split(
        h, ", "
      ) %>% unlist()

      n <- length(h_split)

      hyper_names <- rep(" ", n)
      hyper_vals <- rep(0, n)

      for (i in 1:n) {
        h_split2 <- stringr::str_split(
          h_split[i], "= "
        ) %>% unlist()

        hyper_names[i] <- stringr::str_squish(h_split2[1])
        if (h_split2[2] %in% c("mad", "sample", "huber")) {
          hyper_vals[i] <- stringr::str_squish(h_split2[2])
        }
        else {
          hyper_vals[i] <- as.numeric(h_split2[2])
        }
      }

      hyper_list <- as.list(hyper_vals)
      names(hyper_list) <- hyper_names

      return(hyper_list)
    })

    hypers <- dplyr::bind_rows(hypers)
    hypers <- cbind(dat[, 1], hypers, dat[, (3:n)])
  }
  else {
    if (is.null(summ_stat)){
      ws <- 1
    }
    else{
      ws <- which(dat$stat == summ_stat)
    }
    hyper_list <- as.list(
      stringr::str_split(dat[ws, "hyperparameters"], ", ") %>% unlist()
    )

    hyper_values <- lapply(hyper_list, function(s) {
      hyper <- stringr::str_split(
        s, "= "
      ) %>% unlist()

      if (hyper[2] %in% c("mad", "sample", "huber")) {
        return(hyper[2])
      }
      else {
        return(as.numeric(hyper[2]))
      }
    })

    hyper_names <- lapply(hyper_list, function(s) {
      hyper <- stringr::str_split(
        s, "= "
      ) %>% unlist()

      return(stringr::str_squish(hyper[1]))
    }) %>%
      unlist()

    hypers <- list(
      hyper_names = hyper_names,
      hyper_values = hyper_values
    )
  }
  return(hypers)
}

################################################################################
#' Estimator Attributes Function
#'
#' @description \code{estAttributes()} returns a named list corresponding to the
#'  attributes of a specific estimator implemented in the \code{cvCovEst}
#'  package.
#'
#' @param estimator A \code{character} vector specifying a class of estimator.
#'
#' @return A named \code{list} containing the attributes of the indicated
#'  estimator.
#'
#' @importFrom rlang exec
#' @importFrom stringr str_split str_squish
#' @importFrom dplyr bind_rows
#'
#' @keywords internal
estAttributes <- function(estimator) {
  est_attrs <- list(
    linearShrinkEst = list(has_hypers = TRUE, n_hypers = 1),
    linearShrinkLWEst = list(has_hypers = FALSE, n_hypers = 0),
    thresholdingEst = list(has_hypers = TRUE, n_hypers = 1),
    sampleCovEst = list(has_hypers = FALSE, n_hypers = 0),
    bandingEst = list(has_hypers = TRUE, n_hypers = 1),
    taperingEst = list(has_hypers = TRUE, n_hypers = 1),
    nlShrinkLWEst = list(has_hypers = FALSE, n_hypers = 0),
    denseLinearShrinkEst = list(has_hypers = FALSE, n_hypers = 0),
    scadEst = list(has_hypers = TRUE, n_hypers = 1),
    poetEst = list(has_hypers = TRUE, n_hypers = 2),
    robustPoetEst = list(has_hypers = TRUE, n_hypers = 3),
    adaptiveLassoEst = list(has_hypers = TRUE, n_hypers = 2)
  )


  out <- lapply(estimator, function(e) {
    return(est_attrs[[e]])
  })

  names(out) <- estimator
  return(out)
}
