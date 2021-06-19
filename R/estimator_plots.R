# Estimator Specific Plotting Functions

################################################################################
#' Plot robustPoetEst
#'
#' @description \code{plotRobustPoetEst()} performs actions specific to plotting
#'  the cross-validated risk of the Robust POET estimator.
#'
#' @param dat A data table of cross-validated risks.  Specifically, this is the
#'  \code{risk_df} table output by \code{\link{cvCovEst}()}.
#' @param switch_vars A \code{logical} indicating if the x-axis and factor
#'  variables should be switched.  Default is \code{FALSE}.
#' @param min_max A \code{logical}. Default is \code{FALSE}. If \code{TRUE},
#'  only the minimum and maximum values of the factor hyperparameter will be
#'  used.
#'
#' @importFrom dplyr filter %>%
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot geom_path aes vars facet_wrap scale_color_viridis_d labs
#'
#' @return A list of plots
#'
#' @keywords internal
plotRobustPoetEst <- function(
                              dat,
                              switch_vars = FALSE,
                              min_max = FALSE) {
  f_dat <- dat %>%
    dplyr::filter(.data$estimator == "robustPoetEst")

  # Split hyperparameters into new columns
  f_dat <- getHypers(dat = f_dat, new_df = TRUE)

  f_dat$estimator <- paste("robustPoetEst -", f_dat$var_est)
  f_dat$lambda <- as.numeric(f_dat$lambda)
  f_dat$k <- as.numeric(f_dat$k)

  plot_env <- new.env()

  plots <- lapply(unique(f_dat$estimator), function(u) {
    # Check for sufficient observations to plot
    u_dat <- f_dat %>%
      dplyr::filter(.data$estimator %in% u)

    n_lam <- length(unique(u_dat$lambda))
    n_k <- length(unique(u_dat$k))

    if (n_lam == 1 & n_k == 1) {
      remove_message <- paste(
        "Omitting risk plot for ", u,
        " . Not enough hyperparameter combinations"
      )

      message(remove_message)
      return(NULL)
    }
    if (n_lam == 1 & n_k > 1) {
      switch_vars <- TRUE
    }
    if (n_lam > 1 & n_k == 1) {
      switch_vars <- FALSE
    }

    if (switch_vars) {
      factor_range <- range(u_dat$lambda)
      u_dat$lambda <- factor(u_dat$lambda)
      u_dat <- dplyr::arrange(u_dat, .data$k)

      if (min_max) {
        u_dat <- u_dat[which(u_dat$lambda %in% factor_range), ]
      }

      plot1 <- ggplot(u_dat) +
        geom_path(
          aes(
            x = .data$k,
            y = .data$cv_risk,
            color = .data$lambda
          )
        ) +
        labs(
          title = "Change in Cross-Validated Risk",
          x = "k", y = "Cross-Validated Risk"
        )
    }
    else {
      factor_range <- range(u_dat$k)
      u_dat$k <- factor(u_dat$k)
      u_dat <- dplyr::arrange(u_dat, .data$lambda)

      if (min_max) {
        u_dat <- u_dat[which(u_dat$k %in% factor_range), ]
      }

      plot1 <- ggplot(u_dat) +
        geom_path(
          aes(
            x = .data$lambda,
            y = .data$cv_risk,
            color = .data$k
          )
        ) +
        labs(
          title = "Change in Cross-Validated Risk",
          x = "lambda", y = "Cross-Validated Risk"
        )
    }

    plot <- plot1 +
      scale_color_viridis_d(begin = 0, end = 0.8) +
      facet_wrap(facets = vars(.data$estimator)) +
      theme_cvCovEst(plot_type = "risk")

    assign(u, plot, envir = plot_env)
  })

  plot_env <- as.list(plot_env)
  return(plot_env)
}

################################################################################
#' Plot poetEst
#'
#' @description \code{plotPoetEst()} performs actions specific to plotting
#'  the cross-validated risk of the POET estimator.
#'
#' @param dat A data table of cross-validated risks.  Specifically, this is the
#'  \code{risk_df} table output by \code{\link{cvCovEst}()}.
#' @param switch_vars A \code{logical} indicating if the x-axis and factor
#'  variables should be switched.  Default is \code{FALSE}.
#' @param min_max A \code{logical}. Default is \code{FALSE}. If \code{TRUE},
#'  only the minimum and maximum values of the factor hyperparameter will be
#'  used.
#'
#' @importFrom dplyr filter %>%
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot geom_path aes vars facet_wrap scale_color_viridis_d labs
#'
#' @return A plot object
#'
#' @keywords internal
plotPoetEst <- function(
                        dat,
                        switch_vars = FALSE,
                        min_max = FALSE) {
  f_dat <- dat %>%
    dplyr::filter(.data$estimator == "poetEst")

  # Split hyperparameters into new columns
  f_dat <- getHypers(dat = f_dat, new_df = TRUE)

  n_lam <- length(unique(f_dat$lambda))
  n_k <- length(unique(f_dat$k))

  if (n_lam == 1 & n_k == 1) {
    remove_message <- paste(
      "Omitting risk plot for poetEst. Not enough hyperparameter combinations"
    )

    message(remove_message)
    return(NULL)
  }
  if (n_lam == 1 & n_k > 1) {
    switch_vars <- TRUE
  }
  if (n_lam > 1 & n_k == 1) {
    switch_vars <- FALSE
  }

  if (switch_vars) {
    factor_range <- range(f_dat$lambda)
    f_dat$lambda <- as.factor(f_dat$lambda)
    f_dat <- dplyr::arrange(f_dat, .data$k)

    if (min_max) {
      f_dat <- f_dat[which(f_dat$lambda %in% factor_range), ]
    }

    plot1 <- ggplot(f_dat) +
      geom_path(
        aes(
          x = .data$k,
          y = .data$cv_risk,
          color = .data$lambda
        )
      ) +
      labs(
        title = "Change in Cross-Validated Risk",
        x = "k", y = "Cross-Validated Risk"
      )
  }
  else {
    factor_range <- range(f_dat$k)
    f_dat$k <- factor(f_dat$k)
    f_dat <- dplyr::arrange(f_dat, .data$lambda)

    if (min_max) {
      f_dat <- f_dat[which(f_dat$k %in% factor_range), ]
    }

    plot1 <- ggplot(f_dat) +
      geom_path(
        aes(
          x = .data$lambda,
          y = .data$cv_risk,
          color = .data$k
        )
      ) +
      labs(
        title = "Change in Cross-Validated Risk",
        x = "lambda", y = "Cross-Validated Risk"
      )
  }

  plot <- plot1 +
    scale_color_viridis_d(begin = 0, end = 0.8) +
    facet_wrap(facets = vars(.data$estimator)) +
    theme_cvCovEst(plot_type = "risk")

  out <- list(poetEst = plot)
  return(out)
}

################################################################################
#' Plot adaptiveLassoEst
#'
#' @description \code{plotAdaptiveLassoEst()} performs actions specific to
#'  plotting the cross-validated risk of the Adaptive LASSO estimator.
#'
#' @param dat A data table of cross-validated risks.  Specifically, this is the
#'  \code{risk_df} table output by \code{\link{cvCovEst}()}.
#' @param switch_vars A \code{logical} indicating if the x-axis and factor
#'  variables should be switched.  Default is \code{FALSE}.
#' @param min_max A \code{logical}. Default is \code{FALSE}. If \code{TRUE},
#'  only the minimum and maximum values of the factor hyperparameter will be
#'  used.
#'
#' @importFrom dplyr filter %>%
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot geom_path aes vars facet_wrap scale_color_viridis_d labs
#'
#' @return A plot object
#'
#' @keywords internal
plotAdaptiveLassoEst <- function(
                                 dat,
                                 switch_vars = FALSE,
                                 min_max = FALSE) {
  f_dat <- dat %>%
    dplyr::filter(.data$estimator == "adaptiveLassoEst")

  # Split hyperparameters into new columns
  f_dat <- getHypers(dat = f_dat, new_df = TRUE)

  # Check for sufficient observations to plot
  n_lam <- length(unique(f_dat$lambda))
  n_n <- length(unique(f_dat$n))

  if (n_lam == 1 & n_n == 1) {
    remove_message <- paste(
      "Omitting risk plot for adaptiveLassoEst.
      Not enough hyperparameter combinations"
    )

    message(remove_message)
    return(NULL)
  }
  if (n_lam == 1 & n_n > 1) {
    switch_vars <- TRUE
  }
  if (n_lam > 1 & n_n == 1) {
    switch_vars <- FALSE
  }
  if (switch_vars) {
    factor_range <- range(f_dat$lambda)
    f_dat$lambda <- factor(f_dat$lambda)
    f_dat <- dplyr::arrange(f_dat, .data$n)

    if (min_max) {
      f_dat <- f_dat[which(f_dat$lambda %in% factor_range), ]
    }

    plot1 <- ggplot(f_dat) +
      geom_path(
        aes(
          x = .data$n,
          y = .data$cv_risk,
          color = .data$lambda
        )
      ) +
      labs(
        title = "Change in Cross-Validated Risk",
        x = "n", y = "Cross-Validated Risk"
      )
  }
  else {
    factor_range <- range(f_dat$n)
    f_dat$n <- factor(f_dat$n)
    f_dat <- dplyr::arrange(f_dat, .data$lambda)

    if (min_max) {
      f_dat <- f_dat[which(f_dat$n %in% factor_range), ]
    }

    plot1 <- ggplot(f_dat) +
      geom_path(
        aes(
          x = .data$lambda,
          y = .data$cv_risk,
          color = .data$n
        )
      ) +
      labs(
        title = "Change in Cross-Validated Risk",
        x = "lambda", y = "Cross-Validated Risk"
      )
  }

  plot <- plot1 +
    scale_color_viridis_d(begin = 0, end = 0.8) +
    facet_wrap(facets = vars(.data$estimator)) +
    theme_cvCovEst(plot_type = "risk")

  out <- list(adaptiveLassoEst = plot)
  return(out)
}
