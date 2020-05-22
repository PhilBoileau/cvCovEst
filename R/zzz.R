.onAttach <- function(...) {
  packageStartupMessage(paste0(
    "cvCovEst v", utils::packageDescription("cvCovEst")$Version,
    ": ", utils::packageDescription("cvCovEst")$Title
  ))
}
