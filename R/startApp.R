#' The openPrimeR Tutorial. 
#'
#' Starts a Shiny app containing the openPrimeR tutorial, which
#' was built using the learnr package. The application
#' starts locally and should open a new tab in your default
#' browser. If no browser is opened, please consider the console
#' output to identify the local port on which the server is running.
#'
#' @param dev A logical indicating whether to start the development version of the tutorial (default: \code{FALSE}).
#' @note
#' The Shiny app can be started only if you fulfill all 
#' of the suggested package dependencies for the Shiny framework,
#' so please ensure that you've installed openPrimeR including
#' all suggested dependencies.
#'
#' @export
#' @return Opens the openPrimeR tutorial in a web browser.
#' @examples
#' \dontrun{
#' # Open the tutorial
#' if (interactive()) {
#' runTutorial()
#' }
#' }
runTutorial <- function(dev = FALSE) {
    # n.b.: need to manually call to create html files to be included in pkg
    if (dev) {
        # development
        rmarkdown::run(system.file("inst", "tutorials", "introduction", "introduction.Rmd", package = "openPrimeR"))
    } else {
        # release
        learnr::run_tutorial("introduction", package = "openPrimeR")
    }
}

