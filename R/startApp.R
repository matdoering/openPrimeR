#' The openPrimeR Shiny Application.
#'
#' Starts the openPrimeR Shiny application. 
#' A new tab should open in your default
#' browser. If no browser is opened, please consider the console
#' output to identify the local port on which the server is running
#' and manually open the shown URL.
#'
#' @note
#' The Shiny app can be started only if you fulfill all 
#' of the suggested package dependencies for the Shiny framework,
#' so please ensure that you've installed openPrimeR including
#' all suggested dependencies.
#'
#' @export
#' @return Opens the Shiny app in a web browser.
#' @examples
#' # Start the shiny app
#' \dontrun{
#' startApp()
#' }
startApp <- function() {
    appDir <- system.file("shiny", package = "openPrimeR")
    if (appDir == "") {
        stop(paste("Could not find the directory containing the shiny app: ",
        appDir, "\n", 
        "Try re-installing the 'openPrimeR' package.", call. = FALSE))
    }
    # We need to check for all shiny dependencies here, since they are on the "Suggests" pkg list.
    shiny.dependencies <- c("shiny", "shinyBS",
                            "shinyjs", "DT")
    available.deps <- sapply(shiny.dependencies, function(x)
                        requireNamespace(x))
    if (all(available.deps)) { 
        shiny::runApp(appDir, display.mode = "normal", launch.browser = TRUE)
    } else {
        msg <- paste("Cannot start the shiny App. The following ",
                    "R extensions are missing: ", 
                    paste(shiny.dependencies[!available.deps], collapse = ","),
                    ". Please install them first.", sep = "")
        stop(msg)
    }
} 

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
#' # Open the tutorial
#' \dontrun{
#' runTutorial()
#' }
runTutorial <- function(dev = FALSE) {
    # n.b.: need to manually call to create html files to be included in pkg
    if (dev) {
        # development
        rmarkdown::run("src/openPrimeR/inst/tutorials/introduction/introduction.Rmd")
    } else {
        # release
        learnr::run_tutorial("introduction", package = "openPrimeR")
    }
}

