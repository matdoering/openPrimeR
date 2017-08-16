########
# Custom error handling functions
#########
#' Custom Error
#' 
#' Creates an error with a custom class.
#'
#' @param subclass String giving the specific type of error.
#' @param message Message to be displayed to the user.
#' @param call Environment where the error ocurred.
#' @param ... Other arguments to be passed to  the condition function.
#' @return Generates a custom error.
#' @keywords internal
my.error <- function(subclass, message, call = sys.call(-1), ...) {
    c <- condition(c(subclass, "openPrimeR_Err", "error"), message, call = call, ...)
    stop(c)
}
#' Custom Warning.
#' 
#' Creates a warning with a custom class.
#'
#' @param subclass String giving the specific type of error.
#' @param message Message to be displayed to the user.
#' @param call Environment where the error ocurred.
#' @param ... Other arguments to the condition function.
#' @return Generates a custom warning.
#' @keywords internal
my.warning <- function(subclass, message, call = sys.call(-1), ...) {
  c <- condition(c(subclass, "openPrimeR_Err", "warning"), message, call = call, ...)
  warning(c)
}
#' Condition Constructor
#'
#' Constructs a condition for custom errors.
#'
#' @param subclass String giving the specific error.
#' @param message String giving the user message.
#' @param call Environment object.
#' @param ... Other arguments for the output structure.
#' @return A condition structure.
#' @keywords internal
condition <- function(subclass, message, call = sys.call(-1), ...) {
    # subclass: inherit from warning/error/message
  structure(
    class = c(subclass, "condition"), # always inherit from condition
    list(message = message, call = call),
    ...
  )
}


