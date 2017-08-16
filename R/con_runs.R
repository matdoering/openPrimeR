###########
# Number of runs
###########

#' Number of Runs
#'
#' Computes the longest run of a single character in the input sequence.
#' 
#' @param x Primer character sequences.
#'
#' @return The longest repeat of a single character in \code{x}.
#' @keywords internal
nbr.of.runs <- function(x) {
    if (length(x) == 0 || class(x) != "character") {
        return(NULL)
    }
    # disambiguate primers
    x.iupac <- convert.from.iupac(x)
    result <- rep(0, length(x))
    idx <- which(x != "")
    s <- lapply(x.iupac[idx], function(y) strsplit(y, split = ""))
    i <- NULL
    runs <- foreach(i = seq_along(idx), .combine = c) %dopar% {
        r <- max(sapply(s[[i]], function(y) max(rle(y)$lengths)))
    }
    if (length(idx) != 0) {
        result[idx] <- runs
    }
    return(result)
}
