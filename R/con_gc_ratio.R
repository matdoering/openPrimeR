#########
# Constraint: GC ratio
#########

#' GC ratio
#'
#' Computes the ratio of G/Cs in a sequence. 
#'
#' In case of ambiguities, the mean GC ratio of all possible sequences is computed.
#' 
#' @param x Input sequence.
#' @return The fraction of G/Cs in \code{x}.
#' @keywords internal
compute.gc.ratio <- function(x) {
    x.iupac <- convert.from.iupac(tolower(x))
    result <- rep(0, length(x))
    idx <- which(x != "")
    s <- lapply(seq_along(idx), function(x) strsplit(x.iupac[[idx[x]]], split = ""))
    # check each disambiguated seq
    hits <- lapply(seq_along(s), function(x) sapply(s[[x]], 
                                    function(y) length(which(y == "g" | y == "c")))) 
    if (length(idx) != 0) {
        # compute mean ratio for ambig seqs
        ratios <- sapply(seq_along(idx), function(y) 
                                            mean(hits[[y]]/nchar(x[idx[y]])))
        result[idx] <- ratios
    }
    return(result)
}
