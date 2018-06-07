################
# Constraint: GC clamp
#################

#' GC clamp
#'
#' Determines the number of consecutive G/Cs at the 3' end. 
#'
#' @param y Pirmer sequence from 5' to 3'.
#'
#' @return The length of the GC clamp.
#' @keywords internal
evaluate.GC.clamp <- function(y) {
    y.iupac <- convert.from.iupac(y)
    n <- 5  # consider the last 5 bases
    i <- NULL
    GC.counts <- foreach(i = seq_along(y.iupac), .combine = c) %dopar% {
        x <- y.iupac[[i]]
        K <- rep(n, length(x))
        K[which(nchar(x) < K)] <- nchar(x)[which(nchar(x) < K)]  # nbr of bases to look at
        tails <- sapply(seq_along(x), function(y) substr(x[y], nchar(x[y]) - K[y] + 
            1, nchar(x[y])))
        # count CG hits
        G.hits <- gregexpr("g", tails)
        C.hits <- gregexpr("c", tails)
        hits <- lapply(seq_along(G.hits), function(x) {
            if (G.hits[[x]][1] == -1 && C.hits[[x]][1] == -1) {
                return(NA)
            } else if (G.hits[[x]][1] == -1) {
                return(C.hits[[x]])
            } else if (C.hits[[x]][1] == -1) {
                return(G.hits[[x]])
            } else {
                h <- c(G.hits[[x]], C.hits[[x]])
                h <- h[order(h)]
            }
        })
        # for disambiguated seqs: take the shortest GC clamp found in any seq
        # reason: if we want to have a gc clamp, we dont want to have an ambig
        # seq that represents a seq that actually doesn't have a gc clamp
        gc.counts <- unlist(lapply(seq_along(hits), 
                   function(x) consecutive.GC.count(hits[[x]], nchar(tails)[x])))
        gc.count <- 0
        if (any(gc.counts > 3)) {
            # more than 3 GCs are undesirable
            gc.count <- max(gc.counts)    
        } else {
            # few GCs (e.g. 0) are undesirable
            gc.count <- min(gc.counts)
        }
    }
    return(GC.counts)
}
#' Consecutive GCs
#'
#' Determines the maximum number of consecutive G/Cs
#'
#' @param y Positions where G/C occurs.
#' Positions are numbered from 1 to 5 where 5 is the end of the primer.
#' @param len Is the number of bases from the primer end considered.
#'
#' @return The maximal number of consecutive G/Cs.
#' @keywords internal
consecutive.GC.count <- function(y, len) {
    if (!any(y %in% len)) {
        # no gc clamp because no G/C at the termianl primer pos
        return(0)  
    }
    diffs <- rle(diff(y))
    # select the last consecutive element
    sel <- length(diffs$values)
    if (sel == 0) {
        return(1)
    }
    if (diffs$values[sel] != 1) {
        return(1)
    }
    count <- diffs$lengths[sel] + 1  # x diffs correspond to x+1 numbers
    return(count)
}
