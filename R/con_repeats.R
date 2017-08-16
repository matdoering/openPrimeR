###################
# Number of repeats
###################
#' Number of Repeats
#' 
#' Computes the number of dinucleotide repeats in the input sequences.
#'
#' @param x Input sequence strings.
#' @return The maximal number of dinucleotide repeats in \code{x}.
#' @keywords internal
nbr.of.repeats <- function(x) {
    # check for dinucleotide repeats, e.g. AT AT AT AT (max nbr of repeats: 4) x:
    # vector of primers (character)
    if (length(x) == 0 || class(x) != "character") {
        return(NULL)
    }
    nts <- c("a", "c", "g", "t")
    options <- expand.grid(nts, nts, stringsAsFactors = FALSE)  # all dinucleotides
    colnames <- c("ID", "Dinucleotide", "nbr_repeats")
    result <- data.frame(matrix(rep(NA, length(colnames) * nrow(options) * length(x)), 
        nrow = length(x) * nrow(options), ncol = length(colnames)), stringsAsFactors = FALSE)
    colnames(result) <- colnames
    # disambiguate primers to deal with amibuity codes
    x.iupac <- convert.from.iupac(x)
    i <- NULL
    results <- foreach(i = seq_len(nrow(options)), .combine = "rbind") %dopar% {
        # for every dinuclotide combination, count the number of occurrences in the seqs
        dinuc <- paste(as.character(options[i, ]), collapse = "")
        expr <- paste("(", dinuc, ")+", sep = "")
        hits <- gregexpr(expr, x.iupac)
        nbr.di <- unlist(lapply(hits, function(y) max(attr(y, "match.length"))))
        nbr.di <- sapply(nbr.di, function(x) if (x == -1) {
            0
        } else {
            x/2  # match length divided by 2 because we want to output the count of dinucleotides, not that of single nucleotides
        })
        t.s <- (i - 1) * length(x) + 1
        t.e <- t.s + length(x) - 1
        di.result <- data.frame(ID = seq_along(x), Dinucleotide = dinuc, nbr_repeats = nbr.di, 
            stringsAsFactors = FALSE)
    }
    repeat.table <- ddply(results, c("ID"), summarize, max_repeat = max(substitute(nbr_repeats)))  # max dinucleotide repeat length for every sequence
    return(repeat.table$max_repeat)
}
