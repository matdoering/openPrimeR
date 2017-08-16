##### FUNCTIONS DEALING WITH AMBIGUITIES

#' Sequence complement
#'
#' Complements the input sequence (re-write of seqinr comp function for gap support)
#' 
#' @param seq Input char vector.
#' @param forceToLower if TRUE the input is transformed to lower case.
#' @param ambiguous if TRUE ambiguous IUPAC nucleotides are complemented.
#' @return The complement of \code{seq}.
#' @keywords internal
comp <- function(seq, forceToLower = TRUE, ambiguous = FALSE) {
      
    if (all(seq %in% LETTERS)) {
        isUpper <- TRUE
    } else {
        isUpper <- FALSE
    }
    seq <- tolower(seq)
    result <- as.vector(seqinr::n2s((3 - seqinr::s2n(seq))))
    if (ambiguous) {
        result[which(seq == "b")] <- "v"
        result[which(seq == "d")] <- "h"
        result[which(seq == "h")] <- "d"
        result[which(seq == "k")] <- "m"
        result[which(seq == "m")] <- "k"
        result[which(seq == "s")] <- "s"
        result[which(seq == "v")] <- "b"
        result[which(seq == "w")] <- "w"
        result[which(seq == "n")] <- "n"
        result[which(seq == "y")] <- "r"
        result[which(seq == "r")] <- "y"
    }
    result[which(seq == "n")] <- "n"
    result[which(seq == "-")] <- "-"
    if (isUpper && !forceToLower) {
        result <- toupper(result)
    }
    return(result)
}

#' Conversion from IUPAC nucleotides
#' 
#' Convert sequences with IUPAC ambiguity codes to all possible sequences without ambiguities.
#'
#' @param seqs A vector of strings.
#'
#' @return lists containing the disambiguated input sequences
#' @keywords internal
convert.from.iupac <- function(seqs) {
    if (length(seqs) == 0) {
        return(NULL)
    }
    # TODO: replace with DECIPHER's function, is much faster
    #warning("DEPRECATD: replace with DECIPHER::Disambiguate for speedup")
    gapped.codemap <- Biostrings::IUPAC_CODE_MAP
    gapped.codemap["-"] <- "-"
    s <- strsplit(seqs, split = "")
    i <- NULL # define iteration variable
    result <- foreach(i = seq_along(seqs), .combine = c) %dopar% {
        if (seqs[i] == "") {
            list(seqs[i])
        } else {
            conv <- sapply(s[[i]], function(x) gapped.codemap[toupper(x)])
            # determine all combinations
            C <- sapply(conv, function(x) strsplit(x, split = ""))
            combis <- expand.grid(C)
            combi.str <- list(apply(combis, 1, function(x) tolower(paste(x, collapse = ""))))
        }
    }
    return(result)
}

#' Merge sequences.
#'
#' Merges the input sequences to one sequence containing IUPAC ambiguity codes.
#'
#' @param seqs Vector of strings.
#'
#' @return Consensus sequence of seqs.
#' @keywords internal
convert.to.iupac <- function(seqs) {
     
    if (length(unique(nchar(seqs))) != 1) {
        warning("Cannot convert to IUPAC. Sequences have different lengths.")
        return(seqs[1])
    }
    N <- nchar(seqs)[1]
    s <- strsplit(toupper(seqs), split = "")
    seq <- sapply(1:N, function(x) paste(sapply(seq_along(s), function(y) s[[y]][x]), 
        collapse = ""))
    res <- paste(tolower(Biostrings::mergeIUPACLetters(seq)), collapse = "")
    return(res)
}

#' Sequence complement
#'
#' Computes the complement of the input sequence.
#'
#' @param seq The input sequence strings.
#'
#' @return The complements of the input sequences.
#' @keywords internal
complement.sequence <- function(seq) {
    idx <- which(seq != "")
    s.orig <- strsplit(seq[idx], split = "")
    complements <- lapply(s.orig, function(x) comp(x, ambiguous = TRUE))
    na.idx <- sapply(complements, function(x) which(is.na(x)))
    for (i in seq_along(na.idx)) {
        j <- na.idx[[i]]
        if (length(j) != 0) {
            s <- complements[[i]]
            s[j] <- s.orig[[i]][j]
            complements[[i]] <- s
        }
    }
    s <- sapply(complements, function(x) paste(x, collapse = ""))
    ret <- rep("", length(seq))
    ret[idx] <- s
    return(ret)
}
#' Reversion of a sequence
#'
#' Reverses the input sequences.
#'
#' @param seq the input sequence.
#'
#' @return The input sequence in reverse order.
#' @keywords internal
rev.sequence <- function(seq) {
      
    s <- sapply(strsplit(seq, split = ""), function(x) paste(rev(x), collapse = ""))
    return(s)
}

#' Reverse complement of a sequence
#'
#' Computes the reverse complement of the input sequences.
#'
#' @param seq the input strings
#'
#' @return The reverse complement of the input sequences.
#' @keywords internal
rev.comp.sequence <- function(seq) {
       s <- sapply(strsplit(seq, split = ""), function(x) paste(rev(comp(x, ambiguous = TRUE)), 
        collapse = ""))
    return(s)
}

#' Indices for merging sequences
#' 
#' Identifies the indices of similar input sequences to be merged.
#'
#' @param seqs The input sequence strings.
#' @param max.degeneracy The maximal allowed degeneracy of a merged seq.
#' 
#'@return A list of lists containing the indices of seqs to be merged. 
#'   For example [[1,2,3]] would indicate to merge primers 1, 2, and 3.
#' @keywords internal
get.merge.idx <- function(seqs, max.degeneracy) {
    s <- convert.from.iupac(seqs)
    merge.list <- NULL
    fw.l <- nchar(seqs)
    fw.lu <- unique(fw.l)
    if (0 %in% fw.lu) {
        # don't merge 0 length entries
        fw.lu <- fw.lu[-which(fw.lu == 0)]
    }
    len.idx <- lapply(seq_along(fw.lu), function(x) which(fw.l == fw.lu[x]))
    names(len.idx) <- fw.lu  # len.idx: all indices of primers of a given length
    merge.indices <- NULL
    for (i in seq_along(len.idx)) {
        # for every primer length
        idx <- len.idx[[i]]  # all primers of current len
        if (length(idx) == 1) {
            next
        }
        rnames <- idx[unlist(lapply(seq_along(s[idx]), function(x) rep(x, length(s[[idx[x]]]))))]  # idx of original primer in df
        p.matrix <- ape::as.DNAbin(ape::as.alignment(strsplit(unlist(s[idx]), split = "")))
        # set rownames of p.matrix
        rownames(p.matrix) <- seq_along(rnames)
        tree <- hclust.tree(p.matrix)
        if (length(tree) != 0) {
            t.seqs <- get.tree.seqs(tree, max.degeneracy, p.matrix)
            if (length(t.seqs) != 0) {
                # determine all possible merges of primers
                merge.idx <- lapply(t.seqs$Merge_Idx, function(x) unique(rnames[as.numeric(strsplit(x, 
                  ",")[[1]])]))
                merge.idx <- merge.idx[sapply(merge.idx, function(x) length(x) > 
                  1)]
                merge.sel <- merge.select(merge.idx)  # select only possible merges
                merge.indices <- c(merge.indices, merge.sel)
            }
        }
    }
    return(merge.indices)
}
#' Select merge indices
#'
#' Greedily identifies the smallest number of possible sequences merges that can be performed.
#'
#' @param merge.idx list of lists containing the indices of possible merges
#'
#' @return The smallest number of possible merge operations as an index list.
#' @keywords internal
merge.select <- function(merge.idx) {
      merge.len <- sapply(merge.idx, length)
    merge.idx <- merge.idx[order(merge.len, decreasing = TRUE)]  # sort by size of merge
    merge.out <- NULL  # vector of lists
    for (i in seq_along(merge.idx)) {
        idx <- merge.idx[[i]]
        if (any(idx %in% unlist(merge.out))) {
            # can't merge these primers
            next
        }
        merge.out <- c(merge.out, list(idx))
    }
    return(merge.out)
}
#' Merge similar primers
#'
#' Merges similar primers contained in the input primer data frame.
#'
#' @param primer.df Primer data frame.
#' @param mode.directionality Analysis direction.
#' @param max.degeneracy Maximal degeneracy of merged primers.
#'
#' @return A primer data frame where similar primers are merged into one entry. 
#' @keywords internal
merge.ambig.primers <- function(primer.df, mode.directionality = c("fw", "rev", "both"), max.degeneracy) {
    if (length(mode.directionality) == 0) {
        stop("'mode.directionality' not supplied.")
    }
    mode.directionality <- match.arg(mode.directionality) 
    if (mode.directionality == "fw") {
        merge.idx <- get.merge.idx(primer.df$Forward, max.degeneracy)
    } else if (mode.directionality == "rev") {
        merge.idx <- get.merge.idx(primer.df$Reverse, max.degeneracy)
    } else {
        # mode both merge only if we can merge both, fw, and rev, or if one of them is
        # not present.
        merge.idx.fw <- get.merge.idx(primer.df$Forward, max.degeneracy)
        merge.idx.rev <- get.merge.idx(primer.df$Reverse, max.degeneracy)
        idx.fw <- which(primer.df$Forward != "")
        idx.rev <- which(primer.df$Reverse != "")
        both <- intersect(idx.fw, idx.rev)
        m.fw <- merge.idx.fw[sapply(merge.idx.fw, function(x) !any(x %in% both))]  # merges of primers ONLY fw
        m.rev <- merge.idx.rev[sapply(merge.idx.rev, function(x) !any(x %in% both))]  # merges of primers ONLY rev
        # dont allow merges of primers containing both directions with primers containing
        # only one direction.  we might miss some merges among primers with both
        # directions, but that's ok.
        b.merges <- NULL # merges for pairs of primers
        if (length(merge.idx.fw) != 0 && length(merge.idx.rev) != 0) {
            both.merges.fw <- which(sapply(merge.idx.fw, function(x) all(x %in% both)))
            both.merges.rev <- which(sapply(merge.idx.rev, function(x) all(x %in% 
                both)))
            # if there's primers to be merged of one direction, check whether they perfectly agree with another pair
            if (length(both.merges.fw) != 0) {
                both.merges <- which(sapply(merge.idx.fw[both.merges.fw], function(x) any(sapply(merge.idx.rev[both.merges.rev], 
                  function(y) all(x == y)))))
                b.merges <- merge.idx.fw[both.merges.fw[both.merges]]  # merges of primers for both directions
            }
        } 
        merge.idx <- c(b.merges, m.fw, m.rev)
    }
    result <- merge.primer.entries(primer.df, merge.idx, mode.directionality)
    return(result)
}
#' Merge similar primers
#'
#' Merges the entries of similar entries in the input primer data frame, given a list with merge indices.
#'
#' @param opti.result Input primer data frame.
#' @param merge.idx List of lists with merge indices (get.merge.idx).
#' @param mode.directionality Direction of primers.
#'
#' @return A primer data frame where entries of similar primers are merged.
#' @keywords internal
merge.primer.entries <- function(opti.result, merge.idx, mode.directionality = c("fw", "rev","both")) {
    if (length(mode.directionality) == 0) {
        stop("'mode.directionality' not supplied.")
    }
    mode.directionality <- match.arg(mode.directionality)   
    if (length(opti.result) == 0 || nrow(opti.result) == 0 || length(merge.idx) == 
        0) {
        # nothing to merge
        return(opti.result)
    }
    new.seqs.fw <- NULL
    new.seqs.rev <- NULL
    if (mode.directionality == "fw") {
        new.seqs.fw <- merge.primer.entries.single(opti.result$Forward, merge.idx)
    } else if (mode.directionality == "rev") {
        new.seqs.rev <- merge.primer.entries.single(opti.result$Reverse, merge.idx)
    } else {
        # both
        new.seqs.fw <- merge.primer.entries.single(opti.result$Forward, merge.idx)
        new.seqs.rev <- merge.primer.entries.single(opti.result$Reverse, merge.idx)
    }
    R.idx <- NULL  # idx of merged rows to be removed after the for loop
    for (i in seq_along(merge.idx)) {
        # idea:retain the first primer.df entry per merge group -> repl.idx
        cur.idx <- merge.idx[[i]]
        repl.idx <- cur.idx[1]
        rm.idx <- cur.idx[2:length(cur.idx)]
        # replace the sequences of forward/reverse
        if (length(new.seqs.fw) != 0) {
            opti.result$Forward[repl.idx] <- new.seqs.fw[i]
        }
        if (length(new.seqs.rev) != 0) {
            opti.result$Reverse[repl.idx] <- new.seqs.rev[i]
        }
        # update ID
        opti.result$ID[repl.idx] <- paste(opti.result$ID[cur.idx], collapse = ";", 
            sep = "")
        # remove the merged rows
        R.idx <- c(R.idx, rm.idx)
    }
    if (length(R.idx) != 0) {
        opti.result <- opti.result[-R.idx, ]
    }
    return(opti.result)
}

#' Merge input sequences
#'
#' Merges the input sequences given a list with merge indices.
#'
#' @param seqs The input sequences.
#' @param merge.idx List of list with merge indices.
#'
#' @return Merged input sequences according to the input merge indices.
#' @keywords internal
merge.primer.entries.single <- function(seqs, merge.idx) {
    new.seqs <- rep("", length(merge.idx))
    for (i in seq_along(merge.idx)) {
        idx <- merge.idx[[i]]
        merge.seqs <- seqs[idx]
        if (any(merge.seqs == "")) {
            # merge not indended for this direction
            next
        }
        s <- convert.to.iupac(merge.seqs)
        new.seqs[i] <- s
    }
    return(as.character(new.seqs))
}

#' Disambiguation of Sequences.
#'
#' Disambiguates the input sequences, but does not disambiguate highly generate sequences.
#' @param template.seqs A \code{DNAStringSet} object with sequences to disambiguate.
#' @param gap.char The character indicating gaps in alignments.
#' @param degen.cutoff The maximal degeneration of sequences to be disambiguated.
#' @return A \code{DNAStringSetList} object with disambiguated sequences.
#' @keywords internal
my.disambiguate <- function(template.seqs, gap.char = "-", degen.cutoff = 2^10) {
    # consider ambiguous templates up to a cutoff of degeneration
    degen <- score_degen(strsplit(tolower(as.character(template.seqs)), split = ""), gap.char = gap.char)
    degen.idx <- which(degen <= degen.cutoff) # cutoff for degeneration
    degen.seqs <- DECIPHER::Disambiguate(template.seqs[degen.idx])
    seqs <- Biostrings::DNAStringSetList(as.list(as.character(template.seqs)))
    seqs[degen.idx] <- degen.seqs
    return(seqs)
}
