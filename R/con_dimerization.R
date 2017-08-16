########
# Constraint: cross dimerization
########

#' Self dimerization
#' 
#' Computes possible self-dimers.
#'
#' @param primers.1 Input primers
#' @param primers.2 (Copy/reverse) of the input primers 
#' @param ions Sodium-equivalent ionic concentration.
#' @param annealing.temp The annealing temperature.
#' @param no.structures Whether the dimerization structure shall be computed.
#' @return Possible self-dimer conformations.
#' @keywords internal
get.self.dimers <- function(primers.1, primers.2, ions, annealing.temp,
                            no.structures = FALSE) {
    # pairs (fw - fw, fw-rev, rev-rev)
    if (length(primers.1) != length(primers.2)) {
        stop("Both primer vectors should have the same length.")
    }
    idx <- which(primers.1 != "" & primers.2 != "")
    dimer.data <- get.dimer.data(primers.1[idx], primers.2[idx],
                                  annealing.temp, ions, no.structures = no.structures)
    return(dimer.data)
}

#' Cross dimers
#'
#' Computes all possible primer cross-dimers.
#' 
#' @param primers.1 Input primers.
#' @param primers.2 Input primers.
#' @param ions Sodium-equivalent ionic concentration.
#' @param annealing.temp The PCR annealing temperature.
#' @param check.idx indices of primers for checking cross-dimerization
#' @param no.structures Whether to compute structures of dimers.
#' @param mode 'symmetric', if \code{primers.1} and \code{primers.2}
#' carry the same information (i.e. fw-fw, rev-rev, fw-rev), 'asymetric' else.
#'
#' @return Data frame with potential cross dimers.
#' @keywords internal
get.cross.dimers <- function(primers.1, primers.2, ions, 
                            annealing.temp, check.idx = NULL, 
                            no.structures = FALSE,
                            mode = c("symmetric", "asymmetric")) {
    mode <- match.arg(mode) 
    if (length(primers.1) == 0 || length(primers.2) == 0) {
        return(NULL)
    }
    #message("Creating dimer combinations ...")
    if (mode == "symmetric") {
        # fewer options for combinations: can exclude self-dimers, as well as already
        # considered combinations (e.g. if we have 1-2, we don't need to 2-1 if 2 in the
        # right set and 2 in the left set correspond to each other).
        combis <- do.call("rbind", parallel::mclapply(seq_along(primers.1), function(x) expand.grid(x, 
            seq_along(primers.2)[-(1:x)])))
    } else {
        # all combinations
        combis <- do.call("rbind", parallel::mclapply(seq_along(primers.1), function(x) expand.grid(x, 
            seq_along(primers.2))))
    }
    if (length(check.idx) != 0) {
        combi.idx <- which(apply(combis, 1, function(x) any(x %in% check.idx)))
        combis <- combis[combi.idx, ]
    }
    idx <- which(primers.1[combis[, 1]] != "" & primers.2[combis[, 2]] != "")  # only consider primer pairings representing actual sequences
    if (length(idx) == 0) {
        return(NULL)
    }
    combis <- combis[idx, ]
    result <- get.dimer.data(primers.1[combis[, 1]], primers.2[combis[, 2]], annealing.temp, ions, no.structures)
    result <- result[, c("DeltaG", "Structure", "Idx1", "Idx2")]
    # modify indices from internal index to used index
    result$Idx1 <- combis[, 1][result$Idx1]
    result$Idx2 <- combis[, 2][result$Idx2]
    return(result)
}

#' Create free energy matrix
#'
#' Creates a matrix giving the deltaG values of all primers.
#'
#' @param primer.df Primer data frame.
#' @param G.df Free energy data for the primers.
#' @param primer.df.2 Optional second primer data frame
#'
#' @return Matrix with the smallest free dimerization energy for every primer pair.
#' @keywords internal
create.G.matrix <- function(primer.df, G.df, primer.df.2 = NULL) {
    # asymmetric: when primer.df.2 represents some other primers that aren't matched to primer.df (during optimization)
    mode <- "asymmetric"  
    if (length(primer.df.2) == 0) {
        primer.df.2 <- primer.df
        mode <- "symmetric"
    }
    if (length(G.df) == 0) {
        # no possible dimerizations detected
        out <- matrix(rep(0, nrow(primer.df) * nrow(primer.df.2)), nrow = nrow(primer.df), 
            ncol = nrow(primer.df.2))
        return(out)
    }
    G <- matrix(rep(0, nrow(primer.df) * nrow(primer.df.2)), nrow = nrow(primer.df), 
        ncol = nrow(primer.df.2))
    # write entries
    G[cbind(G.df$Idx1, G.df$Idx2)] <- G.df$DeltaG
    if (mode == "symmetric") {
        G[cbind(G.df$Idx2, G.df$Idx1)] <- G.df$DeltaG
    }
    return(G)
}

#' Dimerization matrix
#'
#' Computes a matrix indicating all dimerizing primers according to a DeltaG cutoff.
#'
#' @param G Matrix with free energies of all considered primer pairs.
#' @param deltaG.cutoff Primers with free energies below the cutoff are considered dimerizing.
#'
#' @return Binary matrix with dimerization events according to the \code{deltaG.cutoff}.
#' Contains '1' if primers (i,j) dimerize and '0' else.
#' @keywords internal
compute.dimer.matrix <- function(G, deltaG.cutoff = -7) {
    D <- G
    D[G <= deltaG.cutoff] <- 1  # primers i and j dimerize
    D[G > deltaG.cutoff] <- 0  # primer i and j don't dimerize
    # print some infos
    nbr.pairs <- dim(G)[1] * dim(G)[1]
    nbr.dimers <- length(which(D == 1))
    ratio.dimers <- nbr.dimers/(dim(G)[1] * dim(G)[1])
    if (nrow(D) == 0 || all(D == 0)) {
        message("No dimerizing primers detected :-)")
    } else {
        dimerizing.primers <- length(unique(which(D == 1, arr.ind = TRUE)[, 1]))
        cat(paste("Dimerization info: \n\to DeltaG cutoff: ", deltaG.cutoff, "\n\to Number of pairs: ", 
            nbr.pairs, "\n\to Number of dimers: ", nbr.dimers, " (", round(ratio.dimers * 
                100, 2), "%)\n\to Number of dimerizing primers: ", dimerizing.primers, 
            "\n", sep = ""))
    }
    return(D)
}
#' Self dimers
#'
#' Computes all possible self dimers for the primers in the input data frame.
#'
#' @param primer.df Input primer data frame.
#' @param primer_conc Primer concentration.
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @param annealing.temp The PCR annealing temperature.
#' @param for.shiny Whether the output is to be formatted for HTML.
#' @param no.structures Whether dimerization structures shall be outputted.
#' @return Data frame with thermodynamic information on all self dimers.
#' @keywords internal
compute.all.self.dimers <- function(primer.df, primer_conc, na_salt_conc, 
                     mg_salt_conc, k_salt_conc, tris_salt_conc, annealing.temp, for.shiny = FALSE, no.structures = FALSE) {
    # computes data frame with deltaG values of possible primer self-dimers considers
    # all directions: a primer can bind fw-fw, fw-rev, rev-rev
    #message("Computing self dimers @ ", annealing.temp)
    ions <- compute.sodium.equivalent.conc(na_salt_conc, mg_salt_conc, 
                                           k_salt_conc, tris_salt_conc)

    fw.fw <- get.self.dimers(primer.df$Forward, primer.df$Forward, 
                             ions, annealing.temp, no.structures = no.structures)  # fw-fw primer dimers
    if (length(fw.fw) != 0) {
        fw.fw$Direction <- "fw-fw"
    }
    rev.rev <- get.self.dimers(primer.df$Reverse, primer.df$Reverse, 
                               ions, annealing.temp, no.structures = no.structures)  # rev-rev primer dimers
    if (length(rev.rev) != 0) {
        rev.rev$Direction <- "rev-rev"
    }
    fw.rev <- get.self.dimers(primer.df$Forward, primer.df$Reverse, 
                              ions, annealing.temp, no.structures = no.structures)  # fw-rev primer dimers
    if (length(fw.rev) != 0) {
        fw.rev$Direction <- "fw-rev"
    }
    results <- rbind(fw.fw, rev.rev, fw.rev)
    if (length(results) != 0) {
        # remove NAs and add ID
        results <- na.omit(results)
        # results$Idx1 is missing :o
        results$ID <- primer.df$ID[results$Idx1]
    }
    if (for.shiny) {
        results <- view.dimer.df(results, "Self")
    }
    return(results)
}

#' Selection of cross dimer index
#'
#' Select the index with the smallest DeltaG value.
#'
#' @param deltaG Data frame with thermodynamic info.
#' @param primers The corresponding primers.
#'
#' @return The indices with smallest DeltaG for every primer.
#' @keywords internal
select.min.cross.idx <- function(deltaG, primers) {
    # select smallest deltaG idx from data frame of interaction terms, in the order
    # of the primers (NA if nothing is available)
    indices <- seq_along(primers)
    sel.idx <- unlist(lapply(indices, function(x) {
        idx <- which(deltaG$Idx1 == x | deltaG$Idx2 == x)
        if (length(idx) == 0) {
            return(NA)
        }
        min.idx <- which.min(deltaG$DeltaG[idx])
        ret <- idx[min.idx]
    }))
    return(sel.idx)
}
#' Cross dimerization
#'
#' Compute DeltaG data frame for possible primer cross-dimers.
#'
#' @param primer.df Input primers data frame.
#' @param primer_conc Primer concentration.
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @param annealing.temp The PCR annealing temperature.
#' @param check.idx Indices of primers for checking cross-dimerization.
#' @param for.shiny Whether to format for HTML output.
#' @param no.structures Whether dimer structures shall not be determined.
#' If \code{TRUE}, structures are not computed resulting in faster runtimes.
#' @return All cross dimers.
#' @keywords internal
# p.df <- do.call(my_rbind, replicate(30, primer.df, simplify = FALSE))
# X <- compute.all.cross.dimers.unfiltered(p.df, primer_conc, na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc, annealing.temp, check.idx = NULL, for.shiny = FALSE, no.structures = FALSE)
compute.all.cross.dimers.unfiltered <- function(primer.df, primer_conc, na_salt_conc, 
    mg_salt_conc, k_salt_conc, tris_salt_conc, annealing.temp, check.idx = NULL, for.shiny = FALSE, no.structures = FALSE) {
     #TIME <- Sys.time() 
     # ion computation takes some time ...
     ions <- compute.sodium.equivalent.conc(na_salt_conc, mg_salt_conc, 
                                           k_salt_conc, tris_salt_conc)
    # computes all combinations of cross dimers: fw-fw, rev-rev, fw-rev
    # message("Computing cross dimers @ ", annealing.temp)
    # check.idx: only compute cross-dimers for a subset of primers
    fw.fw <- get.cross.dimers(primer.df$Forward, primer.df$Forward, 
                              ions, annealing.temp, check.idx, no.structures)  # fw-fw primer dimers
    if (length(fw.fw) != 0) {
        fw.fw$Direction <- "fw-fw"
    }
    rev.rev <- get.cross.dimers(primer.df$Reverse, primer.df$Reverse, 
                                ions, annealing.temp, check.idx, no.structures)  # rev-primer dimers
    if (length(rev.rev) != 0) {
        rev.rev$Direction <- "rev-rev"
    }
    fw.rev <- get.cross.dimers(primer.df$Forward, primer.df$Reverse, 
                               ions, annealing.temp, check.idx, no.structures)  # fw-rev primer dimers
    if (length(fw.rev) != 0) {
        fw.rev$Direction <- "fw-rev"
    }
    results <- rbind(fw.fw, rev.rev, fw.rev)
    if (length(results) != 0) {
        # remove NAs and add ID
        results <- na.omit(results)
        results$Primer_1 <- as.character(primer.df$ID[results$Idx1])
        results$Primer_2 <- as.character(primer.df$ID[results$Idx2])
    }
    if (for.shiny) {
        results <- view.dimer.df(results, "Cross")
    }
    #message("Time was: ", Sys.time() - TIME)
    return(results)
}
#' Cross dimerization
#'
#' Compute worst-case DeltaG data frame with all possible primer cross-dimers.
#'
#' @param primer.df Input primers data frame.
#' @param primer_conc Primer concentration.
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @param annealing.temp The PCR annealing temperature.
#' @param results (optional) Cross dimer data frame (unfiltered)
#' @param check.idx Indices of primers for checking cross-dimerization.
#' @param for.shiny Whether the table is inteded for HTML display.
#' @param no.structures Whether dimerization structures shall not be outputted.
#' @return Worst-case cross dimers.
#' @keywords internal
compute.all.cross.dimers <- function(primer.df, primer_conc, na_salt_conc, mg_salt_conc, 
    k_salt_conc, tris_salt_conc, annealing.temp, results = NULL, check.idx = NULL, for.shiny = FALSE, no.structures = FALSE) {
    if (length(results) == 0) {
        results <- compute.all.cross.dimers.unfiltered(primer.df, primer_conc, 
                        na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc, 
                        annealing.temp, check.idx, no.structures = no.structures)
    }
    if (length(results) == 0) {
        # no cross dimers found
        return(NULL)
    } 
    # filter out duplicate conformations for the same primer pair
    # replaced plyr call with dplyr call for speed for large matrices ..
    #result <- ddply(results, c("Idx1", "Idx2"), function(x) arrange(x, substitute(DeltaG))[1, ])
    #result <- results %>% dplyr::group_by(.dots = c("Idx1", "Idx2")) %>% dplyr::summarise(DeltaG = min(DeltaG))
    result <- as.data.frame(results %>% dplyr::group_by(.dots = c("Idx1", "Idx2")) %>% dplyr::slice(which.min(DeltaG)))
    if (for.shiny) {
        result <- view.dimer.df(result, "Cross")
    }
    return(result)
}

#' Formats a Dimerization Structure for HTML.
#' @param structures A character vector of dimerization structures.
#' @return HTML-formatted character vectors.
#' @keywords internal
html.format.structure <- function(structures) {
    # check if structures are already formatted
    if (length(structures) == 0) {
        return(structures)
    }
    if (grepl("<div", structures[[1]])) {
        # already formatted!
        return(structures)
    }
    # format structures for non-NA indices only
    idx <- which(!is.na(structures))
    pre <- "<div class = 'verbatim'>"
    post <- "</div>"
    s <- strsplit(structures[idx], "\n")
    s <- unlist(lapply(s, function(x) {
        if (length(x) != 0) {
            paste0("<b>5'-</b>", x[[1]], "<b>-3'</b>", "<br>", 
                paste0("   ", x[[2]]), 
                "<br>", "<b>3'-</b>", x[[3]], "<b>-5'</b>")
        } else {
            ""
        }
    }))
    out <- paste0(pre, s, post)
    result <- rep(NA, length(structures))
    result[idx] <- out
    return(out)
}
#' Formatted dimerization data.
#'
#' Format a dimerization data frame for frontend output.
#'
#' @param dimers Dimerization data frame.
#' @param type Type of dimerization.
#'
#' @return A data frame whose columns are formatted in a user-readable way.
#' @keywords internal
view.dimer.df <- function(dimers, type = c("Self", "Cross")) {
    if (length(type) == 0) {
        stop("Please supply a 'type' argument.")
    }
    type <- match.arg(type)
    out <- dimers
    if (type == "Self") {
        struct.col <- "Self_Dimer_Structure"
    } else {
        struct.col <- "Cross_Dimer_Structure"
    }
    if (!struct.col %in% colnames(out)) {
        # name wasn't changed yet ..
        struct.col <- "Structure"
    }
    structures <- out[, struct.col]
    out[, struct.col] <- html.format.structure(structures)
    rm.cols <- c("Ambig_Index1", "Ambig_Index2") # Idx1, Idx2 retained for shiny
    m <- match(rm.cols, colnames(out))
    m <- m[!is.na(m)]
    if (length(m) != 0) {
        out <- out[,-m]
    }
    #out <- out[, "Primer2" != colnames(out)]  # remove primer2 col
    ## re-write Idx1 to contain the idx of the 'other' primer
    #if (type == "Cross") {
    #    idx <- rep(NA, nrow(primer.df))
    #    for (i in seq_along(primer.df$Forward)) {
    #        if (is.na(out$Idx1[i])) {
    #            next
    #        }
    #        if (out$Idx1[i] == i) {
    #            idx[i] <- out$Idx2[i]
    #        } else {
    #            idx[i] <- out$Idx1[i]
    #            
    #        }
    #    }
    #    out$Idx1 <- primer.df$ID[idx]
    #    out <- out[, "Idx2" != colnames(out)]
    #}
    return(out)
}

#' Self Dimerization.
#'
#' Computes all self dimers in a user-formatted way.
#'
#' @param primer.df Input primer data frame.
#' @param primer_conc Primer concentration.
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @param annealing.temp The PCR annealing temperature.
#' @param for.shiny Whether to format the table for HTML output.
#' @param no.structures Whether dimerization structures shall be outputted.
#' @return A formatted data frame with self-dimerization infos
#' @keywords internal
compute.all.self.dimers.frontend <- function(primer.df, primer_conc, na_salt_conc, 
                                    mg_salt_conc, k_salt_conc, tris_salt_conc, 
                                    annealing.temp, for.shiny = FALSE,
                                    no.structures = FALSE) {
    if (Sys.which("hybrid-min") == "") {
        stop("Cannot compute self dimers without an installation of OligoArrayAux for computing the thermoydnamic properties.")
    }
    self.dimers <- compute.all.self.dimers(primer.df, primer_conc, na_salt_conc, 
                                           mg_salt_conc, k_salt_conc, tris_salt_conc, annealing.temp, for.shiny = for.shiny, no.structures = no.structures)
    # select lowest energy binding mode of a self-dimer:
    self.dimers <- ddply(self.dimers, c("Idx1"), function(x) arrange(x, substitute(DeltaG))[1, ]) # i should retain the self primer index, right?
    self.dimers <- self.dimers[, c("DeltaG", "Structure", "Idx1", "Idx2")]
    if (length(self.dimers) != 0) {
        colnames(self.dimers) <- paste0("Self_Dimer_", colnames(self.dimers))
    }
    if (for.shiny) {
        self.dimers <- view.dimer.df(self.dimers, "Self")
    }
    return(self.dimers)
}
#' Cross Dimerization.
#'
#' Computes all cross dimers in a user-formatted way.
#'
#' @param primer.df Input primer data frame.
#' @param primer_conc Primer concentration.
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @param annealing.temp The PCR annealing temperature.
#' @param for.shiny Whether to format the table for HTML output.
#' @param no.structures Whether to compute structures of dimers.
#' @return A formatted data frame with cross-dimerization infos
#' @keywords internal
compute.all.cross.dimers.frontend <- function(primer.df, primer_conc, na_salt_conc, 
    mg_salt_conc, k_salt_conc, tris_salt_conc, annealing.temp, for.shiny = FALSE, no.structures = FALSE) {
    # computes all combinations of cross dimers: fw-fw, rev-rev, fw-rev
     if (Sys.which("hybrid-min") == "") {
        stop("Cannot compute cross dimers without an installation of OligoArrayAux for computing the thermoydnamic properties.")
    }
    results <- compute.all.cross.dimers.unfiltered(primer.df, primer_conc, na_salt_conc, 
        mg_salt_conc, k_salt_conc, tris_salt_conc, annealing.temp, no.structures = no.structures)
    sel.cols <-  c("DeltaG", "Structure", "Idx1", "Idx2")
    result.cols <- paste0("Cross_Dimer_", sel.cols)
    if (length(results) == 0) {
        # no cross dimers found
        count <- nrow(primer.df)
        out <- data.frame(Cross_Dimer_DeltaG = rep(0, count),
                   Cross_Dimer_Structure = rep("", count),
                   Cross_Dimer_Idx1 = rep(NA, count),
                   Cross_Dimer_Idx2 = rep(NA, count),
                   stringsAsFactors =  FALSE)
        return(out)
    }
    sel.idx <- select.min.cross.idx(results, primer.df$Forward)
    result <- results[sel.idx, sel.cols]
    colnames(result) <- result.cols
    rownames(result) <- NULL
    if (for.shiny) {
        result <- view.dimer.df(result, "Cross")
    }
    return(result)
}
#' Parser for OligoArrayAux Dimerization Data.
#'
#' Parses the free energies and structures of OligoArrayAux.
#'
#' @param deltaG.file A path to a file with OligoArrayAux energies.
#' @param struct.file A path to a file with OligoArrayAux structures.
#' @return A data frame with structures and free energies.
#' @keywords internal
parse.oligo.results <- function(deltaG.file, struct.file) {
    time <- Sys.time()
    results <- try(read.delim(deltaG.file, header = TRUE, 
                stringsAsFactors = FALSE), silent = TRUE)
    if (class(results) == "try-error") {
        warning("Oligo error")
    }
    deltaG <- results[,2]
    #message("DeltaG time: ", Sys.time() - time)
    #time <- Sys.time()
    K.eq <- results[,3] # eq constant; 0 if prefiltered -> no structure!
    prefilter.idx <- which(K.eq == 0) 
    non.prefilter.idx <- which(K.eq != 0)
    if (length(prefilter.idx) == 0) {
        deltaG[prefilter.idx] <- 0 # assume it's no dimer
    }
    if (file.exists(struct.file)) {
        results <- try(readLines(struct.file), silent = TRUE)
        if (class(results) == "try-error") {
            stop("OligoArrayAux: output not produced.")
        }
        last.split <- 1
        split.idx <- which(results == "")
        sel.idx <- setdiff(seq_along(results), split.idx)
        structs <- results[sel.idx]
        structs <- align.structures(structs)
        #message("Struct time: ", Sys.time() - time)
        all.structs <- rep("", length(deltaG))
        if (length(non.prefilter.idx) != 0) {
            all.structs[non.prefilter.idx] <- structs
        }
    } else {
        # non-structure mode
        all.structs <- rep("", length(deltaG))
    }
    out <- data.frame("DeltaG" = deltaG, "Structure" = all.structs, stringsAsFactors = FALSE)
    return(out)
}
#' Identification of Sequence Matches.
#'
#' Identifies matches between two strings
#' provided by OligoArrayAux.
#'
#' @param s1 The aligned nucleotide sequence character vector.
#' @param s2 The aligned, matching substring of \code{s1}.
#' @return A match vector (\code{M} for matches, \code{X} for mismatches).
#' @keywords internal
get.matches <- function(s1, s2) {
    if (any(nchar(s1) != nchar(s2))) { 
        stop("Character vectors of same length required.")
    }
    s1 <- strsplit(s1, "")
    s2 <- strsplit(s2, "")
    match.vector <- lapply(seq_along(s1), function(x) {
        ifelse(s1[[x]] == " ", "X", "M")
    })
    return(match.vector)
}
#' Combination of OligoArrayAux Structure Sequences.
#'
#' Combines the input strings.
#' @param s1 A character vector to be combined with \code{s2}.
#' @param s2 A character vector to be included into \code{s1}.
#' @return A character vector.
#' @keywords internal
combine.strings <- function(s1, s2) {
    s1 <- strsplit(s1, split = "")[[1]]
    s2 <- strsplit(s2, split = "")[[1]]
    idx <- which(s2 != " ")
    out <- s1
    out[idx] <- s2[idx]
    out <- paste0(out, collapse = "")
    return(out)
}
#' Formatting of Dimerization Structures.
#'
#' Formats the given dimerization structures nicely.
#'
#' @param s1 strutcs A character vector, where a block of 4 elements
#' contains: sequence 1 (with removed overlaps), part of sequence 1
#' overlappign with sequence 2, part of sequence 2 
#' overapping with sequence 1, and sequence 2 (with removed overlaps).
#' @return A list of two elements givin the conformation of the first
#' and the second sequence, respectively.
#' @keywords internal
align.structures <- function(structs) {
    # returns a list with two elements: conformation of seq1, conformation of seq2
    # every 
    it <- seq(1, length(structs), 4)
    s1 <- unlist(lapply(seq_along(it), function(x) combine.strings(structs[it[x]], structs[it[x]+1])))
    s2 <- unlist(lapply(seq_along(it), function(x) combine.strings(structs[it[x]+3], structs[it[x]+2])))
    struct.matrix <- do.call(rbind, lapply(seq_along(it), function(x) c(structs[it[x]+1], structs[it[x]+2])))
    match.vector <- get.matches(struct.matrix[,1],struct.matrix[,2])
    match.rep <- unlist(lapply(match.vector, function(x) {
        x[x == "M"] <- "|"
        x[x == "X"] <- " "
        paste(x, collapse = "")
    }))
    out <- paste(s1, match.rep, s2, sep = "\n")
    return(out)
}
expand.grid.unique <- function(x, y, include.equals = TRUE)
{
    x <- unique(x)

    y <- unique(y)

    g <- function(i, x, y)
    {
        z <- setdiff(y, x[seq_len(i-include.equals)])
        if(length(z)) cbind(x[i], z, deparse.level=0)
    }
    if (length(x) <= length(y)) {
        do.call(rbind, lapply(seq_along(x), function(i) g(i, x, y)))
    } else {
        res <- do.call(rbind, lapply(seq_along(y), function(i) g(i, y, x)))
        res[, c(1,2)] <- res[, c(2,1)]
    }
}
#' Preparation of Input for Dimerization.
#' @param s1 Nucleotide character vectors (5' to 3')
#' @param s2 Nucleotide character vectors  (5' to 3')
#' @return A list with two fields containing character vectors.
#' @keywords internal
prepare.dimer.seqs <- function(s1, s2) {
    if (length(s1) == 0 || length(s2) == 0) {
        return(NULL)
    }
    # disambiguate sequences
    s1 <- my.disambiguate(Biostrings::DNAStringSet(s1))
    s2 <- my.disambiguate(Biostrings::DNAStringSet(s2))
    # assign indices
    s1.counts <- IRanges::width(s1@partitioning)
    s2.counts <- IRanges::width(s2@partitioning)
    # indices to indicate the original sequence 
    idx.s1 <- parallel::mclapply(seq_along(s1), function(x) rep(x, s1.counts[x]))
    idx.s2 <- parallel::mclapply(seq_along(s2), function(x) rep(x, s2.counts[x]))
    # index to indicate the position in the list of disambiguated sequences:
    idx.a1 <-  parallel::mclapply(seq_along(s1), function(x) seq_len(s1.counts[x]))
    idx.a2 <-  parallel::mclapply(seq_along(s2), function(x) seq_len(s2.counts[x]))
    # get all combinations
    ambig.options <- parallel::mclapply(seq_along(s1), function(x) {
                    combis <- expand.grid.unique(idx.a1[[x]], idx.a2[[x]])})
    combis.ambig <- do.call(rbind, ambig.options) # index for disambiguated seqs
    colnames(combis.ambig) <- c("Ambig_Index1", "Ambig_Index2")
    combis.ori <- parallel::mclapply(seq_along(s1), function(x) {
                            rep(x, nrow(ambig.options[[x]]))
                })# index for original sequences
    combis.ori <- data.frame(Idx1 = unlist(combis.ori), Idx2 = unlist(combis.ori))
    combis <- cbind(combis.ori, combis.ambig) # TODO: there was a bug here when cbinding for a large set (> 20 mio entries, different number of elements ...)
    s1 <- IRanges::CharacterList(s1)
    s2 <- IRanges::CharacterList(s2)
    # suppress warning that number of columns don't match
    # warning doesn't matter since we select by the correct Ambig_Index
    s1.mat <- suppressWarnings(do.call(rbind, as.list(s1)))
    s2.mat <- suppressWarnings(do.call(rbind, as.list(s2)))
    # this extraction part still takes away some time since we need to iterate over the full list of combinations to select the elements we want
    S1 <- tolower(s1.mat[as.matrix(combis[, c("Idx1", "Ambig_Index1")])])
    S2 <- tolower(s2.mat[as.matrix(combis[, c("Idx2", "Ambig_Index2")])])
    out <- list("Seqs1" = S1, "Seqs2" = S2, "combis" = combis)
    return(out)
}
#' Simple Batchification
#' @param tasks The tasks to assign to individual batches.
#' @return A list of lists containing indices corresponding to \code{tasks}, each list gives a batch.
#' @keywords internal
batchify.simple <- function(tasks) {
    tasks.per.job <- ceiling(length(tasks)/foreach::getDoParWorkers())
    no.jobs <- ceiling(length(tasks)/tasks.per.job)
    total.tasks <- 0
    batches <- vector("list", no.jobs)  # idx of the primers in out.primers per file
    for (i in seq_len(no.jobs)) {
        tasks.added <- min(tasks.per.job, length(tasks) - total.tasks)
        s <- NULL
        if (i == 1) {
            s <- 1
        } else {
            s <- max(batches[[i - 1]]) + 1
        }
        e <- s + (tasks.added - 1)
        batches[[i]] <- s:e
        total.tasks <- total.tasks + tasks.added
    }
    return(batches)
}
#' Batchification by Temperature.
#' @param tasks The tasks to assign to individual batches.
#' @param annealing.temp The annealing temperatures according to which batches are to be created.
#' @return A list of lists containing indices corresponding to \code{tasks}, each list gives a batch.
#' @keywords internal
batchify.temp <- function(tasks, annealing.temp) {
    # range of sampled annealing temperatures for batchification
    temps <- seq(min(annealing.temp), max(annealing.temp), 2)
    # assign every task to a temp
    temp.idx <- sapply(annealing.temp, function(x) which.min(abs(x - temps)))
    job.temps <- unique(temps[temp.idx])
    no.jobs <- length(job.temps)
    batches <- vector("list", no.jobs)  # idx of the primers in out.primers per file
    names(batches) <- job.temps
    for (i in seq_len(no.jobs)) {
        cur.temp.idx <- which(job.temps[i] == temps)
        cur.temp <- job.temps[i]
        batches[[as.character(cur.temp)]] <- tasks[temp.idx == cur.temp.idx]
    }
    return(batches)
}
#' Creates multiple Batches for Parallelization.
#' @param tasks An integer vector with indices representing individual computations.
#' @param annealing.temps Temperatures according to which to batchify.
#' @return A list of lists containing indices corresponding to \code{tasks}, each list gives a batch.
#' @keywords internal
batchify <- function(tasks, annealing.temps = NULL) {
    if (length(tasks) == 0) {
        # nothing to batchify
        return(NULL)
    }
    if (length(annealing.temps) == 0) {
       batches <- batchify.simple(tasks)
    } else {
        batches <- batchify.temp(tasks, annealing.temps)
    } 
    return(batches)
}
#' Retrieval of dimerization energies.
#'
#' Uses OligoArrayAux to compute dimerization candidates.
#'
#' @param s1 Nucleotide character vectors (5' to 3')
#' @param s2 Nucleotide character vectors  (5' to 3')
#' @param annealing.temp The PCR annealing temperature in Celsius.
#' @param ions The sodium-equivalent ions used in the PCR.
#' @param no.structures Whether to compute structures of dimers.
#' @return A data frame containing free energies in the field
#' \code{DeltaG} and the dimerization structure in \code{Structure}.
#' @keywords internal
get.dimer.data <- function(s1, s2, annealing.temp, ions, no.structures) {
    #print("preparing dimer seqs")
    seq.data <- prepare.dimer.seqs(s1, s2)
    #print("preparation done")
    s1 <- seq.data$Seqs1
    s2 <- seq.data$Seqs2
    combis <- seq.data$combis
    if (length(s1) == 0 || length(s2) == 0) {
        # nothing to be computed ..
        return(NULL)
    }
    if (length(s1) != length(s2)) {
        stop("Cannot compute dimerization: s1 and s2 have different lengths.")
    }
    # parallelize -> split into batches
    if (length(annealing.temp) == 1)  {
        # split up seqs independent of annealing temp into batches for parallelization
        batches <- batchify(seq_along(s1))
        annealing.temp <- rep(annealing.temp, length(batches))
    } else {
        # split up seqs dependent on annealing temp
        batches <- batchify(seq_along(s1), annealing.temp)
		# adjust to batch temperatures:
		annealing.temp <- as.numeric(names(batches))
    }
    i <- NULL
    #print("TA")
    #print(annealing.temp)
    #for (i in seq_along(batches)) {
    #print("Computing free energies ...")
    results <- foreach(i = seq_along(batches), .combine = rbind) %dopar% {
        batch <- batches[[i]]
        f1 <- tempfile(pattern = "oligo_dimers_1_", fileext = ".txt")
        f2 <- tempfile(pattern = "oligo_dimers_2_", fileext = ".txt")
        out.prefix <- tempfile(pattern = "oligo_dimers_out_", fileext = "")
        # select sequences
        my.s1 <- s1[batch]
        my.s2 <- s2[batch]
        # write sequences to files
        seqinr::write.fasta(as.list(my.s1), as.list(names(my.s1)), f1)
        seqinr::write.fasta(as.list(my.s2), as.list(names(my.s2)), f2)
        if (length(annealing.temp[i]) == 0 || is.na(annealing.temp[i])) {
            warning("Dimerization: Annealing temperature was NA. Using setting of 50.")
            annealing.temp[i] <- 50
        }
        # filter doesn't speed up the computations so much, just a little bit
        # so let's use the default
        filter <- 2
        # hybrid-min -n DNA -t 50 -T 50 -N 0.1778085 -q tccttcctcatcttcctg caggaggaggaagaacca
        call <- paste0("hybrid-min -n DNA -t ",
                    annealing.temp[i],
                    " -T ",
                    annealing.temp[i],
                    " -N ",
                    ions,
                    " -o ",
                    out.prefix,
                    " --prefilter=",
                    filter, " ")
        #print(call)
        if (no.structures) {
            call <- paste(call, "-E")
        }
        call <- paste(call, f1, f2)
        #message(call)
        system(call, ignore.stdout = TRUE)
        # retrieve results
        deltaG.file <- paste0(out.prefix, ".dG")
        struct.file <- paste0(out.prefix, ".asc")
        results <- parse.oligo.results(deltaG.file, struct.file)
        unlink(c(f1, f2)) 
        results
    }
    # annotate with indices of combinations and re-order from batch order
    results <- cbind(results[match(seq_along(s1), unlist(batches)), ], combis)
    return(results)
}
