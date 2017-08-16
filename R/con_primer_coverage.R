#########
# Primer coverage functions
#########

#' Unique Coverage Indices
#"
#' Computes the indices of templates that are covered uniquely 
#' covered by an individual primer.
#'
#' @param primer.df Primer data frame.
#' @param template.df Template data frame.
#'
#' @return Index of templates uniquely covered per input primer.
#' @keywords internal
compute.unique.covered.idx <- function(primer.df, template.df) {
    # computes the unique coverage of a primer (the number of templates that are not
    # covered by the other primers in the set)
    cvd.idx <- covered.seqs.to.idx(primer.df$Covered_Seqs, template.df)
    unique.idx <- lapply(seq_along(cvd.idx), function(x) setdiff(cvd.idx[[x]], intersect(cvd.idx[[x]], 
        unlist(cvd.idx[-x]))))
    return(unique.idx)
}

#' Conversion of Coverage Strings to Indices.
#'
#' Converts the input coverage strings (comma separated template Identifiers) into indices.
#'
#' @param covered.seqs Strings of covered sequences to be converted.
#' @param template.df Template data frame containing the identifiers of templates
#'
#' @return Indices of covered templates.
#' @keywords internal
covered.seqs.to.idx <- function(covered.seqs, template.df) {
    #print("covered.seqs.to.idx")
    #print(covered.seqs)
    idx <- lapply(strsplit(as.character(covered.seqs), split = ","), 
        function(x) {
            if (all(is.na(x))) {
                NULL
            } else {
                match(as.numeric(x), template.df$Identifier)
            }
    })
    return(idx)
}
#' Conversion of Template Coverage Indices to ID string
#'
#' Converts the input coverage indices to a comma-separated string with the template IDs.
#'
#' @param covered.seqs Indices of covered template sequences.
#' @param template.df Template data frame.
#' @return A string containing the covered template IDs.
#' @keywords internal
covered.seqs.to.ID.string <- function(covered.seqs, template.df) {
    if (length(covered.seqs) == 0) {
        return("")
    }
    ID <- lapply(strsplit(as.character(covered.seqs), split = ","), function(x) paste(template.df[match(as.numeric(x), 
        template.df$Identifier), "ID"], collapse = ","))
    return(ID)
}
#' Conversion of Primer Indices to ID string
#'
#' Converts the input coverage indices to a comma-separated string with the template IDs.
#'
#' @param covered.primers Identifiers of primers covering sequences.
#' @param primer.df Primer data frame.
#' @return String containing the covered template IDs.
#'
#' @return A string containing the IDs of covering primers.
#' @keywords internal
covered.primers.to.ID.string <- function(covered.primers, primer.df) {
    if (length(covered.primers) == 0) {
        return("")
    }
    ID <- lapply(strsplit(as.character(covered.primers), split = ","), function(x) paste(primer.df[match(x, 
        primer.df$Identifier), "ID"], collapse = ","))
    return(ID)
}
#' Covered Templates
#'
#' Get the indices of covered templates.
#'
#' @param Tm.set Primer data frame.
#' @param template.df Template data set.
#'
#' @return Index of templates that are covered by the primers in \code{Tm.set}.
#' @keywords internal
get.covered.templates <- function(Tm.set, template.df) {
    # for each template, return the idx of covering primers Tm.set: primer data frame
    # template.df: template data frame
    cvd.IDs <- strsplit(Tm.set$Covered_Seqs, split = ",")
    cvd.idx <- covered.seqs.to.idx(Tm.set$Covered_Seqs, template.df)
    template.coverage <- lapply(seq_along(template.df$Identifier), function(j) which(sapply(cvd.IDs, 
        function(x) template.df$Identifier[j] %in% x)))
    return(template.coverage)
}
#' 3' Mismatch Check.
#'
#' Check for mismatches at primer 3' ends.
#'
#' @param template.df Template data frame.
#' @param primer.df Primer data frame.
#' @param mode.directionality Primer directionality.
#' @return Returns the distance of mismatches from the 3' terminal end of primers. 
#' @keywords internal
check.3prime.mismatches <- function(template.df, primer.df, 
                           mode.directionality = c("fw", "rev", "both")) {
     if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' argument.")
    }
    mode.directionality <- match.arg(mode.directionality)
    fw.idx <- primer.df$Forward != ""
    rev.idx <- primer.df$Reverse != ""
    if (mode.directionality == "fw") {
        mismatches <- get.3prime.mismatch.pos(primer.df$Forward, primer.df$Mismatch_pos_fw)
    } else if (mode.directionality == "rev") {
        mismatches <- get.3prime.mismatch.pos(primer.df$Reverse, primer.df$Mismatch_pos_rev)
    } else {
        fw.mismatches <- get.3prime.mismatch.pos(primer.df$Forward, primer.df$Mismatch_pos_fw)
        rev.mismatches <- get.3prime.mismatch.pos(primer.df$Reverse, primer.df$Mismatch_pos_rev)
        mismatches <- vector("list", length(fw.mismatches))
        for (x in seq_along(fw.mismatches)) {
            if (fw.idx[x] && rev.idx[x]) {
                res <- sapply(seq_along(fw.mismatches[[x]]), function(y) min(c(fw.mismatches[[x]][[y]], 
                  rev.mismatches[[x]][[y]])))
            } else if (fw.idx[x]) {
                res <- fw.mismatches[[x]]
            } else {
                res <- rev.mismatches[[x]]
            }
            mismatches[[x]] <- res
        }
    }
    return(mismatches)
}
#' 3' Hexamer Check.
#'
#' Check whether the 3' hexamer of a primer is fully complementary
#' to the corresponding region in the template.
#'
#' @param template.df Template data frame.
#' @param primer.df Primer data frame.
#' @param mode.directionality Primer directionality.
#' @return Returns \code{TRUE} if the 3' hexamer of a primer is fully complementary
#' to the corresponding template region and \code{FALSE} otherwise.
#' @keywords internal
check.3prime.hexamers <- function(template.df, primer.df, 
                           mode.directionality = c("fw", "rev", "both")) {
    term.pos <- check.3prime.mismatches(template.df, primer.df, 
                           mode.directionality = mode.directionality)
    hexamer.ok <- lapply(term.pos, function(x) x > 6) # TRUE, if hexamer is free of mismatches
    return(hexamer.ok)
}

#' Conversion of Mismatch Postions String to List.
#'
#' @param mismatches A character vector where parenthesis give mismatches for a template binding event.
#' @return A list with the mismatches for every template for every primer.
#' @keywords internal
mismatch.string.to.list <- function(mismatches) {
    pattern <- "(?<=\\()(([0-9]+,*)*(?=\\)))"
    m <- regmatches(mismatches, gregexpr(pattern, mismatches, perl = TRUE))
    mm <- lapply(m, function(x) lapply(strsplit(x, split = ","), function(y) {
        if (length(y) == 0) 
            NA else as.numeric(y)
    }))
    return(mm)
}
#' Identification of 3' Mismatches.
#'
#' Computes the lastmost position of a 3' mismatches of a primer with a template.
#'
#' @param primers Primer sequence strings.
#' @param mismatches Comma-separated strings containing the primer mismatch positions.
#'
#' @return The closest position of a mismatch relative to the 3' end of the primer.
#' Here, 1 indicates the terminal position, 2 the penultimate position, and so on.
#' No mismatch is indicated by an infinite value.
#' @keywords internal
get.3prime.mismatch.pos <- function(primers, mismatches) {
    # mismatches: mismatches string vector from primer.df primers: primer character
    # vector
    mm <- mismatch.string.to.list(mismatches)
    # no mismatch -> position is Inf (needed for filtering)
    res <- lapply(seq_along(mm), function(x) 
        if (length(mm[[x]]) != 0) {
            # select primers with mismatches:
            idx <- which(sapply(mm[[x]], function(y) !is.na(y[1])))
            r <- rep(Inf, length(mm[[x]]))
            r[idx] <- sapply(mm[[x]][idx], function(y) nchar(primers[x]) - max(y) +1 ) # worst-case position (furthest mismatch 3' prime)
            r
        } 
    ) 
    return(res)
}
#' Update Coverage Information.
#'
#' Updates the coverage-related columns in the input primer data frame.
#' Does not modify the entries of template-specific coverage columns
#' such as primer efficiency (comma-separated values).
#'
#' Removes all coverage events of templates whose index is not in \code{sel}.
#'
#' @param filtered.df Primer data frame.
#' @param sel List with indices of covered templates to be retained, one list with template indices to keep per primer.
#' @param template.df Template data frame.
#' @param mode Either \code{on_target} to filter on-target binding events
#' or \code{off_target} to filter off-target binding events. The corresponding
#' \code{sel} argument should be different.
#' @param active.constraints The active coverage constraints.
#' @return A primer data frame with updated coverage information.
#' @keywords internal
update.cvg.data <- function(filtered.df, sel, template.df, mode = c("on_target", "off_target"), active.constraints) {
    # updates all cvg-related columns according to the 'sel' list (covered templates
    # to be removed, e.g. due to primer efficiency filter) sel: list with indices for
    # each covered template of a primer that should be retained. all others are
    # removed (not considered covered anymore!)
    if (length(filtered.df$Forward) == 0 || length(filtered.df$Reverse) == 0) {
        stop("Need some primers to update cvg data.")
    }
    if (length(mode) == 0) {
        stop("Please supply the 'mode' argument.")
    }
    mode <- match.arg(mode)
    updated.df <- filtered.df
    # all cols relating to cvg info were strsplit is necessary need to be updated
    special.cols <- c("Mismatch_pos_fw", "Mismatch_pos_rev")  # () delimiters
    # only update the currently active constraints in order to prevent error when data are outdated for some reason.
    template.constraints <- active.constraints
    #print(template.constraints)
    # off-binding conditions:
    other.cols <- c("Covered_Seqs", "Binding_Position_Start_fw", "Binding_Position_End_fw", 
            "Binding_Position_Start_rev", "Binding_Position_End_rev", "Binding_Region_Allowed", 
            "Nbr_of_mismatches_fw", "Nbr_of_mismatches_rev", "Relative_Forward_Binding_Position_Start_fw", 
            "Relative_Forward_Binding_Position_End_fw", "Relative_Forward_Binding_Position_Start_rev", 
            "Relative_Forward_Binding_Position_End_rev", "Relative_Reverse_Binding_Position_Start_fw", 
            "Relative_Reverse_Binding_Position_End_fw", "Relative_Reverse_Binding_Position_Start_rev", 
            "Relative_Reverse_Binding_Position_End_rev", "Binding_Region_Allowed_fw", 
            "Binding_Region_Allowed_rev")
    if (mode == "off_target") {
        special.cols <- paste0("Off_", special.cols)
        template.constraints <- paste0("Off_", template.constraints)
        other.cols <- paste0("Off_", other.cols)
    }
    cvg.cols <- c(special.cols, template.constraints, other.cols)
    fw.idx <- which(filtered.df$Forward != "")
    rev.idx <- which(filtered.df$Reverse != "")
    for (i in seq_along(cvg.cols)) {
        if (!cvg.cols[i] %in% colnames(filtered.df)) {
            next
        }
        #message(paste('Updating: ', cvg.cols[i]))
        #print(filtered.df[, cvg.cols[i]])
        fw.col <- grepl("fw", cvg.cols[[i]])
        rev.col <- grepl("rev", cvg.cols[[i]])
        special.col <- cvg.cols[[i]] %in% special.cols  # no differentiation necessary (full mismatch bracket is removed)
        if (special.col) {
            s <- stringr::str_extract_all(filtered.df[, cvg.cols[i]], "(\\([^()]*\\))")
        } else {
            s <- strsplit(filtered.df[, cvg.cols[i]], split = ",")
        }
        s.new <- sapply(seq_along(s), function(x) paste(s[[x]][sel[[x]]], collapse = ","))
        if (fw.col) {
            updated.df[fw.idx, cvg.cols[i]] <- s.new[fw.idx]
        } else if (rev.col) {
            updated.df[rev.idx, cvg.cols[i]] <- s.new[rev.idx]
        } else {
            # col is not direction-specific
            updated.df[, cvg.cols[i]] <- s.new
        }
    }
    # recompute summary statistics: Coverage_Ratio, primer_coverage
    if (mode == "on_target") {
        # TODO: need to update mean annealing DeltaG as well?
        if ("primer_coverage" %in% colnames(updated.df)) {
            updated.df$primer_coverage <- sapply(sel, length)
            updated.df$Coverage_Ratio <- updated.df$primer_coverage/nrow(template.df)
        }
        # update mean efficiency
        if ("primer_efficiency" %in% colnames(updated.df) && "primer_efficiency" %in% active.constraints) {
            updated.df$mean_primer_efficiency <- unlist(lapply(strsplit(updated.df$primer_efficiency, split = ","), function(x) mean(as.numeric(x)))) 
            updated.df$mean_primer_efficiency[is.na(updated.df$mean_primer_efficiency)] <- 0
        }
    }
    # update specificity only for off-target update -> true off-target events are known now!
    if (mode == "off_target") {
        # n.b.: specificity is not perfectly accurate since we store only the best binding event in the target region -> the actual specificity could be higher than the computed one.
        off.cvg <- sapply(strsplit(updated.df$Off_Covered_Seqs, split = ","), function(x) length(unique(x))) # unique off-covered seqs to be on the same scale as the on-target coverage, which is also unique
        TP <- sapply(strsplit(updated.df$Binding_Region_Allowed, split = ","), function(x) length(which(as.logical(x)))) # on-target cvg count
        #print("TP")
        #print(TP) 
        #print("off")
        #print(off.cvg)
        updated.df$primer_specificity <- TP / (TP + off.cvg)
        # if there's not binding events -> 0 specificity
        updated.df$primer_specificity[is.na(updated.df$primer_specificity)] <- 0
    }
    return(updated.df)
}
#' Identification of Mismatch Mutations.
#'
#' Identifies primers that induce mutations due to mismatch binding.
#'
#' Checks for one primer and all covered templates whether any templates are bound
#' with mismatches such that a forbidden mutation is induced. A boolean vector indicating 
#' which binding events induce a forbidden mutation is returned.
#'
#' @param primer.seq Primer sequence string.
#' @param pos.start Binding position of primer (start).
#' @param pos.end Binding position of primer (end).
#' @param template.df Template data frame.
#' @param covered.seqs Identifiers of covered templates.
#' @param ORF.data Reading frame information of templates.
#' @param mode.directionality Directionality of primers.
#' @param mutation.types Character vector of the mutation types to be checked for.
#' @return TRUE if the \code{primer.seq} induces a mutation that is
#' forbidden according to the provided \code{mutation.types}.
#' @keywords internal
check.mutations <- function(primer.seq, pos.start, pos.end, template.df, covered.seqs, 
                            ORF.data, mode.directionality = c("fw", "rev"),
                            mutation.types = c("stop_codon", "substitution")) {

    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' argument.")
    }
    if (length(mutation.types) == 0) {
        return(NULL)
    }
    mutation.types <- match.arg(mutation.types, several.ok = TRUE)
    mode.directionality <- match.arg(mode.directionality)
    if (primer.seq == "") {
        return(NULL)
    }
    template.idx <- match(covered.seqs, template.df$Identifier)
    if (length(template.idx) == 0) {
        # nothing to check
        return(NULL)  # no stop codons found
    }
    ORFs <- ORF.data$Frame
    cur.ORF <- ORFs[template.idx]
    my.lex <- template.df[template.idx, ]
    seqs.exon.nt <- my.lex$Sequence
    primer.seqs.nt <- seqs.exon.nt
    # incorporate the changes from the primer seqs into the sequences
    if (mode.directionality == "fw") {
        substr(primer.seqs.nt, pos.start, pos.end) <- primer.seq
    } else {
        substr(primer.seqs.nt, pos.start, pos.end) <- rev.comp.sequence(primer.seq)
    }
    S <- strsplit(seqs.exon.nt, split = "")
    P <- strsplit(primer.seqs.nt, split = "")
    result <- data.frame(matrix(rep(FALSE, length(template.idx) * length(mutation.types)), 
                        nrow = length(template.idx), ncol = length(mutation.types)))
    colnames(result) <- mutation.types
    for (k in seq_along(S)) {
        frame <- cur.ORF[[k]]
        # compare translated amplicon (primer mismatches) and nt seqs
        seq.aa <- paste(seqinr::translate(S[[k]], frame = frame, sens = "F", ambiguous = TRUE), 
                        collapse = "")
        primer.aa <- paste(seqinr::translate(P[[k]], frame = frame, sens = "F", ambiguous = TRUE), 
                                    collapse = "")
        highlight.aa <- highlight.mismatch(tolower(seq.aa), tolower(primer.aa))
        types <- highlight.aa$type
        #print(seq.aa)
        #print(primer.aa)
        #print(unlist(types))
        idx <- which(sapply(types, function(x) mutation.types %in% x)) 
        result[k, idx] <- TRUE # disallowed mutation type found
    }
    return(result)
}

#' Evaluation of Coverage.
#'
#' Re-evaluates the coverage of primers under exclusion of certain templates.
#'
#' This function requires that \code{primers} was already annotated with primer coverage before.
#'
#' @param primers Primer data frame.
#' @param excluded.seqs Identifiers of templates to be excluded.
#' @param template.df Template data frame
#' @return Primer data frame with updated coverage under the exclusion of \code{excluded.seeqs}.
#' @keywords internal
evaluate.diff.primer.cvg <- function(primers, excluded.seqs, template.df) {
    if (!"primer_coverage" %in% colnames(primers)) {
        stop("Differential primer cvg can only be used when cvg has already been computed!")
    }
    if (length(excluded.seqs) == 0) {
        return(primers)  # nothing to adjust
    }
    # need to update the string with covered seqs in the primer data frame
    s <- strsplit(primers$Covered_Seqs, split = ",")
    i <- NULL
    updated.seqs <- foreach(i = seq_along(s), .combine = rbind) %dopar% {
        seqs <- as.numeric(s[[i]])
        m <- match(seqs, excluded.seqs)
        rm.idx <- which(!is.na(m))
        if (length(rm.idx) != 0) {
            seqs <- seqs[-rm.idx]
        }
        # turn vectors into strings again
        cvg.count <- length(seqs)
        seqs <- paste(seqs, collapse = ",", sep = "")
        r <- data.frame(Covered_Seqs = seqs, primer_coverage = cvg.count, Coverage_Ratio = cvg.count/nrow(template.df), 
            stringsAsFactors = FALSE)
    }
    # update the input data frame
    primers[, colnames(updated.seqs)] <- updated.seqs
    # order again by cvg count
    primers <- primers[order(primers$primer_coverage, decreasing = TRUE), ]
    return(primers)
}

#' Selection of Binding Events
#'
#' Selects primer binding events that are within the allowed binding regions.
#'
#' @param bound.fw Indices of covered templates of a single primer.
#' @param bound.to.allowed.region.fw Corresponding allowed binding regions.
#' @param allowed.other.binding.ratio The ratio of other binding events. If
#' this is different from 0, disallowed binding events will also be reported.
#'
#' @return The indices of the allowed binding events.
#' @keywords internal
select.allowed.binding.events <- function(bound.fw, bound.to.allowed.region.fw, allowed.other.binding.ratio) {
    unique.bound.fw <- unique(bound.fw)
    dup.mapping.fw <- lapply(seq_along(unique.bound.fw), function(x) which(bound.fw == 
        unique.bound.fw[x]))
    sel.idx.fw <- vector("list", length(dup.mapping.fw))
    use.disallowed <- allowed.other.binding.ratio != 0
    for (j in seq_along(dup.mapping.fw)) {
        idx <- dup.mapping.fw[[j]]
        if (!use.disallowed) {
            b.fw <- which(bound.to.allowed.region.fw[idx])
        } else {
            b.fw <- seq_along(bound.to.allowed.region.fw[idx])
        }
        if (length(b.fw) != 0) {
            sel.idx.fw[[j]] <- idx[b.fw]
        } 
    }
    return(sel.idx.fw)
}
#' Selection of Individual Binding Events
#'
#' Selects only binding events of interest.
#'
#' @param fw.binding.filtered IRanges with binding events.
#' @param p.idx Index of binding events to keep.
#' @return An IRanges object containing only the selected binding events.
#' @keywords internal
select.binding.events <- function(fw.binding.filtered, p.idx) {
    # select only the target binding events
    my.binding.fw <- fw.binding.filtered[p.idx]
    # adjust metadata
    metadata(my.binding.fw)$primer_idx <- metadata(my.binding.fw)$primer_idx[p.idx]
    metadata(my.binding.fw)$template_idx <- metadata(my.binding.fw)$template_idx[p.idx]
    metadata(my.binding.fw)$allowed_region <- metadata(my.binding.fw)$allowed_region[p.idx]
    return(my.binding.fw)
}
#' Combination of Binding Events.
#'
#' Appends all binding events.
#'
#' @param my.binding.fw Forward binding events of individual primers.
#' @param my.binding.rev Reverse binding events of individual primers.
#' @param fw.m Forward binding events of paired primers.
#' @param rev.m Reverse binding events of paired primers.
#' @return IRanges of all binding events.
#' @keywords internal
combine.binding.events <- function(my.binding.fw, my.binding.rev, fw.m, rev.m) {
    binding <- c(my.binding.fw, my.binding.rev, fw.m, rev.m)
    # integrate metadata for both directions
    metadata(binding)$primer_idx <- c(metadata(my.binding.fw)$primer_idx, metadata(my.binding.rev)$primer_idx, 
        metadata(fw.m)$primer_idx, metadata(rev.m)$primer_idx)
    metadata(binding)$template_idx <- c(metadata(my.binding.fw)$template_idx, 
        metadata(my.binding.rev)$template_idx, metadata(fw.m)$template_idx, metadata(rev.m)$template_idx)
    metadata(binding)$allowed_region <- c(metadata(my.binding.fw)$allowed_region, 
        metadata(my.binding.rev)$allowed_region, metadata(fw.m)$allowed_region, 
        metadata(rev.m)$allowed_region)
    return(binding)
}
#' Merge of Forward/Reverse Binding Information.
#'
#' Determines binding events of individual and pairs of primers.
#'
#' @param primers The primer data frame.
#' @param fw.binding IRanges object with binding events of fw primers.
#' @param rev.binding IRanges object with binding events of rev primers.
#' @param mode.directionality Primer directionality.
#' @param idx.fw Index of fw primers.
#' @param idx.rev Index of rev primers.
#' @return IRanges with correct binding events.
#' @keywords internal
merge.binding.information <- function(primers, fw.binding.filtered, 
    rev.binding.filtered, mode.directionality = c("fw", "rev", "both"), 
    idx.fw, idx.rev) {
    
    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' argument.")
    }
    mode.directionality <- match.arg(mode.directionality)
    if (mode.directionality == "fw") {
        # only forward primers present
        binding <- fw.binding.filtered
        metadata(binding) <- c(metadata(binding), direction = list(rep("fw", length(binding))))
    } else if (mode.directionality == "rev") {
        # only reverse primers present
        binding <- rev.binding.filtered
        metadata(binding) <- c(metadata(binding), direction = list(rep("rev", length(binding))))
    } else {
        # primers of both directions present
        # differentiate between individual primers and primer pairs:
        merge.idx <- intersect(idx.fw, idx.rev) # primer pairs
        keep.idx <- setdiff(seq_len(nrow(primers)), merge.idx) # individual primers
        # get individual binding events:
        my.binding.fw <- select.binding.events(fw.binding.filtered, 
            metadata(fw.binding.filtered)$primer_idx %in% keep.idx)
        my.binding.rev <- select.binding.events(rev.binding.filtered, 
            metadata(rev.binding.filtered)$primer_idx %in% keep.idx)
        # get pair binding events:
        fw.m <- select.binding.events(fw.binding.filtered, 
            metadata(fw.binding.filtered)$primer_idx %in% merge.idx) # fw bindings of pairs
        rev.m <- select.binding.events(rev.binding.filtered, 
            metadata(rev.binding.filtered)$primer_idx %in% merge.idx) # rev bindings of pairs
        # get idx of fw/rev covered templates
        t.fw <- lapply(seq_len(nrow(primers)), function(x) 
                    metadata(fw.m)$template_idx[metadata(fw.m)$primer_idx == x])
        t.rev <- lapply(seq_len(nrow(primers)), function(x) 
                    metadata(rev.m)$template_idx[metadata(rev.m)$primer_idx == x])
        # determine templates covered by both primers of a pair:
        t <- lapply(seq_len(nrow(primers)), function(x) intersect(t.fw[[x]], t.rev[[x]]))
        ####### cave: multiple binding events are retained here ..
        # select fw binding events of pairs:
        rm <- unlist(lapply(seq_along(t), function(x) which(metadata(fw.m)$template_idx %in% 
            t[[x]] & metadata(fw.m)$primer_idx == x)))
        fw.m <- select.binding.events(fw.m, rm) # fw bindings of pairs
        # select rev binding events of pairs:
        rm <- unlist(lapply(seq_along(t), function(x) which(metadata(rev.m)$template_idx %in% 
            t[[x]] & metadata(rev.m)$primer_idx == x)))
        rev.m <- select.binding.events(rev.m, rm) # rev bindings of pairs
        # integrate pair and individual binding data:
        binding <- combine.binding.events(my.binding.fw, my.binding.rev, fw.m, rev.m)
        # add primer directions to metadata:
        direction <- c(rep("fw", length(my.binding.fw)), rep("rev", length(my.binding.rev)), 
            rep("fw", length(fw.m)), rep("rev", length(rev.m)))
        metadata(binding) <- c(metadata(binding), direction = list(direction))
    }
    return(binding)
}
#' Retrieval of Allowed Binding Indices.
#'
#' Retrieves the indices of allowed binding events in \code{binding}
#' for the primer with index \code{x} and type \code{primer.type}.
#'
#' @param binding IRanges binding information.
#' @param primer.type Direction of primer.
#' @param x Index of primer in the primer data frame.
#' @param allowed.other.binding.ratio The ratio of allowed off-target binding events.
#' @return Indices in binding for primer with index code{x} that are allowed.
#' @keywords internal
get.primer.binding.idx <- function(binding, primer.type = c("fw", "rev", "both"), 
                                x, allowed.other.binding.ratio) {
    # determine individual covered sequences
    if (length(primer.type) == 0) {
        stop("Please supply the 'primer.type' argument.")
    }
    primer.type <- match.arg(primer.type)
	fw.idx <- which(metadata(binding)$direction == "fw" & metadata(binding)$primer_idx == x)
	rev.idx <- which(metadata(binding)$direction == "rev" & metadata(binding)$primer_idx == x)
	# in case of duplicate binding events of one primer into the same template, select allowed ones
	fw.idx <- fw.idx[unlist(select.allowed.binding.events(metadata(binding)$template_idx[fw.idx],
												metadata(binding)$allowed_region[fw.idx], 
                                                allowed.other.binding.ratio))]
	rev.idx <- rev.idx[unlist(select.allowed.binding.events(metadata(binding)$template_idx[rev.idx], 
												metadata(binding)$allowed_region[rev.idx], 
                                                allowed.other.binding.ratio))]
	if (primer.type == "both") {
        # require that both primers of pair bind to allowed regions:
		sel.idx.fw <- fw.idx[which(metadata(binding)$template_idx[fw.idx] %in% metadata(binding)$template_idx[rev.idx])]
		sel.idx.rev <- rev.idx[which(metadata(binding)$template_idx[rev.idx] %in% metadata(binding)$template_idx[fw.idx])]
	} else if (primer.type == "fw") {
		sel.idx.fw <- fw.idx
		sel.idx.rev <- NULL
	} else {
		sel.idx.fw <- NULL
		sel.idx.rev <- rev.idx
	}
    result <- list("fw" = sel.idx.fw, "rev" = sel.idx.rev)
	return(result)
}
#' Retrieval of Relative Binding Positions.
#'
#' Retrieves primer binding position relative to allowed regions of either
#' forward or reverse primers, as specified by \code{direction}.
#'
#' @param allowed Positions where binding is allowed in the templates.
#' @param primer.pos Binding position of primer (absolute).
#' @param direction Direction (either fw/rev).
#' @param covered.seqs.idx Indices of covered templates.
#' @return Numeric of relative binding position to allowed region.
#' @keywords internal
get.relative.binding.pos <- function(allowed, primer.pos, direction, covered.seqs.idx) {
    if (direction == "fw") {
        # relative to fw primer binding region
        rel.pos <- sapply(seq_along(primer.pos), function(x) primer.pos[[x]] - 
            allowed[covered.seqs.idx[x]] - 1)
    } else {
        # relative to rev primer binding region
        rel.pos <- sapply(seq_along(primer.pos), function(x) -(primer.pos[[x]] - 
            allowed[covered.seqs.idx[x]] + 1))
    }
    return(rel.pos)
}
#' Selection of Best (smallest number of mismatches) Binding Event per Template Coverage Event.
#'
#' @param binding Binding information.
#' @param fw.mm.info Info about mismatches.
#' @return A list with entries 'fw' and 'rev' giving the best indices of primers.
#' @keywords internal
select_best_binding <- function(binding, fw.mm.info) {
    result <- vector("list", 2)
    unique.t <- unique(metadata(binding)$template_idx)
    sel.idx.fw.rel <- unlist(lapply(seq_along(unique.t), function(z) {
                    my.t.idx <- which(metadata(binding)$template_idx == unique.t[z])
                    nbr.mm <- fw.mm.info$Nbr[my.t.idx]
                    my.t.idx[which.min(nbr.mm)]}))
    return(sel.idx.fw.rel)
}
#' Evaluation of Primer Coverage.
#'
#' Evaluates the coverage of a set of primers.
#'
#' @param template.df Template data frame.
#' @param primers Primer data frame.
#' @param mode.directionality Primer directionality.
#' @param allowed.mismatches The number of allowed mismatches per binding event.
#' @param allowed.other.binding.ratio Ratio of primers that are allowed to bind to non-allowed regions.
#' If \code{allowed.other.binding.ratio} >0 primers are allowed to bind at any location within the templates.
#' However, a warning is given if the ratio of primers binding to non-target regions exceeds the \code{allowed.other.binding.ratio}.
#' @param allowed.region.definition Definition of the target binding sites used for evaluating the coverage.
#' If \code{allowed.region.definition} is \code{within}, primers have to lie within the allowed binding region.
#' If \code{allowed.region.definition} is \code{any}, primers have to overlap with the allowed binding region.
#' The default is that primers have to bind within the target binding region.
#' @param updateProgress Progress callback function for shiny.
#' @return Primer data frame with information on the covered template sequences.
#' @keywords internal
evaluate.basic.cvg <- function(template.df, primers, mode.directionality = c("fw", "rev", "both"), 
                               allowed.mismatches, allowed.other.binding.ratio, 
                               allowed.region.definition = c("within", "any"), updateProgress = NULL) {

    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' argument.")
    }
    if (length(allowed.region.definition) == 0) {
        stop("Please supply the 'allowed.region.definition' argument.")
    }
    mode.directionality <- match.arg(mode.directionality)
    allowed.region.definition <- match.arg(allowed.region.definition)
    if (length(template.df) == 0 || nrow(template.df) == 0 || length(primers) == 0  || nrow(primers) == 0) {
        return(NULL)
    }
    #message("computing cvg for: ", nrow(primers), " primers.")
    # disregard empty primers
    idx.fw <- which(primers$Forward != "")
    idx.rev <- which(primers$Reverse != "")
    seqs <- Biostrings::DNAStringSet(template.df$Sequence)
    seqs.rc <- Biostrings::DNAStringSet(Biostrings::reverseComplement(seqs))
    # should we consider only as 'on-target', the events in the target region or just consider all events?
    doFilter <- ifelse(allowed.other.binding.ratio == 0, TRUE, FALSE)
    ### fw primers
    fw.binding <- evaluate.cvg(seqs, Biostrings::DNAStringSet(primers$Forward), "fw", allowed.mismatches, 
        updateProgress) # binding events of forward primers
    # determine the allowed binding region:
    allowed.fw <- IRanges(start = template.df$Allowed_Start_fw, end = template.df$Allowed_End_fw)
    # determine allowed binding events for forward primers
    fw.binding.data <- annotate.binding.events(fw.binding, allowed.fw, nrow(primers), allowed.region.definition)
    ### rev primers: analogously
    rev.binding <- evaluate.cvg(seqs.rc, Biostrings::DNAStringSet(primers$Reverse), "rev", allowed.mismatches, 
       updateProgress)
    allowed.rev <- IRanges(start = template.df$Allowed_Start_rev, end = template.df$Allowed_End_rev)
    rev.binding.data <- annotate.binding.events(rev.binding, allowed.rev, nrow(primers), allowed.region.definition)
    if (doFilter) {
        # only retain allowed binding events
        fw.binding.filtered <- fw.binding.data$on_target
        rev.binding.filtered <- rev.binding.data$on_target
    } else {
        # consider all binding events
        fw.binding.filtered <- fw.binding.data$all_binding
        rev.binding.filtered <- rev.binding.data$all_binding
    }
    # require paired primers to bind with both directions, merge results:
    binding <- merge.binding.information(primers, fw.binding.filtered, rev.binding.filtered, mode.directionality, idx.fw, idx.rev)
    # off-target binding:
    off.binding <- merge.binding.information(primers, fw.binding.data$off_target, rev.binding.data$off_target, mode.directionality, idx.fw, idx.rev)
    x <- NULL
    #for (x in seq_along(primers$Identifier)) { # TODO: remove, only for debug 
    on.df <- compute.basic.details(binding, "on_target", template.df, primers, mode.directionality,
                               allowed.mismatches, allowed.other.binding.ratio, 
                               allowed.region.definition, updateProgress)
    off.df <- compute.basic.details(off.binding, "off_target", template.df, primers, mode.directionality,
                               allowed.mismatches, allowed.other.binding.ratio, 
                               allowed.region.definition, updateProgress)
    colnames(off.df) <- paste0("Off_", colnames(off.df))
    cvg.df <- cbind(on.df, off.df)
    ################
    # E: Specificity
    ################
    off.cvg <- unlist(lapply(covered.seqs.to.idx(cvg.df$Off_Covered_Seqs, template.df), function(x) length(unique(x))))
    TP <- unlist(lapply(strsplit(cvg.df$Binding_Region_Allowed, split = ","), function(x) length(which(as.logical(x)))))
    specificity <- TP / (TP + off.cvg)
    specificity[is.na(specificity)] <- 0
    cvg.df$primer_specificity <- specificity
    return(cvg.df)
}
#' Computation of Coverage Details
#'
#' Determines binding properties of primers.
#'
#' @param binding An \code{IRanges} object with primer binding information.
#' @param mode Either \code{on_target} for on-target binding or \code{off_target} for off-target binding.
#' @param template.df Template data frame.
#' @param primers Primer data frame.
#' @param mode.directionality Primer directionality.
#' @param allowed.mismatches The number of allowed mismatches per binding event.
#' @param allowed.other.binding.ratio Ratio of primers that are allowed to bind to non-allowed regions.
#' If \code{allowed.other.binding.ratio} >0 primers are allowed to bind at any location within the templates.
#' However, a warning is given if the ratio of primers binding to non-target regions exceeds the \code{allowed.other.binding.ratio}.
#' @param allowed.region.definition Definition of the target binding sites used for evaluating the coverage.
#' If \code{allowed.region.definition} is \code{within}, primers have to lie within the allowed binding region.
#' If \code{allowed.region.definition} is \code{any}, primers have to overlap with the allowed binding region.
#' The default is that primers have to bind within the target binding region.
#' @param updateProgress Progress callback function for shiny.
#' @return Primer data frame with information on the covered template sequences.
#' @keywords internal
compute.basic.details <- function(binding, mode = c("on_target", "off_target"), template.df, 
                               primers, mode.directionality = c("fw", "rev", "both"), 
                               allowed.mismatches, allowed.other.binding.ratio, 
                               allowed.region.definition = c("within", "any"), updateProgress = NULL) {

    mode <- match.arg(mode)
    x <- NULL
    message("Basic Coverage: ", mode)
    ##########
    # Determine indices of coverage in 'binding' for every primer:
    # motivation: reduce memory consumption of 'foreach' call (doesn't require full 'binding' object)
    #############
    fw.bindings <- vector("list", nrow(primers))
    rev.bindings <- vector("list", nrow(primers))
    meta.cols <- names(binding@metadata)
    for (x in seq_len(nrow(primers))) {
        primer.type <- primers$Direction[x]
        primer.indices <- get.primer.binding.idx(binding, primer.type, x,
                                ifelse(mode == "on_target", allowed.other.binding.ratio, 1))
        cur.bindings <- list(fw = binding[primer.indices$fw, ], rev = binding[primer.indices$rev, ])
        # adjust metadata to selection
        for (y in seq_along(cur.bindings)) {
            cur.binding <- cur.bindings[[y]]
            dir <- names(cur.bindings)[y]
            dir.idx <- primer.indices[[dir]] # select only metadata relating to current direction
            for (col in meta.cols) {
                metadata(cur.bindings[[y]])[[col]] <- metadata(cur.binding)[[col]][dir.idx]
            }
        }
        fw.bindings[[x]] <- cur.bindings[["fw"]]
        rev.bindings[[x]] <- cur.bindings[["rev"]]
    }
    # define iteration vars
    fw.binding <- NULL
    rev.binding <- NULL
    x <- NULL
    cvg.df <- foreach(fw.binding = fw.bindings, rev.binding = rev.bindings, x = iterators::icount(nrow(primers)), .combine = "rbind") %dopar% {
        if (is.function(updateProgress)) {
            detail <- ""
            updateProgress(x/nrow(primers), detail, "set")
        }
        # blueprint output
        p.result <- data.frame(primer_coverage = integer(1), Coverage_Ratio = numeric(1), 
            Binding_Position_Start_fw = character(1), Binding_Position_End_fw = character(1), 
            Binding_Position_Start_rev = character(1), Binding_Position_End_rev = character(1), 
            Binding_Region_Allowed = character(1), Binding_Region_Allowed_fw = character(1),
            Binding_Region_Allowed_rev = character(1), Nbr_of_mismatches_fw = character(1), 
            Nbr_of_mismatches_rev = character(1), Mismatch_pos_fw = character(1), 
            Mismatch_pos_rev = character(1), Covered_Seqs = character(1), 
            # relative binding positions relative to start/end of primers:
            Relative_Forward_Binding_Position_Start_fw = character(1),
            Relative_Forward_Binding_Position_End_fw = character(1), 
            Relative_Forward_Binding_Position_Start_rev = character(1),
            Relative_Forward_Binding_Position_End_rev = character(1), 
            Relative_Reverse_Binding_Position_Start_fw = character(1),
            Relative_Reverse_Binding_Position_End_fw = character(1),
            Relative_Reverse_Binding_Position_Start_rev = character(1),
            Relative_Reverse_Binding_Position_End_rev = character(1),
            stringsAsFactors = FALSE)
        primer.type <- primers$Direction[x]
        #######
        # A: determine covered templates
        ########
		if (primer.type == "both") {
            # templates should be covered by both directions
			covered.seqs.idx <- intersect(metadata(fw.binding)$template_idx,
			                    metadata(rev.binding)$template_idx)
		} else if (primer.type == "fw") {
            covered.seqs.idx <- metadata(fw.binding)$template_idx
		} else {
            covered.seqs.idx <- metadata(rev.binding)$template_idx
		}
        covered.seqs <- template.df$Identifier[covered.seqs.idx]
        ###############
        # Mismatch Info
        ################
        fw.mm.info <- NULL
        rev.mm.info <- NULL
        # on-target:
        subject <- Biostrings::DNAStringSet(template.df$Sequence[covered.seqs.idx])
        w <- Biostrings::width(subject)  # template lengths
        # cumulative index in the concatenation of template strings
        idx <- c(0, cumsum(head(w, length(w) - 1)))
        s <- unlist(subject)
        if (primer.type == "fw" || primer.type == "both") {
            f <- IRanges::Views(s, fw.binding@start + idx, (fw.binding@start + fw.binding@width - 1) + idx)  
            fw.mm.info <- mismatch.info(primers$Forward[x], f) 
        } 
        if (primer.type == "rev" || primer.type == "both") {
            # rev
            f <- IRanges::Views(s, rev.binding@start + idx, 
                        (rev.binding@start + rev.binding@width - 1) + idx)
            rev.mm.info <- mismatch.info(primers$Reverse[x], Biostrings::reverseComplement(f))
        }
        # select only the best binding mode in each template target region: select binding w/ smallest nbr of mismatches
        ##########
        sel.idx.fw <- select_best_binding(fw.binding, fw.mm.info)
        sel.idx.rev <- select_best_binding(rev.binding, rev.mm.info)
        ########
        # update the mismatch info to retain only the selected events
        #########
        # update mismatch info: retain only best mismatch info
        fw.mm.info$Pos <- fw.mm.info$Pos[sel.idx.fw]
        fw.mm.info$Nbr <- fw.mm.info$Nbr[sel.idx.fw]
        rev.mm.info$Pos <- rev.mm.info$Pos[sel.idx.rev]
        rev.mm.info$Nbr <- rev.mm.info$Nbr[sel.idx.rev]
        # update covered seqs: unique for every template
        if (primer.type == "both") {
            # templates should be covered by both directions
			covered.seqs.idx <- intersect(metadata(fw.binding)$template_idx[sel.idx.fw], 
			metadata(rev.binding)$template_idx[sel.idx.rev])
		} else if (primer.type == "fw") {
            covered.seqs.idx <- metadata(fw.binding)$template_idx[sel.idx.fw]
		} else {
            covered.seqs.idx <- metadata(rev.binding)$template_idx[sel.idx.rev]
		}
        covered.seqs <- template.df$Identifier[covered.seqs.idx]
        #####
        # B: determine binding positions in the templates
        ######
        if (length(fw.binding) != 0) {
            # on-target binding: select best binding mode
            primer.exon.start <- fw.binding@start[sel.idx.fw]  # binding pos in the seq
            primer.exon.end <- (fw.binding@start + fw.binding@width - 1)[sel.idx.fw]
       } else {
            primer.exon.start <- NULL
            primer.exon.end <- NULL
       } 
       if (length(rev.binding) != 0) {
            primer.exon.start.rev <- rev.binding@start[sel.idx.rev]
            primer.exon.end.rev <- (rev.binding@start + rev.binding@width - 1)[sel.idx.rev]
        } else {
            primer.exon.start.rev <- NULL
            primer.exon.end.rev <- NULL
        }
        #########
        # C: Relative binding positions
        #########
        # n.b.: not using _initial_here!
        # relative to fw leader:
        rel.primer.start.f <- get.relative.binding.pos(template.df$Allowed_End_fw, primer.exon.start, "fw", covered.seqs.idx)
        rel.primer.end.f <- get.relative.binding.pos(template.df$Allowed_End_fw, primer.exon.end, "fw", covered.seqs.idx)
        rel.primer.start.rev.f <- get.relative.binding.pos(template.df$Allowed_End_fw, primer.exon.start.rev, "fw", covered.seqs.idx)
        rel.primer.end.rev.f <- get.relative.binding.pos(template.df$Allowed_End_fw, primer.exon.end.rev, "fw", covered.seqs.idx)
        # relative to rev leader:
        rel.primer.start.r <- get.relative.binding.pos(template.df$Allowed_Start_rev, primer.exon.start, "rev", covered.seqs.idx)
        rel.primer.end.r <- get.relative.binding.pos(template.df$Allowed_Start_rev, primer.exon.end, "rev", covered.seqs.idx)
        rel.primer.start.rev.r <- get.relative.binding.pos(template.df$Allowed_Start_rev, primer.exon.start.rev, "rev", covered.seqs.idx)
        rel.primer.end.rev.r <- get.relative.binding.pos(template.df$Allowed_Start_rev, primer.exon.end.rev, "rev", covered.seqs.idx)
        ############
        # D: Allowed binding regions check
        ###########
        bound.to.allowed.region.fw <- metadata(fw.binding)$allowed_region[sel.idx.fw]
        bound.to.allowed.region.rev <- metadata(rev.binding)$allowed_region[sel.idx.rev]
        if (primer.type == "both") {
            if (length(bound.to.allowed.region.fw) == 0 || length(bound.to.allowed.region.rev) == 
              0) {
                bound.to.allowed.region <- NULL
            } else {
              bound.to.allowed.region <- bound.to.allowed.region.fw & bound.to.allowed.region.rev
            }
        } else if (primer.type == "fw") {
            bound.to.allowed.region <- bound.to.allowed.region.fw
        } else {
            bound.to.allowed.region <- bound.to.allowed.region.rev
        }
        ############
        # Output results
        #############
        coverage.count <- length(covered.seqs)
        p.result$primer_coverage <- coverage.count
        p.result$Coverage_Ratio <- coverage.count/nrow(template.df)
        # binding position in seqs
        p.result$Binding_Position_Start_fw <- paste(primer.exon.start, collapse = ",")
        p.result$Binding_Position_End_fw <- paste(primer.exon.end, collapse = ",")
        p.result$Binding_Position_Start_rev <- paste(primer.exon.start.rev, collapse = ",")
        p.result$Binding_Position_End_rev <- paste(primer.exon.end.rev, collapse = ",")
        # relative position with regard to fw binding region: fw primer:
        p.result$Relative_Forward_Binding_Position_Start_fw <- paste(rel.primer.start.f, 
            collapse = ",")
        p.result$Relative_Forward_Binding_Position_End_fw <- paste(rel.primer.end.f, 
            collapse = ",")
        # rev primer
        p.result$Relative_Forward_Binding_Position_Start_rev <- paste(rel.primer.start.rev.f, 
            collapse = ",")
        p.result$Relative_Forward_Binding_Position_End_rev <- paste(rel.primer.end.rev.f, 
            collapse = ",")
        # rel to rev binding region:
        p.result$Relative_Reverse_Binding_Position_Start_fw <- paste(rel.primer.start.r, 
            collapse = ",")
        p.result$Relative_Reverse_Binding_Position_End_fw <- paste(rel.primer.end.r, 
            collapse = ",")
        p.result$Relative_Reverse_Binding_Position_Start_rev <- paste(rel.primer.start.rev.r, 
            collapse = ",")
        p.result$Relative_Reverse_Binding_Position_End_rev <- paste(rel.primer.end.rev.r, 
            collapse = ",")
        # allowed binding?
        p.result$Binding_Region_Allowed_fw <- paste(bound.to.allowed.region.fw, 
            collapse = ",")
        p.result$Binding_Region_Allowed_rev <- paste(bound.to.allowed.region.rev, 
            collapse = ",")
        p.result$Binding_Region_Allowed <- paste(bound.to.allowed.region, collapse = ",")
        p.result$Nbr_of_mismatches_fw <- paste(fw.mm.info$Nbr, collapse = ",")
        p.result$Nbr_of_mismatches_rev <- paste(rev.mm.info$Nbr, collapse = ",")
        if (length(fw.mm.info$Pos) != 0) {
            mm.info <- paste(sapply(fw.mm.info$Pos, function(x) paste("(", paste(x, 
              collapse = ","), ")", sep = "")), collapse = ",")
        } else {
            mm.info <- ""
        }
        p.result$Mismatch_pos_fw <- mm.info
        if (length(rev.mm.info$Pos) != 0) {
            mm.info <- paste(sapply(rev.mm.info$Pos, function(x) paste("(", paste(x, 
              collapse = ","), ")", sep = "")), collapse = ",")
        } else {
            mm.info <- ""
        }
        p.result$Mismatch_pos_rev <- mm.info
        p.result$Covered_Seqs <- paste(covered.seqs, collapse = ",")
        #message(paste0(x, "/", nrow(primers))) # iteration status
        p.result
    }
    return(cvg.df)
}
#' Evaluation of Primer Coverage.
#'
#' Evaluates the coverage of a set of primers.
#'
#' @param template.df Template data frame.
#' @param primers Primer data frame.
#' @param mode.directionality Primer directionality.
#' @param settings A \code{DesignSettings} object.
#' @param updateProgress Progress callback function for shiny.
#' @return Primer data frame with information on the covered template sequences.
#' @keywords internal
evaluate.primer.cvg <- function(template.df, primers, mode.directionality = c("fw", "rev", "both"), 
                                settings, updateProgress = NULL) {

    # determine basic coverage: string matching with mismatches
    basic.cvg.df <- evaluate.basic.cvg(template.df, primers, mode.directionality = mode.directionality, 
                                       conOptions(settings)$allowed_mismatches,
                                       conOptions(settings)$allowed_other_binding_ratio,
                                       conOptions(settings)$allowed_region_definition,
                                       updateProgress = updateProgress)
    out.df <- basic.cvg.df
    colnames(out.df) <- paste0("Basic_", colnames(out.df))
    # annotate primers with cvg:
    primers <- update.constraint.values(primers, basic.cvg.df)
    constrained.cvg.df <- evaluate.constrained.cvg(template.df, primers, basic.cvg.df, 
                            mode.directionality, settings, updateProgress = updateProgress)
    out.df <- cbind(out.df, constrained.cvg.df)
    return(out.df)
}
#' Evaluation of Primer Coverage.
#'
#' Evaluates the coverage of a set of primers.
#'
#' @param template.df Template data frame.
#' @param primer.df Primer data frame.
#' @param cvg.df Data frame with basic coverage entries.
#' @param mode.directionality Primer directionality.
#' @param settings A \code{DesignSettings} object.
#' @param updateProgress Progress callback function for shiny.
#' @return Primer data frame with information on the covered template sequences.
#' @keywords internal
evaluate.constrained.cvg <- function(template.df, primer.df, cvg.df, mode.directionality = c("fw", "rev", "both"), 
                                    settings, updateProgress = NULL) {

    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' argument.")
    }
    mode.directionality <- match.arg(mode.directionality)
    if (length(template.df) == 0 || nrow(template.df) == 0 || length(primer.df) == 0  || nrow(primer.df) == 0) {
        return(NULL)
    }
    #t <- proc.time()["elapsed"]
    message("Computing constrained coverage ...")
    #######
    # Selection of binding events
    #######
    cvg.constraints <- cvg_constraints(settings)
    active.constraints <- names(cvg.constraints)
    if (any(grepl("hexamer_coverage", names(conOptions(settings))))) {
        active.constraints <- c(active.constraints, "hexamer_coverage")
    }
    new.df <- check_cvg_constraints(primer.df, template.df, settings, 
                                active.constraints = active.constraints,
                                updateProgress = updateProgress)
    filter.res <- filter.by.constraints(new.df, new.df, cvg.constraints, 
                                      names(cvg.constraints), mode.directionality, template.df)$Filtered
    # only keep the columns relating to coverage:
    keep.cols.basic <- colnames(cvg.df)
    keep.cols.cvg <- colnames(new.df)[unlist(lapply(names(cvg_constraints(settings)), function(x) grep(x, colnames(new.df))))]
    keep.cols <- union(keep.cols.basic, keep.cols.cvg)
    out <- filter.res[, keep.cols[keep.cols %in% colnames(filter.res)]]
    return(out)
}
get_duplex_events <- function(fw.df, annealing.temp, ions) {
    if (length(fw.df) == 0 || nrow(fw.df) == 0) {
        return(NULL)
    }
    #print("get_duplex_Events:")
    #print("fw.df:")
    #print(fw.df)
    # select only unique events
    #print(fw.df)
    dup <- duplicated(fw.df[, c("Primer", "Template")])
    unique.df <- fw.df[!dup,]
    fw.df$Index <- seq_len(nrow(fw.df))
    unique.df$UniqueIndex <- seq_len(nrow(unique.df))
    combi.df <- merge(unique.df, fw.df, by = c("Primer", "Template"))
    # rename identifiers:
    colnames(combi.df) <- gsub("TemplateIdentifier.x", "TemplateIdentifier", colnames(combi.df))
    colnames(combi.df) <- gsub("PrimerIdentifier.x", "PrimerIdentifier", colnames(combi.df))
    combi.df <- combi.df[, !colnames(combi.df) %in% c("TemplateIdentifier.y", "PrimerIdentifier.y")]
    ## restore original order for accessing the unique index:
    combi.df <- combi.df[order(combi.df$Index),]
    #m <- match(combi.df$Template, unique.df$Template)
    duplex.result.fw <- get.dimer.data(unique.df$Primer, unique.df$Template, annealing.temp, ions, no.structures = TRUE) 
    if (length(duplex.result.fw) != 0) {
        duplex.result.fw <- plyr::ddply(duplex.result.fw, c("Idx1"), function(x) arrange(x, substitute(DeltaG))[1, ])
        #duplex.result.fw <- data.frame(PrimerIdentifier = unique.df$PrimerIdentifier, TemplateIdentifier = unique.df$TemplateIdentifier, duplex.result.fw, stringsAsFactors = FALSE)
        duplex.result.fw <- cbind(combi.df, DeltaG = duplex.result.fw[combi.df$UniqueIndex, c("DeltaG")])
    }
    return(duplex.result.fw)
}
#' Determination of the Free Binding Energy.
#'
#' Computest the free energy of annealing between primers and templates. If the \code{mode} is set to "on_target", 
#' the free energies of binding events in the allowed region are computed, while if the \code{mode} is set to
#' "off_target", the free energies of off-target events are computed.
#'
#' @param primer.df A \code{Primers} object.
#' @param template.df A \code{Templates} object.
#' @param annealing.temp The vector of optimal annealing temperatures of the primers.
#' @param settings A \code{DesignSettings} object.
#' @param mode  If the \code{mode} is set to "on_target", 
#' the free energies of binding events in the allowed region are computed, while if the \code{mode} is set to
#' "off_target", the free energies of off-target events are computed.
#' @return A list of lists containing the numeric free energies of the annealing events for every primer.
#' @keywords internal
get.duplex.energies <- function(primer.df, template.df, annealing.temp, settings, mode = c("on_target", "off_target")) {
    ions <- compute.sodium.equivalent.conc(PCR(settings)$Na_concentration, PCR(settings)$Mg_concentration, 
                                           PCR(settings)$K_concentration, PCR(settings)$Tris_concentration)
    if (length(mode) == 0) {
        stop("Please provide the 'mode' argument.")
    }
    mode <- match.arg(mode)
    if (mode == "on_target") {
        # on target binding events
        all.covered.seq.idx <- covered.seqs.to.idx(primer.df$Covered_Seqs, template.df)
        primer.starts.fw <- lapply(strsplit(primer.df$Binding_Position_Start_fw, split = ","), as.numeric)
        primer.ends.fw <- lapply(strsplit(primer.df$Binding_Position_End_fw, split = ","), as.numeric)
        primer.starts.rev <- lapply(strsplit(primer.df$Binding_Position_Start_rev, split = ","), as.numeric)
        primer.ends.rev <- lapply(strsplit(primer.df$Binding_Position_End_rev, split = ","), as.numeric)
    } else {
        # off target binding events
        all.covered.seq.idx <- covered.seqs.to.idx(primer.df$Off_Covered_Seqs, template.df)
        primer.starts.fw <- lapply(strsplit(primer.df$Off_Binding_Position_Start_fw, split = ","), as.numeric)
        primer.ends.fw <- lapply(strsplit(primer.df$Off_Binding_Position_End_fw, split = ","), as.numeric)
        primer.starts.rev <- lapply(strsplit(primer.df$Off_Binding_Position_Start_rev, split = ","), as.numeric)
        primer.ends.rev <- lapply(strsplit(primer.df$Off_Binding_Position_End_rev, split = ","), as.numeric)
    }
    all.covered.templates <- unlist(lapply(all.covered.seq.idx, function(x) template.df$Sequence[x]))
    subject <- Biostrings::DNAStringSet(all.covered.templates)
    w <- Biostrings::width(subject)  # template lengths
    # cumulative index in the concatenation of template strings
    idx <- c(0, cumsum(head(w, length(w) - 1)))
    s <- unlist(subject)
    template.indices <- c(0, cumsum(sapply(all.covered.seq.idx, length)))
    # don't parallelize: dimerization is already parallelized!
    fw.p <- vector("list", nrow(primer.df))
    fw.t <- vector("list", nrow(primer.df))
    rev.p <- vector("list", nrow(primer.df))
    rev.t <- vector("list", nrow(primer.df))
    for (i in seq_len(nrow(primer.df))) {
        primer.type <- primer.df$Direction[i]
        if (template.indices[[i]] == template.indices[[i+1]]) {
            # no coverage
            next
        }
        cur.template.indices <- seq(template.indices[[i]] + 1, template.indices[[i+1]])
        cur.idx.adjustment <- idx[cur.template.indices]
        f.fw <- NULL
        f.rev <- NULL
        if (primer.type == "fw") {
            # fw
            f.fw <- IRanges::Views(s, primer.starts.fw[[i]] + cur.idx.adjustment, primer.ends.fw[[i]] + cur.idx.adjustment)  
        } else if (primer.type == "rev") {
            # rev
            f.rev <- IRanges::Views(s, primer.starts.rev[[i]] + cur.idx.adjustment, primer.ends.rev[[i]] + cur.idx.adjustment)  
        } else {
            # both
            f.fw <- IRanges::Views(s, primer.starts.fw[[i]] + cur.idx.adjustment, primer.ends.fw[[i]] + cur.idx.adjustment)  
            f.rev <- IRanges::Views(s, primer.starts.rev[[i]] + cur.idx.adjustment, primer.ends.rev[[i]] + cur.idx.adjustment)  
        }
        # forward deltaG:
        fw.templates <- rev.comp.sequence(tolower(as.character(f.fw)))
        fw.primers <- rep(primer.df$Forward[i], length(fw.templates))
        # reverse deltaG:
        rev.templates <- tolower(as.character(f.rev))
        rev.primers <- rep(primer.df$Reverse[i], length(rev.templates))
        # write-out:
        fw.p[[i]] <- fw.primers
        fw.t[[i]] <- fw.templates
        rev.p[[i]] <- rev.primers
        rev.t[[i]] <- rev.templates
    }
    fw.ids <- unlist(lapply(seq_along(fw.p), function(x) rep(as.character(primer.df$Identifier[x]), length(fw.p[[x]]))))
    t.ids <- as.character(template.df$Identifier[unlist(all.covered.seq.idx)])
    if (length(fw.ids) != 0) {
        fw.t.ids <- t.ids
    } else {
        fw.t.ids <- NULL
    }
    fw.df <- data.frame(PrimerIdentifier = fw.ids, TemplateIdentifier = fw.t.ids, Primer = unlist(fw.p), Template = unlist(fw.t), stringsAsFactors = FALSE)
    duplex.data.fw <- get_duplex_events(fw.df, annealing.temp[i], ions)
    rev.ids <- unlist(lapply(seq_along(rev.p), function(x) rep(as.character(primer.df$Identifier[x]), length(rev.p[[x]]))))
    if (length(rev.ids) != 0) {
        rev.t.ids <- t.ids
    } else {
        rev.t.ids <- NULL
    }
    rev.df <- data.frame(PrimerIdentifier = rev.ids, TemplateIdentifier = rev.t.ids, Primer = unlist(rev.p), Template = unlist(rev.t), stringsAsFactors = FALSE)
    duplex.data.rev <- get_duplex_events(rev.df, annealing.temp[i], ions)
    out <- vector("list", nrow(primer.df))
    for (i in seq_len(nrow(primer.df))) {
        if (length(duplex.data.fw) != 0) {
            cur.fw.data <- duplex.data.fw[duplex.data.fw$PrimerIdentifier == primer.df$Identifier[i], ]
        } else {
            cur.fw.data <- NULL
        }
        if (length(duplex.data.rev) != 0) {
            cur.rev.data <- duplex.data.rev[duplex.data.rev$PrimerIdentifier == primer.df$Identifier[i], ]
        } else {
            cur.rev.data <- NULL
        }
        if (length(cur.fw.data) != 0 && length(cur.rev.data) != 0) {
            out[[i]] <- sapply(seq_len(nrow(cur.fw.data)), function(x) min(cur.fw.data$DeltaG[x], cur.rev.data$DeltaG))
        } else if (length(cur.fw.data) != 0) {
            out[[i]] <- cur.fw.data$DeltaG
        } else if (length(cur.rev.data) != 0) {
            out[[i]] <- cur.rev.data$DeltaG
        }
    }
    # control:
    # ctrl.df <- cbind(fw.df, "DeltaG" = unlist(out))
    # print(ctrl.df)
    return(out)
}
#' Identification of Mutations Induced by Mismatch Binding Events.
#'
#' Identifies whether mutations are induced by mismatch binding events.
#'
#' Checks for one primer and all covered templates whether any templates are bound
#' with mismatches such that mismatches are induced. A numeric vector indicating 
#' which binding events induce a forbidden mismatch according to 
#' \code{mutation.types} is returned such that \code{1} indicates forbidden events
#' and \code{0} allowed events.
#'
#' @param primer.df A \code{Primers} object.
#' @param template.df A \code{Template} object.
#' @param mutation.types Character vector of the mutation types to be checked for.
#' @return A list containing data frames where an entry of 1 is present if the \code{primer.seq} induces a mutation that is
#' forbidden according to the provided \code{mutation.types}, otherwise 0.
#' @keywords internal
mismatch.mutation.check <- function(primer.df, template.df, mutation.types = c("stop_codon", "substitution")) {
    mutation.types <- match.arg(mutation.types, several.ok = TRUE) 
    all.covered.seqs <- strsplit(primer.df$Covered_Seqs, split = ",")
    ORF.data <- get.ORFs(template.df)
    primer.starts.fw <- lapply(strsplit(primer.df$Binding_Position_Start_fw, split = ","), as.numeric)
    primer.ends.fw <- lapply(strsplit(primer.df$Binding_Position_End_fw, split = ","), as.numeric)
    primer.starts.rev <- lapply(strsplit(primer.df$Binding_Position_Start_rev, split = ","), as.numeric)
    primer.ends.rev <- lapply(strsplit(primer.df$Binding_Position_End_rev, split = ","), as.numeric)
    x <- NULL
    #for (x in seq_i
    stop.list <- foreach(x = seq_len(nrow(primer.df)), .combine = "c") %dopar% {
        covered.seqs <- all.covered.seqs[[x]]
        stop.check.fw <- check.mutations(primer.df$Forward[x], primer.starts.fw[[x]], primer.ends.fw[[x]],
                                           template.df, covered.seqs, ORF.data, "fw", mutation.types)
        stop.check.rev <- check.mutations(primer.df$Reverse[x], primer.starts.rev[[x]], primer.ends.rev[[x]], 
                                            template.df, covered.seqs, ORF.data, "rev", mutation.types)
        primer.type <- primer.df$Direction[x]
        if (primer.type == "both") {
            # any stop codon in a pair -> exclusion
            if (length(stop.check.fw) != 0 && length(stop.check.rev) != 0) {
                stop.check <- stop.check.fw | stop.check.rev
            } else {
                stop.check <- NULL
           }
        } else if (primer.type == "fw") {
            stop.check <- stop.check.fw
        } else if (primer.type == "rev") {
            stop.check <- stop.check.rev
        }
        # output 1 for stop codon and 0 for no stop codon
        stop.check <- 1 * stop.check # numeric
        list(stop.check)
    }
    return(stop.list)
}
get_feature_matrix <- function(primer.df, template.df, mode = c("on_target", "off_target")) {
    # transform to a long-data frame of primer-coverage events
    mode <- match.arg(mode)
    if (!"primer_coverage" %in% colnames(primer.df)) {
        stop("Primer coverage required for model.")
    }
    if (mode == "on_target") {
        if (!"annealing_DeltaG" %in% colnames(primer.df)) {
            stop("Annealing free energy required for model.")
        }
    } else {
        if (!"off_annealing_DeltaG" %in% colnames(primer.df)) {
            stop("Annealing free energy required for model.")
        }
    }
    full.df <- prepare_mm_plot(primer.df, template.df, mode = mode)
    # consider all identified events (this is constrained here, as we use the 'basic coverage entries' (not filtered yet!)
    full.df <- full.df[full.df$Coverage_Type == "constrained",]
    # select unique binding event for all mismatch contacts:
    # the minimal deltaG is selected to ensure that we consider the
    # disambiguated primer with the least nbr of mismatches
    df <- plyr::ddply(full.df, c("Primer", "Template", "Group"), plyr::summarize,
                        annealing_DeltaG = min(substitute(annealing_DeltaG)),
                        Number_of_mismatches = min(substitute(Number_of_mismatches)),
                        Position_3terminus = min(substitute(Position_3terminus)),
                        Position_3terminusLocal = max(substitute(Position_3terminusLocal)))
    if (length(df) != 0 && nrow(df) != 0 && (any(is.na(df$annealing_DeltaG) || is.na(df$Position_3terminusLocal)))) {
        warning("Some values of the coverage model feature matrix were NA; this shouldn't happen!")
    }
    return(df)
}
#' Prediction of Primer Coveragee.
#'
#' Predicts primer coverage using a logistic regression model.
#' Converts coverage probabilities to expected false positive rate
#' for a given probability.
#'
#' @param primer.df A \code{Primers} data frame.
#' @param template.df A \code{Templates} data frame.
#' @param settings A \code{DesignSettings} object.
#' @param mode Whether on-target or off-target events shall be considered.
#' @return The predictions for primer coverage
#' @keywords internal
predict_coverage <- function(primer.df, template.df, settings, mode = c("on_target", "off_target"), updateProgress = NULL) {
    ##################
    # sysdata used:
        # CVG_MODEL (the predictive model)
        # FPR_TABLE (probability to FPR conversion for predictions)
    #########
    if (length(mode) == 0) {
        mode <- "on_target"
    } else {
        mode <- match.arg(mode)
    }
    if (!is(primer.df, "Primers")) {
        stop("Please input a primer data frame.")
    }
    if (!is(template.df, "Templates")) {
        stop("Please provide a valid template data frame.")
    }
    mode.directionality <- get.analysis.mode(primer.df)
    # compute annealingDeltaG:
    cvg_constraints(settings)$annealing_DeltaG <- c("max" = 0)
    if (mode == "on_target") {
        p.df <- compute.constraints(primer.df, mode.directionality, 
                    template.df, settings, 
                    active.constraints = "annealing_DeltaG")
    } else {
        p.df <- compute.constraints(primer.df, mode.directionality, 
                    template.df, settings, 
                    active.constraints = "off_annealing_DeltaG")
    }
    primer.df <- update.constraint.values(primer.df, p.df)
    pred.matrix <- get_feature_matrix(primer.df, template.df, mode = mode)
    if (length(pred.matrix) != 0 && nrow(pred.matrix) != 0) {
        pred <- try(stats::predict(CVG_MODEL, newdata = pred.matrix, type = "response"))
    } else {
        # if 'newdata' is NULL, the training data would be predicted ..
        pred <- NULL
    }
    if (class(pred) == "try-error") {
        stop("Could not predict coverage using the logistic model. Maybe some columns are missing in the data frame? Provided columns were: ", colnames(pred.matrix), ". Dimension of matrix:", dim(pred.matrix))
    }
    #print(pred)
    # transform coverage probabilities to FPR
    idx <- unlist(lapply(pred, function(x) which.min(abs(x - FPR_TABLE$Cutoff))))
    fpr <- FPR_TABLE$FPR[idx]
    # check function
    #summary <- data.frame("Coverage_Probability" = round(pred, 3), "FPR" = round(fpr, 3))
    #print(summary)
    # need to output values as a list, one for each primer; one entry for every coverage event
    if (mode == "on_target") {
        idx <- covered.seqs.to.idx(primer.df$Covered_Seqs, template.df)
    } else {
        idx <- covered.seqs.to.idx(primer.df$Off_Covered_Seqs, template.df)
    }
    pred.idx <- match(pred.matrix$Template, template.df$ID) # idx of templates in the template data frame 
    # annotate each primer with its FPR vector
    if (length(pred.matrix) != 0) {
        out <- lapply(seq_len(nrow(primer.df)), function(x) {
            p.idx <- which(as.character(pred.matrix$Primer) == as.character(primer.df$ID[x]))
            t.idx.pred <- pred.idx[p.idx] # template idx
            t.idx.p <- idx[[x]]
            m <- match(t.idx.pred, t.idx.p)
            if (any(is.na(m))) {
                stop("Fatal error in predict_coverage: could not find coverage event! Maybe the template IDs were non-unique?!")
            }
            fpr[p.idx][m]
        })
    } else {
        out <- vector("list", nrow(primer.df))
    }
    return(out)
}
#' Evaluation of Coverage.
#'
#' Evaluates primer coverage.
#' 
#' @param template.seqs Template sequences as a \code{DNAStringSet}.
#' @param primers Primer sequences as a \code{DNAStringSet}.
#' @param mode.directionality Directionality of primres
#' @param allowed.mismatches Allowed number of mismatches between a primer and a template.
#' @param updateProgress Progress function for shiny
#' @return IRanges object with primer coverage information.
#' @keywords internal
evaluate.cvg <- function(template.seqs, primers, 
                    mode.directionality = c("fw", "rev"), 
                    allowed.mismatches, updateProgress = NULL) {

    # todo: implement progress function?
    if (length(mode.directionality) == 0) {
        stop("Please provide the 'mode.directionality' arg.")
    }
    mode.directionality <- match.arg(mode.directionality)
    if (length(template.seqs) == 0 || length(primers) == 0) {
        return(NULL)
    }
    if (!is(primers, "DNAStringSet")) {
        stop("primers must be a DNAStringSet.")
    }
    if (!is(primers, "DNAStringSet")) {
        stop("seqs must be a DNAStringSet.")
    }
    # consider ambiguous templates up to a cutoff of degeneration
    seqs <- my.disambiguate(template.seqs)
    # determine mapping from disambiguated seqs to input seqs
    l.p <- unlist(lapply(primers, length))
    l.s <- unlist(lapply(seqs, length))
    l.p <- rep(seq_along(primers), l.p)
    l.s <- rep(seq_along(seqs), l.s)
    seqs <- unlist(seqs)
    seqDB <- unlist(seqs)  # one continuous string of all template seqs
    # determine mapping of pos in seqDB to template id
    w <- Biostrings::width(seqs)  # template lengths
    idx <- c(0, cumsum(w))
    w <- c(0, w)
    template.mapping <- IRanges(start = idx[1:(length(idx) - 1)] + 1, width = w[2:length(w)])
    names(template.mapping) <- l.s
    f <- IRanges()  # forward start/end positions in the templates
    fp <- integer()  # fp: forward primer indices, one index for every template hit
    non.empty.idx <- which(Biostrings::width(primers) != 0)
    i <- NULL
    f <- foreach(i = seq_along(non.empty.idx), .combine = c) %dopar% {

        idx <- non.empty.idx[i]
        # matchPattern returns all matches below the max mismatch cutoff! need to select the best match later on
        temp <- Biostrings::matchPattern(primers[[idx]], seqDB, max.mismatch = allowed.mismatches, 
                                        with.indels = FALSE, fixed = FALSE) # fixed FALSE to allow ambigs
        temp <- as(temp, "IRanges")
        # add additional info
        fp <- rep(idx, length(temp))
        # store the primer idx as metadata using GRanges
        if (length(temp) != 0) {
            temp <- GenomicRanges::GRanges(seqnames = mode.directionality, ranges = temp, primer_idx = fp)  # changed from RangedData (deprecated)
        } else {
            temp <- GenomicRanges::GRanges()
        }
    }
    if (length(f) == 0) {
        return(IRanges())
    }
    fp <- as.numeric(as.character(f$primer_idx))  # primer indices
    f <- IRanges::Views(seqDB, start(f), end(f))  # extract the sub-sequence of the primer matches
    # try to map from 'f' Iranges to template.mapping to find which template seqs
    # were bound
    overlaps <- IRanges::findOverlaps(f, template.mapping)  # f is query, template.mapping is the subject
    fw.table <- f[as.matrix(overlaps)[, 1]] # as.matrix from IRanges ..
    fw.idx <- fp[as.matrix(overlaps)[, 1]]
    fw.template.table <- template.mapping[as.matrix(overlaps)[, 2]]  # names: identifier of seq
    # extract regions of binding: subtract fw.table from fw.template.table
    fw.tab <- as(fw.table, "IRanges")
    fw.starts <- fw.tab@start - fw.template.table@start + 1
    fw.ends <- fw.starts + fw.tab@width - 1
    template.ends <- fw.template.table@width
    r <- fw.starts > 0 & fw.ends <= template.ends
    df <- data.frame(Primer = fw.idx, Seq = fw.template.table@NAMES)
    if (mode.directionality == "rev") {
        binding.regions <- IRanges(start = fw.template.table@width[r] - fw.starts[r] + 
                                    1 - fw.tab@width[r] + 1, width = fw.tab@width[r])
    } else {
        binding.regions <- IRanges(start = fw.starts[r], width = fw.tab@width[r])
    }
    metadata(binding.regions) <- list(primer_idx = fw.idx[r], template_idx = as.numeric(fw.template.table@NAMES[r]))
    return(binding.regions)
}
#' Information about Mismatches.
#' 
#' Computes information about mismatch binding events.
#'
#' @param primer Primer character vector.
#' @param seqs Template binding sequences of primers as a \code{XStringsView} object.
#' @return List with positions and number of mismatches of the \code{primer} in the 
#' \code{seqs}. The list contains the field \code{mm.pos} containing a list
#' with the positions of the mismatches and the field \code{Nbr} containing
#' a numeric vector with the number of mismatches per template binding event.
#' @keywords internal
mismatch.info <- function(primer, seqs) {
    if (length(seqs) == 0) {
        return(NULL)
    }
    p <- my.disambiguate(DNAStringSet(primer))[[1]]
    s <- my.disambiguate(DNAStringSet(seqs))
    # select template with minimal mismatches
    As <- lapply(p, function(p.seq) {
        modes <- lapply(s, function(s.seq) Biostrings::compareStrings(rep(tolower(as.character(p.seq)), length(s.seq)), tolower(as.character(s.seq))))
        mm.pos <- lapply(modes, function(x) lapply(strsplit(x, split = ""), function(y) which(y == "?")))
        mm.count <- lapply(mm.pos, function(x) sapply(x, length))
        idx <- unlist(lapply(mm.count, function(x) which.min(x)))
        res <- lapply(seq_along(idx), function(x) mm.pos[[x]][[idx[x]]])
    })
    # select primer with minimal mismathches
    mm.pos <- lapply(seq_along(seqs), function(x) {
        nbr <- sapply(As, function(p.seqs) length(p.seqs[[x]]))
        sel <- which.min(nbr)
        As[[sel]][[x]]
    })
    mm.count <- sapply(mm.pos, length) # number of mismatches
    result <- list(Pos = mm.pos, Nbr = mm.count)
    return(result)
}
#' Annotation of Primer Binding Events.
#'
#' Annotates whether primer binding events are in the allowed binding region or not. 
#'
#' @param fw.binding IRanges with coverage information.
#' @param allowed.range IRanges of the allowed binding ranges in the templates.
#' @param nbr.primers Number of primers to consider.
#' @param allowed.region.definition Definition of the allowed binding region
#' @return IRanges with annotations of (preliminary) specificity and allowed binding. 
#' The field \code{all_binding} contains all binding regions, \code{on_target} contains
#' all events in the target region, and \code{off_target} contains all off-target
#' binding events.
#' @keywords internal
annotate.binding.events <- function(fw.binding, allowed.range, nbr.primers,
                                 allowed.region.definition = c("within", "any")) {
    if (length(allowed.region.definition) == 0) {
        stop("Please provide the 'allowed.region.definition' argument.")
    }
    allowed.region.definition <- match.arg(allowed.region.definition)
    out <- list("on_target" = IRanges::IRanges(), "all_binding" = IRanges::IRanges(), "off_target" = IRanges::IRanges())
    if (length(fw.binding) == 0) {
        return(out)  # nothing to filter/annotate
    }
    t.idx <- metadata(fw.binding)$template_idx
    p.idx <- metadata(fw.binding)$primer_idx
    t.u <- unique(t.idx)
    # determine allowed binding events
    p.u <- unique(p.idx)
    # determine the allowed range corresponding to each binding event in 'fw.binding' -> replicate allowed range
    cur.bindings <- vector("list", nbr.primers)
    template.indices <- vector("list", nbr.primers)
    for (i in seq_along(p.u)) {
        idx <- which(p.idx == p.u[i])
        template.indices[[i]] <- t.idx[idx]  # current set of templates to be checked for primer binding
        cur.bindings[[i]] <- fw.binding[idx]
    }
    # unify all bindings
    all.bindings <- do.call(c, cur.bindings)
    all.allowed <- allowed.range[unlist(template.indices)]
    # either checks for 'any' overlaps or 'within' allowed region overlaps of primers with template regions
    allowed <- IRanges::overlapsAny(all.bindings, all.allowed, type = allowed.region.definition)
    metadata(fw.binding) <- c(metadata(fw.binding), allowed_region = list(allowed))
    # select on-target binding events
    allowed.binding <- fw.binding[allowed]
    # update existing metadata
    for (col in names(allowed.binding@metadata)) {
        metadata(allowed.binding)[[col]] <- metadata(fw.binding)[[col]][allowed]
    }
    # select off-target binding events: 
    disallowed.binding <- fw.binding[!allowed]
    # update existing metadata
    for (col in names(disallowed.binding@metadata)) {
        metadata(disallowed.binding)[[col]] <- metadata(fw.binding)[[col]][!allowed]
    }
    out <- list("on_target" = allowed.binding, "all_binding" = fw.binding, "off_target" = disallowed.binding)
    return(out)
}
