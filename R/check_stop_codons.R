#######
# Functions for checking templates for introduced stop codons due to mismatches
#######

#' Highlight mismatches
#' 
#' Collects information on the mutations present in the input and highlights the
#' mutations in the sequence.
#' 
#' @param seq character vector of the original sequence
#' @param mm.seq  character vector of the mutated sequence 
#' 
#' @return A list highlighting the mutations and additional information 
#' (mutation type, number, etc.)
#' @keywords internal
highlight.mismatch <- function(seq, mm.seq) {
    # search for introduced mutations / stop codons
    mutation.types <- vector("list", length(seq))
    nbr.mismatches <- rep(NA, length(seq))
    pos <- vector("list", length(seq))  # mismatch posis
    for (i in seq_along(seq)) {
        s <- unlist(strsplit(seq[i], split = ""))
        mm <- unlist(strsplit(mm.seq[i], split = ""))
        if (length(s) != length(mm)) {
            warning("Sequence lengths do not agree. Mutation analysis might be incorrect.")
        }
        idx <- which(s != mm)
        # find identical stop codons:
        stop.idx <- intersect(which(s == "*"), which(mm == "*"))
        idx <- c(idx, stop.idx)
        pos[[i]] <- idx
        if (length(idx) > 0) {
            s[idx] <- toupper(s[idx])
            mm[idx] <- toupper(mm[idx])
        }
        seq[i] <- paste(s, collapse = "")
        mm.seq[i] <- paste(mm, collapse = "")
        # determine mutation type
        if (length(idx) == 0) {
            # no AA mutation found -> silent genomic mutation
            mutation.types[[i]] <- "silent"
        } else {
            stop.idx <- which(mm[idx] == "*")
            substitution.idx <- setdiff(idx, stop.idx)
            if (length(stop.idx) != 0) {
                # stop codon in AA seq
                mutation.types[[i]] <- "stop_codon"
            } 
            if (length(substitution.idx) != 0) {
                # DNA substitution 
                # else in order to call stop codons as stop codons
                # and not as substitutions (subs are not as bad!)
                mutation.types[[i]] <- c(mutation.types[[i]], "substitution")
            }
        }
        nbr.mismatches[i] <- length(idx)
    }
    result <- list(seq = seq, mm = mm.seq, type = mutation.types, nbr = nbr.mismatches, pos = pos)
    #print(result)
    return(result)
}
#' Identification of ORFs.
#' 
#' Given a template data frame, identify the exon reading frames in the sequences.
#' 
#' @param template.df template data frame.
#' 
#' @return Returns a data frame containing the shift of the ORF (either 0,1, or 2) 
#' for every sequence, as well as a comment in case of problems.
#' @keywords internal
get.ORFs <- function(template.df) {
    nSeqs <- seq_len(nrow(template.df))
    S <- strsplit(template.df$Sequence, split = "")
    ORFs <- rep(NA, length(nSeqs))
    comments <- rep("", length(ORFs))
    frames <- rep(FALSE, length(nSeqs))
    # consider all possible reading frames:
    frames <- 0:2
    k <- NULL
    reading.frames <- foreach (k = seq_along(frames), .combine = "c") %dopar% {
        f <- frames[k] # current reading frame
        seq.aa <- unlist(lapply(S, function(x) paste(seqinr::translate(x, frame = f, sens = "F", 
                                    NAstring = "X", ambiguous = TRUE), collapse = "")))
        sel <- !grepl("\\*", seq.aa) # true if there's no stop codon in the translated seq
        list(data.frame("No_Stop_Codon" = sel, "Translation" = seq.aa, stringsAsFactors = FALSE))
    }
    frame.df <- do.call(cbind, lapply(reading.frames, function(x) x$No_Stop_Codon)) # one column for every ORF
    non.stop.counts <- apply(frame.df, 1, function(x) length(which(x)))
    # if non.stop.count is 1 -> only possible frame found
    # if non.stop.count is 0 -> non-functional seq -> determine right framee
    # if non.stop.count is greater than 1 -> multiple possible frames -> determine right frame
    translations <- rep(NA, nrow(template.df))
    idx <- which(non.stop.counts == 1)
    ORFs[idx] <- apply(frame.df[idx, , drop = FALSE], 1, which) - 1
    # determine the correct ORF for sequences with stop codons / multiple non-stop codon translations
    na.idx <- which(non.stop.counts > 1)
    na.idx.na <- which(non.stop.counts == 0) 
    if (length(na.idx) != 0) {
        comments[na.idx] <- "Multiple frames without stop codons found. Translation frame was chosen based on similarity."
    }
    if (length(na.idx.na) != 0) {
        comments[na.idx.na] <- "All template translation frames contained stop codons. Check your input sequence. Translation frame was chosen based on similarity."
        warning("Templates with IDs: ", paste(template.df$ID[na.idx.na], collapse = ","), ". ", paste0(unique(comments[na.idx.na], collapse = ",")))
    }
    missing.idx <- c(na.idx, na.idx.na) # idx of missing ORF seqs
    # use the existing translations for determining the ORFs for the missing seqs
    seq.aa <- sapply(seq_along(idx), function(x) reading.frames[[ORFs[idx[x]] + 1]]$Translation[idx[x]])
    translations[idx] <- seq.aa
    if (length(seq.aa) != 0) {
        for (k in seq_along(missing.idx)) {
            cur.idx <- missing.idx[k]
            cur.translations <- sapply(reading.frames, function(x) x$Translation[cur.idx])
            # determine best fit to any other correctly translated sequence
            dists <- unlist(lapply(cur.translations, function(test.seq) min(as.vector(stringdist::stringdistmatrix(test.seq, seq.aa, method = "lv")))))
            #print(dists) # for stop codon seq, it may not be clear which frame to choose ..
            min.idx <- which.min(dists)
            ORFs[cur.idx] <- min.idx - 1
            # store best translations:
            if (k %in% na.idx) {
                # multiple frames without stop codon found
                translations[cur.idx] <- cur.translations[min.idx]
            } else {
                # no translation without stop found: give all translations
                translations[cur.idx] <- paste0(cur.translations, collapse = "|")
            }
        }
    }
    # there should be no remaining ORFs with NA or -1, except if seq.aa is empty
    idx <- which(is.na(ORFs))
    if (length(idx) != 0) {
        ORFs[idx] <- 0  # just set first ORF arbitrarily
        warning(paste0(template.df$ID[idx], ". Could not determine the ORFs for all sequences."))
        comments[idx] <- paste0(comments[idx], "ORF unknown: no available seqs for comparison.")
    }
    result <- data.frame(ID = template.df$ID, Frame = ORFs, Comment = comments, Sequence = template.df$Sequence, Translation = translations)
    return(result)
}

#' Adjust ORF position
#' 
#' Adjusts the reading frame according to the position at which we consider a
#' subsequence. 
#' 
#' @param ORFs the reading frames (0,1,2).
#' @param seq.start the position where a sequence is extracted .
#' 
#' @return The adjusted reading frames for the given start positions.
#' @keywords internal
adjust.ORF.start <- function(ORFs, seq.start) {
    frame.offset <- (seq.start - 1 - ORFs)%%3
    new.frame <- (3 - frame.offset)%%3
    return(new.frame)
}

#' Split a sequence 
#'
#' Splits a sequence at a specified positions
#'
#' @param target The target string.
#' @param index The position for the split.
#'
#' @return List with splitted strings
#' @keywords internal
split_str_by_index <- function(target, index) {
    index <- sort(index)
    substr(rep(target, length(index) + 1), start = c(1, index), stop = c(index - 
        1, nchar(target)))
}

#' Interleave strings

#' Combines the input vectors in an interleaved fashion. 
#' @param v1 Input string.
#' @param v2 Input string.
#'
#' @return The interleaved combination of \code{v1} and \code{v2}.
#' @keywords internal
interleave <- function(v1, v2) {
    ord1 <- 2 * (1:length(v1)) - 1
    ord2 <- 2 * (1:length(v2))
    result <- c(v1, v2)[order(c(ord1, ord2))]
    return(result)
}
#' String Insertion.
#'
#' Inserts a string into another string at the speficied position.
#'
#' @param target The string to be modified.
#' @param insert The string to be inserted.
#' @param index The position where the insertion should take place.
#'
#' @return A string where \code{insert} is inserted into \code{target} at position \code{index}.
#' @keywords internal
insert_str <- function(target, insert, index) {
      insert <- insert[order(index)]
    index <- sort(index)
    return(paste(interleave(split_str_by_index(target, index), insert), collapse = ""))
}
#' Format mismatches
#' 
#' Formats a sequence for highlighting mismatches in an alignment.
#'
#' @param seq The input sequence.
#' @param pos  The mismatch positions to be formatted.
#' @param format.type Vector of giving the style (bold/italics) for each pos.
#'
#' @return The input sequence with highlighted mismatch positions.
#' @keywords internal
format.seq.ali <- function(seq, pos, format.type) {
       s <- toupper(seq)
    if (length(pos) == 0) {
        # no mismatches
        return(s)
    }
    if (length(format.type) != length(pos)) {
        stop("style length should be the same as format pos length")
    }
    style.start <- rep(NA, length(format.type))
    style.end <- rep(NA, length(format.type))
    idx <- format.type == "bold"
    style.start[idx] <- "<b>"
    style.end[idx] <- "</b>"
    idx <- format.type == "italics"
    style.start[idx] <- "<i>"
    style.end[idx] <- "</i>"
    ins <- unlist(lapply(seq_along(style.start), function(x) c(style.start[x], style.end[x])))
    positions <- unlist(lapply(pos, function(x) c(x, x + 1)))
    s <- insert_str(s, ins, positions)
    return(s)
}
#' Mismatch overview table
#'
#' Computes a table summarizing all of the mismatches caused by the primers in the
#' input data frame.
#' @param primer.df Primer data frame.
#' @param template.df Template data.
#' @param mode.directionality Direction of primers.
#'
#' @return: A data frame summarizing all mismatches of the input primers with the input templates.
#' @keywords internal
compute.mismatch.table <- function(primer.df, template.df, mode.directionality = c("fw", "rev")) {

    if (length(primer.df) == 0 || nrow(primer.df) == 0) {
        return(NULL)
    }
    mismatch.table <- NULL
    seq.idx <- covered.seqs.to.idx(primer.df$Covered_Seqs, template.df)  # seq idx isn't correct for j= 2
    ORF.data <- get.ORFs(template.df)
    ORFs <- ORF.data$Frame
    ali.style <- "<div style = 'font-family: Courier,courier;'>"
    for (j in seq_len(nrow(primer.df))) {
        # message(j)
        primer <- primer.df[j, ]
        primer.seq <- NA
        if (mode.directionality == "fw") {
            primer.seq <- primer$Forward
        } else {
            primer.seq <- primer$Reverse
        }
        if (primer.seq == "" || primer$primer_coverage == 0) {
            # nothing to check :-)
            next
        }
        idx <- seq.idx[[j]]
        ids <- template.df$ID[idx]
        cur.ORF <- ORFs[idx]
        comments <- ORF.data$Comment[idx]
        # retrieve primer-specific ORFs for ORFs that are NA (stop codon might be in
        # another region than the primer-binding region)
        if (mode.directionality == "fw") {
            start.pos <- as.numeric(unlist(strsplit(primer$Binding_Position_Start_fw, 
                split = ",")))
            end.pos <- as.numeric(unlist(strsplit(primer$Binding_Position_End_fw, 
                split = ",")))
            extract.pos <- rep(1, length(start.pos))  # extract string from first position
        } else {
            start.pos <- as.numeric(unlist(strsplit(primer$Binding_Position_Start_rev, 
                split = ",")))
            end.pos <- as.numeric(unlist(strsplit(primer$Binding_Position_End_rev, 
                split = ",")))
            extract.pos <- start.pos  # extract string from start.pos of binding
        }
        seqs.exon.nt <- template.df$Sequence[idx]
        primer.seqs.nt <- seqs.exon.nt
        # incorporate the changes from the primer seqs into the sequences
        is.indel <- ifelse(nchar(substr(primer.seqs.nt, start.pos, end.pos)) != nchar(primer.seq), 
            TRUE, FALSE)  # with indels: this is different than 
        if (mode.directionality == "fw") {
            substr(primer.seqs.nt, start.pos, end.pos) <- primer.seq
        } else {
            substr(primer.seqs.nt, start.pos, end.pos) <- rev.comp.sequence(primer.seq)
        }
        sel <- which(primer.seqs.nt != seqs.exon.nt)  # mismatches were introduced
        if (length(sel) == 0) {
            # nothing to analyze
            next
        }
        primer.seqs.nt <- primer.seqs.nt[sel]
        seqs.exon.nt <- seqs.exon.nt[sel]
        cur.ORF <- cur.ORF[sel]
        ids <- ids[sel]
        primer.id <- primer$ID
        comments <- comments[sel]
        # see if translation has an asterisk
        S <- strsplit(gsub("-", "", seqs.exon.nt), split = "")
        P <- strsplit(gsub("-", "", primer.seqs.nt), split = "")
        seq.aa <- lapply(seq_along(cur.ORF), function(x) seqinr::translate(S[[x]], frame = cur.ORF[x], 
            sens = "F", NAstring = "X", ambiguous = TRUE))
        seq.aa.string <- sapply(seq.aa, function(x) paste(x, collapse = ""))
        idx <- which(sapply(seq.aa, function(x) "*" %in% x))
        cur.ORF <- as.list(cur.ORF)
        if (length(idx) != 0) {
            comments[idx] <- paste(comments[idx], "Best primer binding region translation (", 
                paste(seq.aa.string[idx], sep = ""), ") had a stop codon. Displaying all three translation frames.", 
                sep = "")
            for (x in idx) {
                cur.ORF[[x]] <- 0:2
            }
        }
        highlight.nt <- highlight.mismatch(tolower(seqs.exon.nt), tolower(primer.seqs.nt))
        for (k in seq_along(S)) {
            primer.pos <- regexpr(primer.seq, tolower(highlight.nt$mm[k]))  # pos of primer binding
            primer.pos <- as.numeric(primer.pos):(as.numeric(primer.pos) + 
                          attr(primer.pos, "match.length") - 1)
            style.info <- rep("bold", length(highlight.nt$pos[[k]]))
            form.pos <- highlight.nt$pos[[k]]
            # form.pos <- c(form.pos, primer.pos) style.info <- c(style.info, rep('italics',
            # length(primer.pos)))
            mm.format <- format.seq.ali(highlight.nt$mm[k], form.pos, style.info)
            NT.alignment <- paste(ali.style, format.seq.ali(highlight.nt$seq[k], 
                highlight.nt$pos[[k]], rep("bold", length(highlight.nt$pos[[k]]))), 
                "<br>", mm.format, "</div>", sep = "")
            frames <- cur.ORF[[k]]
            seq.aa <- unlist(lapply(seq_along(frames), function(x) paste(seqinr::translate(S[[k]], 
                frame = frames[x], sens = "F", NAstring = "X", ambiguous = TRUE), 
                collapse = "")))
            primer.aa <- unlist(lapply(seq_along(frames), function(x) paste(seqinr::translate(P[[k]], 
                frame = frames[x], sens = "F", NAstring = "X", ambiguous = TRUE), 
                collapse = "")))
            highlight.aa <- highlight.mismatch(tolower(seq.aa), tolower(primer.aa))
            if (length(highlight.aa$seq) == 1) {
                AA.alignment <- paste(ali.style, format.seq.ali(highlight.aa$seq[1], 
                                highlight.aa$pos[[1]], rep("bold", length(highlight.aa$pos[[1]]))), 
                                "<br>", format.seq.ali(highlight.aa$mm[[1]], highlight.aa$pos[[1]], 
                                rep("bold", length(highlight.aa$pos[[1]]))), "</div>", sep = "")
            } else {
                AA.alignment <- NA
            }
            types <- highlight.aa$type
            out.types <- sapply(types, function(x) paste(x, collapse = ";")) 
            output <- data.frame(Primer = primer.id, Primer_Identifier = as.character(primer$Identifier), 
                                Template = ids[k], Primer_Seq = primer$Forward, 
                                Mutation_Type = out.types,
                                NT_nbr_mm = highlight.nt$nbr[k], AA_nbr_mm = paste(highlight.aa$nbr, 
                                collapse = ";", sep = ""), Seq_NT = highlight.nt$seq[k], 
                                Primer_NT = highlight.nt$mm[k], 
                                Seq_AA = paste(highlight.aa$seq, collapse = ";", sep = ""), 
                                Primer_AA = paste(highlight.aa$mm, collapse = ";", sep = ""), 
                                Comment = comments[k], Alignment_NT = NT.alignment, 
                                Alignment_AA = AA.alignment, stringsAsFactors = FALSE)
            mismatch.table <- rbind(mismatch.table, output)
        }
    }
    return(mismatch.table)
}

