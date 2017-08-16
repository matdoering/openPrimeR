########
# Functions for formatted output.
###########

#' View the Input Primers.
#'
#' Creates a formatted primers table.
#'
#' @param primer.df A \code{Primers} object.
#' @param mode.directionality The direction of the primers.
#' @param for.shiny Whether output is intended for Shiny.
#' @return A formatted primer table.
#' @keywords internal
view.input.primers <- function(primer.df, mode.directionality, for.shiny = TRUE) {
    if (nrow(primer.df) == 0 || length(primer.df) == 0) {
        return(NULL)
    }
    view.df <- asS3(primer.df) # modifying columns here -> cannot guarantee object validity
    # remove all columns where all values are missing
    excl.idx <- which(unlist(lapply(seq_len(ncol(view.df)), function(x) 
                all(view.df[,x] == "" | is.na(view.df[,x])))))
    # additionally remove columns depending on run mode
    excl <- NULL
    if (mode.directionality == "fw") {
        excl <- grep("_rev", colnames(view.df))
    } else if (mode.directionality == "rev") {
        excl <- grep("_fw", colnames(view.df))
    }
    excl.idx <- unique(c(excl.idx, excl))
    if (length(excl.idx) != 0) {
        view.df <- view.df[, -excl.idx]
    }
    excl.col <- c("Identifier", "Run")
    view.df <- exclude.cols(excl.col, view.df)
    view.df <- modify.col.rep(view.df, for.shiny)
    return(view.df)
}
#' View the Evaluated Primers.
#'
#' Creates a formatted primers table.
#'
#' @param primer.df A \code{Primers} object.
#' @param template.df A \code{Templates} object.
#' @param mode.directionality The direction of the primers.
#' @param view.cvg.individual Whether information on individual coverage events should be retrained.
#' @param for.shiny Whether the table is intended for Shiny (HTML) or not.
#' @return A formatted primer table.
#' @keywords internal
view.cvg.primers <- function(primer.df, template.df, mode.directionality, 
                            view.cvg.individual = c("active", "inactive"), 
                            for.shiny = TRUE) {
    if (length(primer.df) == 0) {
        return(NULL)
    } else if (nrow(primer.df) == 0) {
        return(primer.df)
    }
    if (length(view.cvg.individual) == 0) {
        stop("Please supply the 'view.cvg.individual' arg.")
    }
    view.cvg.individual <- match.arg(view.cvg.individual)
    view.df <- asS3(primer.df)
    # convert binding positions to range
    view.df$Binding_Position_Start_fw <- pos.to.range(view.df$Binding_Position_Start_fw, 
        view.df$Binding_Position_End_fw)
    view.df$Binding_Position_Start_rev <- pos.to.range(view.df$Binding_Position_Start_rev, 
        view.df$Binding_Position_End_rev)
    view.df$Relative_Forward_Binding_Position_Start_fw <- pos.to.range(view.df$Relative_Forward_Binding_Position_Start_fw, 
        view.df$Relative_Forward_Binding_Position_End_fw)
    view.df$Relative_Forward_Binding_Position_Start_rev <- pos.to.range(view.df$Relative_Forward_Binding_Position_Start_rev, 
        view.df$Relative_Forward_Binding_Position_End_rev)
    view.df$Relative_Reverse_Binding_Position_Start_fw <- pos.to.range(view.df$Relative_Reverse_Binding_Position_Start_fw, 
        view.df$Relative_Reverse_Binding_Position_End_fw)
    view.df$Relative_Reverse_Binding_Position_Start_rev <- pos.to.range(view.df$Relative_Reverse_Binding_Position_End_rev, 
        view.df$Relative_Reverse_Binding_Position_Start_rev)
    # reorder columns:
    col.order <- c("ID", "Forward", "Reverse", "Coverage_Ratio",
                   "Nbr_of_mismatches_fw", "Nbr_of_mismatches_rev", 
                   "Relative_Forward_Binding_Position_Start_fw", 
                   "Relative_Reverse_Binding_Position_Start_rev",
                  "melting_temp", "melting_temp_diff", "Self_Dimer_DeltaG",
                  "Cross_Dimer_DeltaG", "Structure_deltaG", 
                  "primer_specificity", "mean_primer_efficiency",
                  "gc_ratio_fw", "gc_ratio_rev", "gc_clamp_fw",
                  "gc_clamp_rev", "no_repeats_fw", "no_repeats_rev", 
                  "no_runs_fw", "no_runs_rev",
                  "Covered_Seqs",  # primer coverage not used here
                  "Self_Dimer_Structure", "Cross_Dimer_Structure",
                  "Structure_fw", "Structure_rev", 
                  "Degeneracy_fw", "Degeneracy_rev",
                  "primer_length_fw", "primer_length_rev", 
                  "Direction")
    cur.order <- col.order[col.order %in% colnames(view.df)]
    other.cols <- setdiff(colnames(view.df), cur.order)  # don't care about these cols in ordering
    view.df <- view.df[, c(cur.order, other.cols)]
    # don't order!
    #view.df <- view.df[order(view.df$primer_coverage, decreasing = TRUE), ]
    # remove unimportant cols from table
    view.df <- view.df[, colnames(view.df) %in% col.order]
    # percent format some columns:
    percent.format.cols <- c("Coverage_Ratio", "gc_ratio_fw", "gc_ratio_rev",
                             "mean_primer_efficiency", "primer_specificity")
    for (i in seq_along(percent.format.cols)) {
        col <- percent.format.cols[i]
        if (col %in% colnames(view.df)) {
            view.df[, col] <- paste0(round(view.df[, col], 3) * 100, "%")
        }
    }
        # convert covered template seqs to group representation
    cvd <- covered.seqs.to.ID.string(as.character(view.df$Covered_Seqs), template.df)
    if (length(unique(template.df$Group)) >= 2 && view.cvg.individual == "inactive") {
        # show gene group instead of identifiers
        idx <- covered.seqs.to.idx(as.character(view.df$Covered_Seqs), template.df)
        cvd <- string.list.format(sapply(seq_along(idx), function(x) paste(template.df[idx[[x]], 
            "Group"], collapse = ",")))
    }
    view.df$Covered_Seqs <- unlist(cvd)
    format.cols <- c("Binding_Position_Start_fw", "Binding_Position_Start_rev", 
        "Relative_Forward_Binding_Position_Start_fw", "Relative_Reverse_Binding_Position_Start_rev", "Nbr_of_mismatches_fw", "Nbr_of_mismatches_rev", 
        "Binding_Region_Allowed", "Binding_Region_Allowed_fw", "Binding_Region_Allowed_rev", 
        "T_EVAL_primer_efficiency", "constraints_passed_T")
    if (view.cvg.individual == "inactive") {
        for (i in seq_along(format.cols)) {
            if (format.cols[i] %in% colnames(view.df)) {
                if (format.cols[i] %in% c("Nbr_of_mismatches_fw", "Nbr_of_mismatches_rev")) {
                    view.df[, format.cols[i]] <- string.list.format(view.df[, format.cols[i]], order.mode = "value")
                } else {
                    view.df[, format.cols[i]] <- string.list.format(view.df[, format.cols[i]])
                }
            }
        }
        mismatch.cols <- c("Mismatch_pos_fw", "Mismatch_pos_rev")
        for (j in seq_along(mismatch.cols)) {
            col <- mismatch.cols[j]
            if (col %in% colnames(view.df)) {
                view.df[, col] <- string.list.format(gsub("[()]", "", 
                    gsub(",()", "", 
                        gsub("(),", "", view.df[, col], fixed = TRUE), 
                    fixed = TRUE)))
            }
        }
    }
    # html representation of dimer structures:
    structure.cols <- c("Self_Dimer_Structure", "Cross_Dimer_Structure")
    for (j in seq_along(structure.cols)) {
        col <- structure.cols[j]
        if (col %in% colnames(view.df)) {
            view.df[, col] <- html.format.structure(as.character(view.df[, col]))
        }
    }
    #####
    # if we don't have primers of both directions, merge fw and rev annotations
    both.mode <- any(view.df$Forward != "" & view.df$Reverse != "")
    fw.idx <- which(view.df$Forward != "")
    rev.idx <- which(view.df$Reverse != "")
    if (!both.mode) {
        # merge constraint columns to single column
        cols <- c("gc_ratio", "gc_clamp", "Nbr_of_mismatches",
                  "no_repeats", "no_runs", "primer_length",
                  "Structure", "Sequence", "Relative_Binding_Position",
                  "Degeneracy")
        special.cols <- list("Sequence" = c("Forward", "Reverse"),
                            "Relative_Binding_Position" = 
                            c("Relative_Forward_Binding_Position_Start_fw",
                            "Relative_Reverse_Binding_Position_Start_rev"))
        for (i in seq_along(cols)) {
            col <- cols[i]
            if (col %in% names(special.cols)) {
                col.dir <- special.cols[[col]]
            } else {
                col.dir <- paste0(col, c("_fw", "_rev"))
            }
            if (all(col.dir %in% colnames(view.df))) {
                merge.col <- unlist(lapply(seq_len(nrow(view.df)), function(x) 
                                if (x %in% fw.idx) {
                                    view.df[x, col.dir[1]] 
                                } else {
                                    view.df[x, col.dir[2]]
                                }))
                view.df[, col.dir[1]] <- merge.col
                # rename col
                colnames(view.df)[colnames(view.df) == col.dir[1]] <- col
                # remove old columns
                view.df <- view.df[,-which(colnames(view.df) == col.dir[2])]
            }
        }
    }
    ########
    excl.col <- c("Binding_Position_End_fw", "Binding_Position_End_rev", "Relative_Forward_Binding_Position_Start_rev", 
        "Relative_Reverse_Binding_Position_Start_fw", "Relative_Forward_Binding_Position_End_fw", 
        "Relative_Reverse_Binding_Position_End_rev", "Relative_Reverse_Binding_Position_End_fw", 
        "Relative_Forward_Binding_Position_End_rev")
    view.df <- exclude.cols(excl.col, view.df)

    # round all numeric columns:
    nums <- vapply(view.df, is.numeric, FUN.VALUE = logical(1))
    view.df[,nums] <- round(view.df[,nums], 2)
    ############
    # rename columns for frontend (format/units)
    name.change.cols <- c("melting_temp", "melting_temp_diff", "Self_Dimer_DeltaG", 
                          "Cross_Dimer_DeltaG", "Structure_deltaG",
                          "Coverage_Ratio", 
                          "gc_ratio", "gc_ratio_fw", "gc_ratio_rev", 
                          "gc_clamp", "gc_clamp_fw", "gc_clamp_rev", 
                          "Nbr_of_mismatches", 
                          "Nbr_of_mismatches_fw", "Nbr_of_mismatches_rev", 
                          "no_repeats", "no_repeats_fw", "no_Repeats_rev",
                          "no_runs", "no_runs_fw", "no_runs_rev", 
                          "Structure", "Structure_fw", "Structure_rev",
                          "mean_primer_efficiency", "primer_specificity",
                          "Binding_Position_Start",
                          "Binding_Position_Start_fw", 
                          "Binding_Position_Start_rev", 
                          "Relative_Binding_Position",
                          "Relative_Forward_Binding_Position_Start_fw", 
                          "Relative_Reverse_Binding_Position_Start_rev")
    if (for.shiny) {
        new.names <- c("T<sub>m</sub> [&deg;C]", "&Delta;T<sub>m</sub> [&deg;C]",
                    "Self Dimer &Delta;G[<sup>kcal</sup>&frasl;<sub>mol</sub>]", 
                    "Cross Dimer &Delta;G[<sup>kcal</sup>&frasl;<sub>mol</sub>]",
                    "Structure &Delta;G[<sup>kcal</sup>&frasl;<sub>mol</sub>]",
                    "Coverage", "GC_ratio", "GC_ratio_fw", "GC_ratio_rev", 
                    "GC_clamp", "GC_clamp_fw", "GC_clamp_rev", 
                    "Mismatches", "Mismatches_fw", "Mismatches_rev", 
                    "Repeat count", "Repeat count_fw", "Repeat count_rev", 
                    "Run count", "Run count_fw", "Run count_rev", 
                    "Secondary structure", "Secondary structure_fw",
                    "Secondary structure_rev", 
                    "efficiency", "specificity",
                    "Template_Binding_Range", 
                    "Template_Binding_Range_fw", "Template_Binding_Range_rev", 
                    "Relative_Binding_Range",
                    "Relative_fw_Binding_Range", "Relative_rev_Binding_Range")
    } else {
        new.names <- c("Tm [C]", "Delta Tm [C]",
                    "Self Dimer", "Cross Dimer", "Structure",
                    "Coverage", "GC_ratio", "GC_ratio_fw", "GC_ratio_rev", 
                    "GC_clamp", "GC_clamp_fw", "GC_clamp_rev", 
                    "Mismatches", "Mismatches_fw", "Mismatches_rev", 
                    "Repeat count", "Repeat count_fw", "Repeat count_rev", 
                    "Run count", "Run count_fw", "Run count_rev", 
                    "Secondary structure", "Secondary structure_fw",
                    "Secondary structure_rev", 
                    "efficiency", "specificity",
                    "Template_Binding_Range", 
                    "Template_Binding_Range_fw", "Template_Binding_Range_rev", 
                    "Relative_Binding_Range",
                    "Relative_fw_Binding_Range", "Relative_rev_Binding_Range")
    }
    idx <- match(colnames(view.df), name.change.cols)
    colnames(view.df)[!is.na(idx)] <- new.names[idx[!is.na(idx)]]
    view.df <- view.input.primers(view.df, mode.directionality, for.shiny)
    return(view.df)
}

#' View the Evaluated Primers.
#'
#' Creates a formatted primers table.
#'
#' @param primer.df A \code{Primers} object.
#' @param template.df A \code{Templates} object.
#' @return A formatted primer table.
#' @keywords internal
view.primers <- function(primer.df, template.df) {
    mode.directionality <- get.analysis.mode(primer.df)
    # pair primers for a more compact representation
    primer.df <- pair_primers(primer.df, template.df)
    out.df <- view.cvg.primers(primer.df, template.df, 
                            mode.directionality, view.cvg.individual = "inactive", 
                            for.shiny = FALSE)
    # select a column subset:
    out.df <- out.df # not here ...
    # format relative binding range: select the majority binding region
    range.cols <- c("Relative Binding Range", "Relative Binding Range (fw)", "Relative Binding Range (rev)")
    for (i in seq_along(range.cols)) {
        col <- range.cols[i]
        if (col %in% colnames(out.df)) {
            out.df[, col] <- unlist(lapply(strsplit(out.df[, col], split = ","),
          function(x) {
            str <- ifelse(length(x) != 0, x[[1]], "")
            out <- ""
            if (str != "") {
                # the relative binding range: e.g. 0 to 30
                conv <- stringr::str_extract_all(str,"\\(?[0-9.%-]+\\)?")[[1]]
                if (as.numeric(conv[2]) > 0) { # end of binding
                    conv[2] <- paste0("+", conv[2])
                }
                if (as.numeric(conv[1]) > 0) { # start of binding
                    conv[1] <- paste0("+", conv[1])
                }
                out <- paste0(paste0(conv[1:2], collapse = " to "), " ", conv[3], collapse = " ")
            }
            return(out)
            }))
        }
    }
    # abbreviate IDs:
    out.df$ID <- abbreviate(out.df$ID, getOption("openPrimeR.plot_abbrev"))
    return(out.df)
}
#' Format a Sequence for LateX output.
#'
#' Formats a sequence for LateX report output.
#'
#' @param seqs Character vector of sequences.
#' @return Formatted sequences.
#' @keywords internal
format.seqs.tex <- function(seqs) {
    seqs <- seqs
    # italicize ambiguous bases 
    ambig.bases <- tolower(setdiff(names(Biostrings::IUPAC_CODE_MAP), Biostrings::DNA_BASES))
    replacement <- paste0("\\\\textit{", ambig.bases, "}")
    names(replacement) <- ambig.bases
    seqs <- stringr::str_replace_all(seqs, replacement)
    return(seqs)
}
#' View the Evaluated Primers in the Report.
#'
#' Creates a formatted primers table for the report PDF.
#'
#' @param primer.df A \code{Primers} object.
#' @param template.df A \code{Templates} object.
#' @return A formatted primer table.
#' @keywords internal
view.primers.report <- function(primer.df, template.df) {
    out.df <- view.primers(primer.df, template.df)
    # determine nbr of columns
    if ("Sequence" %in% colnames(out.df)) {
        # single direction: more space available
        out.df <- out.df[,1:6]
    } else {
        # pairs of primers: less space available
        out.df <- out.df[,1:4]
    }
    col.names <- colnames(out.df) 
    # change mismatch representation to IQR representation:
    # TODO error for rev primers ...
    mm.values <- sapply(seq_len(nrow(primer.df)), function(x) paste0(primer.df[x, "Nbr_of_mismatches_fw"], 
                                                                primer.df[x, "Nbr_of_mismatches_rev"], collapse = ","))
    #print(mm.values)
    #print("MM values:")
    mm.rep <- string.to.IQR(mm.values)
    # TODO: integrate value columns when merging primers in 'view.primers()'
    if (length(mm.rep) == nrow(out.df)) {
        # TODO: for both, the size of the primer data frame is reduced
        # and we haven't defined the values of the merged columns!!!
        # -> number of mismatches is NA at the moment, also, we don't show more than 4 columns in the report anyway.
        out.df$Mismatches <- mm.rep
    }
    # tidy up the column names that are displayed only (the first couple of columns):
    col.names[col.names == "Tm [C]"] <- "$\\text{T}_\\text{m} [^{\\circ}\\text{C}]$"
    col.names[col.names == "Delta Tm [C]"] <- "$\\Delta\\text{T}_\\text{m} [^{\\circ}\\text{C}]$"
    col.names[col.names == "Relative Binding Range"] <- "Position"
    col.names[col.names == "Relative Binding Range (fw)"] <- "Position (rev)"
    col.names[col.names == "Relative Binding Range (rev)"] <- "Position (fw)"
    col.names[col.names == "Self Dimer"] <- "$\\text{Self Dimer} \\Delta\\text{G}$"
    col.names[col.names == "Cross Dimer"] <- "$\\text{Cross Dimer } \\Delta\\text{G}$"
    # sequence can either be in 'Sequence' or in 'Forward' and 'Reverse' for 'both', right? TODO. test report for 'both'.
    out.df[col.names == "Sequence"] <- format.seqs.tex(out.df$Sequence)
    out.df[col.names == "Forward"] <- format.seqs.tex(out.df$Forward)
    out.df[col.names == "Reverse"] <- format.seqs.tex(out.df$Reverse)
    colnames(out.df) <- col.names
    return(out.df)
}

#' Exclusion of Columns
#'
#' Removes columns from a data frame.
#'
#' @param excl.col Names of columns in \code{template.df} to be removed.
#' @param template.df Data frame for which columns in \code{excl.col} should be removed.
#' @return \code{template.df} with removed columns as specified in \code{excl.col}.
#' @keywords internal
exclude.cols <- function(excl.col, template.df) {
    m <- match(excl.col, colnames(template.df))
    if (any(!is.na(m))) {
        view.df <- template.df[, -m[!is.na(m)]]
    } else {
        view.df <- template.df
    }
    return(view.df)
}

#' Modification of Column Names.
#'
#' Modifies column names for frontend output.
#' 
#' @param template.df The data frame whose column names are to be modified.
#' @param for.shiny Whether formatting should be for shiny.
#' @return \code{template.df} with modified column names.
#' @keywords internal
modify.col.rep <- function(template.df, for.shiny = TRUE) {
    # identicate direction of constraint with brackets:
    colnames(template.df) <- gsub("_fw", " (fw)", colnames(template.df))
    colnames(template.df) <- gsub("_rev", " (rev)", colnames(template.df))
    # replace underscores in columns:
    colnames(template.df) <- gsub("_", " ", colnames(template.df))
    # upper case first letter
    colnames(template.df) <- Hmisc::capitalize(colnames(template.df))
    return(template.df)
}
#' Format a String List.
#'
#' Formats a string list, summarizing values with percentages.
#'
#' @param values The string list to format.
#' @param order.mode How the result should be ordered.
#' For "percentage", the strings are ordered by their percentages, while
#' for "value", the strings are ordered by their values.
#' @return A formatted string with percentage annotations.
#' @keywords internal
string.list.format <- function(values, order.mode = c("percentage", "value")) {
    # values in comma separated string will be formatted
    order.mode <- match.arg(order.mode)
    v <- strsplit(values, split = ",")
    n.idx <- which(sapply(v, function(x) length(x)) == 0)
    # for every list element determine its frequency distribution
    dist <- lapply(v, function(x) table(x))
    if (order.mode == "percentage") {
        # order by largest percentage
        dist <- lapply(dist, function(x) x[order(x, decreasing = TRUE)])
    } else if (order.mode == "value") {
        values <- unique(unlist(lapply(dist, function(x) as.numeric(names(x)))))
        values <- as.character(values[order(values)])
        dist <- lapply(dist, function(x) {
                        x <- x[values]
                        x[is.na(x)] <- 0
                        names(x) <- values
                        return(x)
                       })
    } else {
        stop("Unknown order mode.")
    }
    f <- unlist(lapply(dist, function(x) paste(names(x), 
        " (", round((x/sum(x, na.rm = TRUE)), 
        2) * 100, "%)", collapse = ",", sep = "")))
    if (length(n.idx) != 0) {
        f[n.idx] <- ""
    }
    return(f)
}

