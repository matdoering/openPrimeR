#' Overview of Primer Set Properties.
#'
#' Creates an overview of the properties of multiple primer sets by providing
#' the inter-quartile range of primer properties in bracket notation.
#'
#' @param template.data List with \code{Template} objects corresponding
#' to \code{primer.data}. 
#' @param primer.data List with evaluated \code{Primers} objects whose properties
#' are to be summarized.
#' @param sample.name Either a single identifier or identifiers for every
#' \code{Templates} object in \code{template.data}.
#' By default, \code{sample.name} is \code{NULL} such that the
#' \code{Run} annotations in the \code{Templates} objects provided
#' by \code{template.data} are used.
#' @return A data frame summarizing the properties of each primer set.
#' @export
#' @examples
#' data(Comparison)
#' tab <- get_comparison_table(template.data[1:3], primer.data[1:3], "IGH")
get_comparison_table <- function(template.data, primer.data, sample.name = NULL) {
    if (length(template.data) != 0) {
        if (length(sample.name) == 0) {
            sample.name <- sapply(template.data, function(x) unique(x$Run))
        } else if (length(sample.name) == 1) {
            sample.name <- rep(sample.name, length(template.data))
        } else if (length(sample.name) != length(template.data)) {
            stop("Incorrect number of sample names provided.")
        }
        if (length(template.data) != length(primer.data)) {
            stop("Please ensure that you provide the same number of templates as primers.")
        }
    }
    template.classes <- sapply(template.data, function(x) class(x))
    if (any(template.classes != "Templates")) {
        stop("Input templates should be a list with objects ",
             "of class 'Templates'.")

    }
    primer.classes <- sapply(primer.data, function(x) class(x))
    if (any(primer.classes != "Primers")) {
        stop("Input primers should be a list with objects ",
             "of class 'Primers'.")
    }
    all.results <- vector("list", length(primer.data))
    for (j in seq_along(primer.data)) {
        primer.df <- primer.data[[j]]
        lex.seq <- template.data[[j]]
        set.name <- unique(primer.df$Run)
        mode.directionality <- get.analysis.mode(primer.df)
        set.stats <- primer.set.parameter.stats(primer.df, mode.directionality, lex.seq)
        if (length(set.stats) != 0) {
            cur.result <- data.frame(Primers = set.name, Templates = sample.name[j], stringsAsFactors = FALSE)
            set.stats <- cbind(cur.result, set.stats, row.names = NULL, stringsAsFactors = FALSE)
        }
        all.results[[j]] <- set.stats
    }
    sample.result <- do.call(plyr::rbind.fill, all.results)
    # sort sample result by cvg
    if ("Coverage" %in% colnames(sample.result)) {
        hits <- regexpr("[0-9\\.]+%", sample.result$Coverage)
        cvg <- regmatches(sample.result$Coverage, hits)
        cvg <- substr(cvg, 1, nchar(cvg) - 1)
        o <- order(as.numeric(cvg), decreasing = TRUE)
    }
    if (length(sample.result) != 0) {
        # nice labels for columns without underscores:
        colnames(sample.result) <- gsub("_", " ", colnames(sample.result))
    }
    return(sample.result)
}
#' Preparation of Comparison Plot for Evaluation.
#'
#' @param primer.data List with objects of class \code{Primers}. Each 
#' list entry corresponds to a single primer set.
#' @param constraint.settings List with settings for each constraint.
#' If \code{NULL} (the default), use the available evaluation results
#' in \code{primer.data}.
#' @param plot.p.vals Whether p-values from Fisher's exact test
#' should be annotated for every primer set.
#' @return Plot indicating the ratio of primers
#' fulfilling the constraints specified in \code{constraint.settings}
#' for each primer set in \code{primer.data}. 
#' @keywords internal
prepare.constraint.plot <- function(primer.data, 
        constraint.settings, plot.p.vals = FALSE) {
    if (length(primer.data) == 0) {
        return(NULL)
    }
    if (length(constraint.settings) == 0) {
        return(NULL)
    }
    # check type of primers
    primer.classes <- sapply(primer.data, function(x) class(x))
    if (any(primer.classes != "Primers")) {
        stop("Input primers should be lists with objects ",
             "of class 'Primers' ...")
    }
    # determine available columns of the primer sets
    eval.cols <- get.eval.cols(primer.data, constraint.settings)
    if (length(eval.cols) == 0) {
        # nothing to plot
        return(NULL)
    }
    if (length(constraint.settings) != 0) {
        # re-evaluate primer.data with the available input constraints
        eval.names <- gsub("^EVAL_", "", eval.cols)
        primer.data <- eval.comparison.primers(primer.data, constraint.settings[names(constraint.settings) %in% eval.names])
    } 
    # determine frequency of fulfilled constraints
    new.df <- NULL
    p.vals <- rep(NA, length(primer.data))
    for (i in seq_along(primer.data)) {
        primer.df <- asS3(primer.data[[i]])
        if (nrow(primer.df) != 0) {
            if (plot.p.vals) {
                # compute p-vals
                p.data <- primer_significance(primer.df, active.constraints = names(constraint.settings))
                if (length(p.data) != 0) {
                    p.val <- as.numeric(p.data)
                } else {
                    p.val <- NA
                }
                p.vals[i] <- p.val
            } else {
                p.val <- NA
            }
            eval.df <- primer.df[, eval.cols, drop = FALSE]
            if (nrow(eval.df) != 0) {
                df <- data.frame(t(apply(eval.df, 2, function(x) length(which(x))/nrow(eval.df))))
                df <- cbind(Run = unique(primer.df$Run), P_value = p.val, df)
                m.df <- melt(df, c("Run", "P_value"), variable.name = "Constraint", value.name = "Ratio")
            } else {
                m.df <- NULL
            }
            new.df <- rbind(new.df, m.df)
        }
    }
    if (length(new.df) != 0 && length(p.vals) != 0) {
        # adjust/assign p-vals for multiple hypothesis testing
        m <- match(new.df$P_value, p.vals) 
        p.vals <- stats::p.adjust(p.vals, method = "bonf")
        new.df$P_value <- p.vals[m]
    }
    # re-name constraints
    if (length(new.df) != 0) {
        new.df$Constraint <- gsub("EVAL_", "", new.df$Constraint)
        levels <- unique(new.df$Constraint)[order(as.character(unique(new.df$Constraint)))]
        new.df$Constraint <- factor(new.df$Constraint, levels = levels)
        # order Runs by name
        new.df$Run <- factor(new.df$Run, levels = unique(new.df$Run)[order(as.character(unique(new.df$Run)))])
        if (plot.p.vals) {
            p.rep <- format(new.df$P_value, scientific = TRUE, digits = 2)
            new.df$Run <- paste(new.df$Run, " (p-val: ", p.rep, ")", sep = "")
        }
    }
    return(new.df)
}
#' Plot Primer Mismatches
#'
#' Plots primer mismatches for every set.
#'
#' @param primer.data List with primer data frames.
#' @param template.data List with template data frames.
#' @param allowed.mismatches Allowed mismatches.
#' @param highlight.set Primer sets to highlight in the plot.
#' @return Plot of mismatches for comparison.
#' @keywords internal
plot_primer.comparison.mismatches <- function(primer.data, template.data, allowed.mismatches,
        highlight.set = NULL) {
    con.cols <- list("Nbr_of_mismatches_fw", "Nbr_of_mismatches_rev")
    names(con.cols) <- unlist(con.cols)
    con.identifier <- "Number of mismatches"
    new.data <- primer.data
    for (i in seq_along(primer.data)) {
        primer.df <- primer.data[[i]]
        if (any(!"Nbr_of_mismatches_fw" %in% colnames(primer.df))) {
            return(NULL)
        }
        mm.fw <- strsplit(as.character(primer.df$Nbr_of_mismatches_fw), split = ",")
        mm.rev <- strsplit(as.character(primer.df$Nbr_of_mismatches_rev), split = ",")
        r.fw <- do.call(rbind, lapply(seq_along(mm.fw), function(x) primer.df[rep(x, 
            length(mm.fw[[x]])), ]))
        if (length(r.fw) != 0) {
            r.fw$Nbr_of_mismatches_fw <- as.numeric(unlist(mm.fw))
            r.fw$Nbr_of_mismatches_rev <- rep(NA, nrow(r.fw))
        }
        r.rev <- do.call(rbind, lapply(seq_along(mm.rev), function(x) primer.df[rep(x, 
            length(mm.rev[[x]])), ]))
        if (length(r.rev) != 0) {
            r.rev$Nbr_of_mismatches_rev <- as.numeric(unlist(mm.rev))
            r.rev$Nbr_of_mismatches_fw <- rep(NA, nrow(r.rev))
        }
        df <- rbind(r.fw, r.rev)
        if (length(df) != 0 && nrow(df) != 0) {
            new.data[[i]] <- df
        } else {
            # data frame was empty
            new.data[[i]] <- primer.df
        }
    }
    boundaries <- list("Nbr_of_mismatches_fw" = allowed.mismatches,
                        "Nbr_of_mismatches_rev" = allowed.mismatches)
    p <- plot_primer.comparison.box(new.data, con.identifier, con.cols, 
        boundaries, highlight.set = highlight.set)
    return(p)
}

#' Boxplot for Primer Comparison
#'
#' Constructs a box plot showing constraint values for each primer set.
#'
#' @param primer.data List with primer data frames.
#' @param con.identifier Identifier of constraint to be plotted.
#' @param con.cols Column names with the constraint values in the primer data frames.
#' @param boundaries List with constraint settings.
#' @param y.limits Limits for the extent of the y-axis.
#' @param show.points If \code{TRUE} (the default), individual data points
#' are visualized in the boxplot, otherwise they are not shown.
#' @param highlight.set The identifier of a primer set to highlight in the plot.
#' @return A boxplot for primer comparison.
#' @keywords internal
plot_primer.comparison.box <- function(primer.data,
    con.identifier, con.cols, boundaries, y.limits = NULL, 
    show.points = TRUE, highlight.set = NULL) {

    if (length(primer.data) == 0) {
        return(NULL)
    }
    if (length(con.cols) == 0) {
        # no constraints to be plotted according to 'con.cols'
        return(NULL)
    }
    plot.df <- rbind.primer.data(primer.data)
    if (length(plot.df) == 0) {
        return(NULL)
    }
    if (any(unlist(lapply(con.cols, function(x) !x %in% colnames(plot.df))))) {
        return(NULL)
    }
    plot.df <- reshape2::melt(plot.df, id.vars = c("ID", "Run", "Forward", "Reverse"), 
        measure.vars = unlist(con.cols),
        variable.name = "Constraint", 
        value.name = "Value")
    # check whether there's any variable available for plotting ...
    if (length(highlight.set) != 0) {
        # check whether highlight set is specified correctly
        m <- match(highlight.set, plot.df$Run)
        na.idx <- which(is.na(m))
        if (length(na.idx) != 0) {
            msg <- paste("Highlight set not available in data:",
                paste(highlight.set[na.idx], collapse = ","))
            warning(msg)
            highlight.set <- highlight.set[!is.na(m)]
        }
    }
    # Select only fw values for fw primers etc.
    fw.idx <- which(plot.df$Forward != "")
    rev.idx <- which(plot.df$Reverse != "")
    plot.df$Direction <- ifelse(grepl("_fw", plot.df$Constraint), "Forward", ifelse(grepl("_rev", plot.df$Constraint), "Reverse", "Overall"))
    plot.df.fw <- plot.df[plot.df$Direction == "Forward" & plot.df$ID %in% plot.df$ID[fw.idx],]
    plot.df.rev <- plot.df[plot.df$Direction == "Reverse" & plot.df$ID %in% plot.df$ID[rev.idx],]
    plot.df <- rbind(plot.df.fw, plot.df.rev, plot.df[plot.df$Direction == "Overall",])
    ########
    con.names <- unlist(lapply(plot.df$Constraint, function(x) names(con.cols)[sapply(seq_along(con.cols), function(y) any(x %in% con.cols[[y]]))]))
    con.levels <- unique(con.names)[order(unique(con.names))]
    constraint.names <- unlist(constraints_to_unit(con.levels, TRUE))
    plot.df$Constraint <- factor(con.names, levels = con.levels,
                                labels = sapply(constraint.names, deparse))
    #######
    #plot.df$Constraint <- gsub("_fw", "", plot.df$Constraint)
    #plot.df$Constraint <- gsub("_rev", "", plot.df$Constraint)
    title <- "Constraints"
    boundary.data <- NULL # store plot boundaries
    cons <- as.character(levels(plot.df$Constraint))
    #mod.names <- unlist(lapply(con.cols, function(x) gsub("_fw", "", x[1])))
    for (i in seq_along(cons)) {
        con <- cons[i] 
        con.name <- con.levels[i]
        bound <- boundaries[[con.name]]
        if (length(bound) == 0) {
            warning(paste("Constraint boundary for ", con.name, " is not available.", sep = ""))
            next
        }
        cur.bound <- data.frame(Constraint = rep(con, length(bound)),
                                Z = bound)
        boundary.data <- rbind(boundary.data, cur.bound)
    }
    plot.df$Run <- abbreviate(plot.df$Run, getOption("openPrimeR.plot_abbrev"))
    plot.df$Run <- factor(plot.df$Run, levels = unique(plot.df$Run)[order(as.character(unique(plot.df$Run)))])
    p <- ggplot(plot.df, aes_string(x = "Run", y = "Value", colour = "Run")) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle(title) + ylab("Value") +
        geom_boxplot(outlier.shape = NA) +
        guides(colour = FALSE) # remove legend: is redundant here
    if (length(con.identifier) > 1) {
        # multiple constraints are plotted
        p <- p + facet_wrap(stats::as.formula("~Constraint"), 
                scales = "free_y", ncol = 2, labeller = ggplot2::label_parsed)
    } else {
        # a single constraint is plotted
        p <- p + ylab(con.identifier[[1]])
    }
    if (show.points) {
        p <- p + geom_point(alpha = 0.25, position = position_jitter(width = 0.25, 
            height = 0))
    }
    # set colors:
    pal <- getOption("openPrimeR.plot_colors")["Run"] # the RColorBrewer palette to use
    colors <- colorRampPalette(brewer.pal(8, pal))(length(unique(plot.df$Run)))
    p <- p + scale_colour_manual(values = colors)
    if (length(boundary.data) != 0) {
        p <- p + geom_hline(data = boundary.data, aes_string(yintercept = "Z"),
                    size = 0.25, 
                    colour = "red", linetype = 2)
    }
    if (length(y.limits) != 0) {
        p <- p + coord_cartesian(ylim = y.limits)  # don't lose bars outside of limits
    }
    annotation.constraints <- c("cross_dimerization", "self_dimerization", "secondary_structure")
    if (length(con.cols) == 1 && names(con.cols)[1] %in% annotation.constraints) {
        # requires library grid
        if (names(con.cols)[1] %in% c("cross_dimerization", "self_dimerization")) {
            text.high <- "Dimerization\nunlikely"
            text.low <- "Dimerization\nlikely"
        } else if (names(con.cols)[1] == "secondary_structure") {
            text.high <- "Secondary structure\nunlikely"
            text.low <- "Secondary structure\nlikely"
        }
        label.col <- "grey20"
        text_high <- grid::textGrob(text.high, gp=grid::gpar(fontsize=13, fontface="bold", col = label.col ))
        text_low <- grid::textGrob(text.low, gp=grid::gpar(fontsize=13, fontface="bold", col = label.col))
        x.pos <- -1 # position on the x axis of annotation
        high.pos <- max(plot.df$Value, na.rm = TRUE) + 0.1 * min(plot.df$Value, na.rm=TRUE)
        low.pos <- min(plot.df$Value, na.rm = TRUE) - 0.1 * min(plot.df$Value, na.rm= TRUE)
        #######
        p <- p + theme(plot.margin = unit(c(1,1,1,6), "lines")) +
        annotation_custom(text_high, ymin = high.pos, ymax = high.pos, xmin=x.pos,xmax=x.pos) + 
        annotation_custom(text_low,ymin = low.pos, ymax = low.pos,xmin=x.pos,xmax=x.pos) +
        scale_x_discrete(expand = c(0.1,0.0)) +
        geom_segment(aes(x = x.pos, y = high.pos + 0.25 * min(plot.df$Value, na.rm = TRUE), xend = x.pos, 
                    yend = low.pos - 0.25 * min(plot.df$Value, na.rm = TRUE)), color = label.col, size = 1.25,
                    arrow = arrow(length = unit(0.03, "npc")))
    }
    if (length(highlight.set) != 0) {
        # highlight selected sets
        sel <- levels(plot.df$Run) %in% highlight.set
        p <- p + theme(
                axis.text.x = element_text(
                    face=ifelse(sel, "bold","plain"),
                    colour = ifelse(sel, 
                            "grey20", "grey30")))
                    #size = ifelse(sel, 
                            #x.labels.size + 1, x.labels.size)))
    }
    return(p)
}
#' Comparison Coverage Stats.
#'
#' Computes coverage stats for primer comparison.
#'
#' @param primer.data List with primer data frames.
#' @param template.data List with template data frames.
#' @return Coverage statistics for comparing primers.
#' @keywords internal
comparison.cvg <- function(primer.data, template.data) {
    info <- comparison.stats.raw(primer.data, template.data)
    if (length(info) != 0) {
        cvg <- paste(round(info$Coverage * 100, 2), "%", sep = "")
        cvg.string <- paste(info$Nbr_Primers, " primers cover ", 
                    info$Nbr_Covered, 
                    " of ", info$Nbr_Templates, " (", cvg, 
                    ") templates", sep = "")
        result <- data.frame(Run = info$Run, 
                    Coverage = round(info$Coverage, 3), 
                    Info = cvg.string, stringsAsFactors = FALSE)
        return(result)
    } else {
        return(NULL)
    }
}
#' Computation of Raw Stats for Primer Comparison
#' 
#' Computes raw stats for primer comparison.
#'
#' @param primer.data List with primer data frames.
#' @param template.data List with template data frames.
#' @return Raw statistics for primer comparison.
#' @keywords internal
comparison.stats.raw <- function(primer.data, template.data) {
    if (length(primer.data) == 0 || length(template.data) == 0 ||
        is.null(template.data[[1]])) {
        return(NULL)
    }
    run.names <- get.run.names(primer.data)
    N <- unlist(lapply(primer.data, function(x) nrow(x)))
    N.templates <- unlist(lapply(template.data, function(x) nrow(x)))
    nbr.covered <- unlist(lapply(seq_along(primer.data), function(x) length(unique(unlist(covered.seqs.to.idx(primer.data[[x]]$Covered_Seqs, 
        template.data[[x]]))))))
    cvg <- unlist(lapply(seq_along(primer.data), function(x) get_cvg_ratio(primer.data[[x]], 
        template.data[[x]])))
    o <- order(cvg, decreasing = TRUE)
    if (length(cvg) == length(run.names)) {
        # if we have enough data (primers available for each set)
        info <- data.frame(Run = run.names, Coverage = cvg, Nbr_Primers = N, Nbr_Templates = N.templates, 
            Nbr_Covered = nbr.covered, stringsAsFactors = FALSE)
        info <- info[o, ]
        return(info)
    } else {
        return(NULL)
    }
}
build.comparison.table <- function(primers, templates) {
    primers <- set.run.names(primers) # ensure that primer run names are unique
    primer.runs <- get.run.names(primers)
    template.runs <- get.run.names(templates)
    m <- which.max(c(length(primer.runs), length(template.runs)))
    if (m == 1) { # more primer sets than templates
        template.runs <- c(template.runs, rep(NA, length(primer.runs) - length(template.runs)))
    } else { # more templates than primer sets
        primer.runs <- c(primer.runs, rep(NA, length(template.runs) - length(primer.runs)))
    }
    out <- data.frame(PrimerSet = primer.runs, 
                      Templates = template.runs) 
    # numeric index is important for selection!
    rownames(out) <- seq_len(nrow(out))
    if (length(templates) != length(primers)) {
        templates <- replicate(length(primers), templates[1])
    }
    cvg.data <- openPrimeR:::comparison.cvg(primers, templates)
    my.ids <- seq_along(primers)
    if (length(cvg.data) != 0) {
        # match out to cvg.data -> sorted
        m <- match(as.character(out$PrimerSet), as.character(cvg.data$Run)) # TODO: assume that IDs are unique!!!
        add.df <- cvg.data[m,]
        out <- cbind(out, Coverage = add.df$Coverage, Info = add.df$Info)
        # sort by highest coverage
        o <- order(out$Coverage, decreasing = TRUE)
        out <- out[o, ]
        out$Coverage <- paste0(round(out$Coverage * 100, 2), "%")
    }
    return(out)
}
