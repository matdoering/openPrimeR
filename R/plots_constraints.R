####################
# Constraint plots
####################
#' Plot of Template Folding Energies.
#'
#' Plots the DeltaDeltaG of template folding, which is 
#' the difference between the free energy change of
#' the unconstrained folding and the free energy change
#' of the constrained folding.
#'
#' @param fold.df A data frame with free energies for the
#' template regions.
#' @return A plot of DeltaDeltaG.
#' @keywords internal
plot_template_structure <- function(fold.df) {
    if (length(fold.df) == 0) {
        # nothing to plot
        return(NULL)
    }
    # discretize regions:
    interval.df <- fold.df[, c("Constraint_Region_Start", "Constraint_Region_End")]
    intervals <- unlist(lapply(seq_len(nrow(interval.df)), function(x) paste0("[", interval.df[x, 1], ",", interval.df[x, 2], "]")))
    # introduce an interval order:
    o.intervals <- order(interval.df$Constraint_Region_Start)
    intervals <- factor(intervals, levels = unique(intervals[o.intervals]))
    fold.df$Region <- intervals
    ylab <- expression(paste(Delta, " G (no structure) - ", Delta, " G (structure) [kcal/mol]", sep = ""))
    p <- ggplot(fold.df) + 
        geom_boxplot(aes_string(x = "Region", y = "Delta_Delta_G", fill = "Group")) +
        facet_wrap(~Direction, scales = "free_x") +
        ggtitle("Template secondary structures") + 
        ylab(ylab)
    return(p)
}
#' Retrieve data for Constraint Deviations.
#'
#' @param constraint.df An evaluated object of class \code{Primers}.
#' @param constraint.settings A list with settings for the constraints
#' that are to be evalated.
#' @return A data frame providing primer-specific 
#' information on deviations of primer properties from
#' the desired properties.
#' @keywords internal
get_constraint_deviation_data <- function(constraint.df, constraint.settings) {
    if (length(constraint.df) == 0 || nrow(constraint.df) == 0) {
        return(NULL)
    }
    if (!is(constraint.df, "Primers")) {
        stop("Please input a valid primer data frame.")
    }
    constraints <- names(constraint.settings)
    # remove efficiency constraint
    constraints <- constraints[constraints != "primer_efficiency"]
    if (length(constraints) == 0) {
        return(NULL)
    }
    con.idx <- get.constraint.value.idx(constraints, constraint.df)
    con.cols <- lapply(con.idx, function(x) names(constraint.df)[x])
    boundaries <- constraint.settings[constraints]
    constraints <- constraints_to_unit(constraints)
    ######
    plot.df <- reshape2::melt(constraint.df, id.vars = c("ID"), 
        measure.vars = unlist(con.cols),
        variable.name = "Constraint", 
        value.name = "Value")
    # annotate plot.df with 'real constraint names'
    # Select only fw values for fw primers etc.
    fw.idx <- which(constraint.df$Forward != "")
    rev.idx <- which(constraint.df$Reverse != "")
    # determine whether an entry in plot df concerns forward or reverse primers or any type:
    plot.df$Direction <- ifelse(grepl("_fw", plot.df$Constraint), "Forward", 
                            ifelse(grepl("_rev", plot.df$Constraint), "Reverse", "Overall"))
    plot.df.fw <- plot.df[plot.df$Direction == "Forward" & plot.df$ID %in% constraint.df$ID[fw.idx],]
    plot.df.rev <- plot.df[plot.df$Direction == "Reverse" & plot.df$ID %in% constraint.df$ID[rev.idx],]
    plot.df <- rbind(plot.df.fw, plot.df.rev, plot.df[plot.df$Direction == "Overall",])
    con.names <- unlist(lapply(plot.df$Constraint, function(x) names(con.cols)[sapply(seq_along(con.cols), function(y) any(x %in% con.cols[[y]]))]))
    plot.df$Constraint <- factor(con.names, levels = unique(con.names)[order(unique(con.names))])
    # convert to percentage deviations
    myfun <- function(value, con.name, constraint.settings) {
        con.name <- as.character(con.name)
        setting <- constraint.settings[[con.name]]
        if ("min" %in% names(setting) && value < setting["min"]) {
            return((value - setting["min"]) / abs(setting["min"]))
        } else if ("max" %in% names(setting) && value > setting["max"]) {
            return((value - setting["max"]) / abs(setting["max"]))
        } else {
            return(0)
        }
    }
    df <- ddply(plot.df, c("ID", "Constraint", "Value"), plyr::here(summarize),            
                Deviation = myfun(substitute(Value), substitute(Constraint), 
                              constraint.settings))
    df$Run <- constraint.df$Run[1]
    return(df)
}
#' Histogram of Efficiencies
#'
#' Plots a histogram of primer efficiencies.
#'
#' @param primer.df Primer data frame.
#' @param opti.constraints List with constraint settings.
#' @return A plot of primer efficiencies.
#' @keywords internal
plot_constraint.histogram.primer.efficiencies <- function(primer.df, opti.constraints) {
    con.cols <- list("primer_efficiency" = "primer_efficiency")
    con.identifier <- "Primer efficiency"
    if (!"primer_efficiency" %in% colnames(primer.df)) {
        return(NULL)
    }
    eff <- strsplit(primer.df$primer_efficiency, split = ",")
    df <- do.call(rbind, lapply(seq_along(eff), function(x) primer.df[rep(x, length(eff[[x]])), 
        ]))
    df$primer_efficiency <- as.numeric(unlist(eff))
    boundaries <- opti.constraints["primer_efficiency"]
    p <- plot_constraint.histogram(df, con.cols, con.identifier, boundaries)
    return(p)
}
#' Histogram of Number of Mismatches.
#'
#' Plots a histogram of mismatches.
#'
#' @param primer.df Primer data frame.
#' @param allowed.mismatches Number of allowed mismatches.
#' @return A plot of the number of primer mismatches.
#' @keywords internal
plot_constraint.histogram.nbr.mismatches <- function(primer.df, allowed.mismatches) {
    if (length(primer.df) == 0 || nrow(primer.df) == 0) {
        return(NULL)
    }
    con.cols <- c("Nbr_of_mismatches_fw", "Nbr_of_mismatches_rev")
    con.cols <- list("allowed_mismatches" = con.cols)
    con.identifier <- "Number of mismatches"
    if (any(!"Nbr_of_mismatches_fw" %in% colnames(primer.df))) {
        return(NULL)
    }
    mm.fw <- strsplit(as.character(primer.df$Nbr_of_mismatches_fw), split = ",")
    mm.rev <- strsplit(as.character(primer.df$Nbr_of_mismatches_rev), split = ",")
    # could keep identifiers for new primer.df
    r.fw <- do.call(rbind, lapply(seq_along(mm.fw), function(x) primer.df[rep(x, 
        length(mm.fw[[x]])), ]))
    r.fw$Nbr_of_mismatches_fw <- as.numeric(unlist(mm.fw))
    r.fw$Nbr_of_mismatches_rev <- rep(NA, nrow(r.fw))
    r.rev <- do.call(rbind, lapply(seq_along(mm.rev), function(x) primer.df[rep(x, 
        length(mm.rev[[x]])), ]))
    r.rev$Nbr_of_mismatches_rev <- as.numeric(unlist(mm.rev))
    r.rev$Nbr_of_mismatches_fw <- rep(NA, nrow(r.rev))
    df <- rbind(r.fw, r.rev)
    boundaries <- list("allowed_mismatches" = allowed.mismatches)
    p <- plot_constraint.histogram(df, con.cols, con.identifier, boundaries)
    return(p)
}
#' Histogram of Constraints.
#'
#' Plots a histogram of constraint values.
#'
#' @param primer.df Primer data frame, not necessarily a \code{Primers} object.
#' @param con.cols Constraint identifiers in \code{primer.df} to plot.
#' @param con.identifier Name of the constraint to plot.
#' @param boundaries List with constraint settings.
#' @param x.limits Interval limiting the extent of the x-axis.
#' @return A constraint histogram plot.
#' @keywords internal
plot_constraint.histogram <- function(primer.df, con.cols, con.identifier, boundaries = NULL, 
    x.limits = NULL) {

    if (length(primer.df) == 0 || length(con.cols) == 0 || length(con.identifier) == 0) {
        return(NULL)
    }
    if (any(unlist(lapply(con.cols, function(x) !x %in% colnames(primer.df))))) {
        return(NULL)
    }
    plot.df <- reshape2::melt(primer.df, id.vars = c("ID"), 
        measure.vars = unlist(con.cols),
        variable.name = "Constraint", 
        value.name = "Value")
    # Select only fw values for fw primers etc.
    fw.idx <- which(primer.df$Forward != "")
    rev.idx <- which(primer.df$Reverse != "")
    # determine whether an entry in plot df concerns forward or reverse primers or any type:
    plot.df$Direction <- ifelse(grepl("_fw", plot.df$Constraint), "Forward", 
                            ifelse(grepl("_rev", plot.df$Constraint), "Reverse", "Overall"))
    plot.df.fw <- plot.df[plot.df$Direction == "Forward" & plot.df$ID %in% primer.df$ID[fw.idx],]
    plot.df.rev <- plot.df[plot.df$Direction == "Reverse" & plot.df$ID %in% primer.df$ID[rev.idx],]
    plot.df <- rbind(plot.df.fw, plot.df.rev, plot.df[plot.df$Direction == "Overall",])
    con.names <- unlist(lapply(plot.df$Constraint, function(x) names(con.cols)[sapply(seq_along(con.cols), function(y) any(x %in% con.cols[[y]]))]))
    con.levels <- unique(con.names)[order(unique(con.names))]
    constraint.names <- unlist(constraints_to_unit(con.levels, TRUE))
    plot.df$Constraint <- factor(con.names, levels = con.levels,
                                labels = sapply(constraint.names, deparse))
    plot.colors <- c(Forward = brewer.pal(8, "Accent")[5], Reverse = "#0B4D46", Overall = brewer.pal(8, 
        "Accent")[2])
    boundary.data <- NULL # store plot boundaries
    cons <- as.character(levels(plot.df$Constraint))
    for (i in seq_along(cons)) {
        con <- cons[i] 
        con.name <- con.levels[i]
        bound <- boundaries[[con.name]] 
        if (length(bound) == 0) {
            warning(paste("Bound for ", con.name, " not available."))
            next
        }
        cur.bound <- data.frame(Constraint = rep(con, length(bound)),
                                Z = bound)
        boundary.data <- rbind(boundary.data, cur.bound)
    }
    p <- ggplot(plot.df, aes_string(x = "Value", fill = "Direction")) +
        xlab("Value") + 
        ylab("Count") + ggtitle("Constraints") + 
        scale_fill_manual(values = plot.colors) + 
        geom_histogram()
    if (length(con.identifier) > 1) {
        # create multiple facets
        p <- p + facet_wrap(stats::as.formula("~Constraint"), scales = "free_x", ncol = 2,
                            labeller = ggplot2::label_parsed)

    } else {
        # change y axis label to con.identifier
        p <- p + xlab(con.identifier[[1]])
    }
    if (length(boundary.data) != 0) {
        p <- p + geom_vline(data = boundary.data, 
                    aes_string(xintercept = "Z"),
                    size = 0.25, colour = "red", linetype = 2)
    }
    if (length(x.limits) != 0) {
        # don't lose bars outside of limits
        p <- p + coord_cartesian(xlim = x.limits)  
    }
    return(p)
}

#' Retrieval of Evaluation Columns.
#'
#' Retrieves the evaluation columns by intersecting
#' the already evaluated constraints in \code{primer.data}
#' as well as the constraints specified in input constraint settings.
#'
#' @param primer.data A list with \code{Primers} objects.
#' @param constraint.settings A list with constraint settings.
#' @return A character vector with EVAL-columns.
#' @keywords internal
get.eval.cols <- function(primer.data, constraint.settings) {
    available.cols <- unique(unlist(lapply(primer.data, function(x) 
                            colnames(x)[grep("^EVAL_", colnames(x))])))
    eval.names <- gsub("^EVAL_", "", available.cols)
    eval.cols <- intersect(eval.names, names(constraint.settings))
    if (length(eval.cols) != 0) {
        eval.cols <- paste0("EVAL_", eval.cols)
    }
    return(eval.cols)
}

#' Primer Set Statistics
#'
#' Creates an overview of all primer set constraint values.
#'
#' @param primer.df Primer data frame.
#' @param mode.directionality Direction of primers.
#' @param lex.seq Template data frame.
#' @return A data frame with statistics.
#' @keywords internal
primer.set.parameter.stats <- function(primer.df, mode.directionality, lex.seq) {
    # compute summary stats of an evaluated primer set check which stats are
    # available
    if (length(primer.df) == 0 || nrow(primer.df) == 0) {
        return(NULL)
    }
    used.measure <- "IQR"  # IQR (inter-quartile range) or CI (confidence interval)
    eval.col.names <- colnames(primer.df)[grep("^EVAL_", colnames(primer.df))]
    eval.cols <- gsub("EVAL_", "", eval.col.names)
    col.idx <- get.constraint.value.idx(eval.cols, primer.df)
    if (length(col.idx) == 0) {
        warning("No properties available to display stats for: ", primer.df$Run[1])
        return(NULL)
    }
    val.name <- lapply(col.idx, function(x) colnames(primer.df)[x])
    # use frontend identifiers for columns:
    names(val.name) <- constraints_to_unit(names(val.name), use.unit = FALSE)
    single.idx <- which(sapply(val.name, function(x) length(x) == 1))
    double.idx <- which(sapply(val.name, function(x) length(x) == 2))
    # select which direction features to show:
    sel.cols <- val.name[single.idx] # always show features for both directions
    # select by directionality
    if (mode.directionality == "fw") {
        fw.cols <- lapply(seq_along(double.idx), function(x) val.name[[double.idx[x]]][1])
        names(fw.cols) <- names(double.idx)
        sel.cols <- c(sel.cols, fw.cols)
    } else if (mode.directionality == "rev") {
        rev.cols <- lapply(seq_along(double.idx), function(x) val.name[[double.idx[x]]][2])
        names(rev.cols) <- names(double.idx)
        sel.cols <- c(sel.cols, rev.cols) 
    } else { # all columns
        sel.cols <- val.name
    }
    sel.cols <- unlist(sel.cols)  # columns where we can extract the values / create statistics
    names(sel.cols) <- gsub(" ", "_", names(sel.cols)) # gaps are ugly in data frame ..
    sel.cols <- sel.cols[sel.cols %in% colnames(primer.df)]
    # remove "Basic_*" columns
    sel.cols <- sel.cols[!grepl("Basic_", sel.cols)]
    nbr.primers <- nrow(primer.df)
    if (used.measure == "CI") {
        # CI not so appropriate here imo
        means <- plyr::ddply(asS3(primer.df)[, sel.cols], plyr::.(), plyr::numcolwise(mean, na.rm = TRUE))[-1]
        sds <- plyr::ddply(asS3(primer.df)[, sel.cols], plyr::.(), plyr::numcolwise(sd))[-1]  # numcolwise leads to problems with read.csv with NA entries when field is not converted to numeric by read.csv function
        na.idx <- which(is.na(sds))
        if (length(na.idx) != 0) {
            sds[na.idx] <- 0
        }
        error <- qnorm(0.975) * sds/sqrt(nbr.primers)
        conf.left <- round(means - error, 2)
        conf.right <- round(means + error, 2)
        entries <- paste("[", conf.left, ",", conf.right, "]", sep = "")
    } else if (used.measure == "IQR") {
        IQR.data <- plyr::ddply(asS3(primer.df)[, sel.cols], plyr::.(), plyr::numcolwise(quantile, na.rm = TRUE))[-1]
        IQR.l <- IQR.data[2, ]  # 1st quantile
        IQR.r <- IQR.data[4, ]  # 3rd quantile
        entries <- paste("[", round(IQR.l, 2), ",", round(IQR.r, 2), "]", sep = "")
    } else {
        stop("Unknown measure for creating primer stats!")
    }
    #num.cols <- sapply(sel.cols, function(x) is.numeric(primer.df[, x]))
    entry.df <- do.call(cbind, as.list(entries))
    colnames(entry.df) <- names(sel.cols)
    if (length(lex.seq) != 0) {
        cvg <- paste(round(get_cvg_ratio(primer.df, lex.seq), 4) * 
            100, "%", sep = "")
        nbr.total <- nrow(lex.seq)
        nbr.covered <- length(unique(unlist(covered.seqs.to.idx(primer.df$Covered_Seqs, 
            lex.seq))))
        cvg <- paste(nbr.covered, " of ", nbr.total, " (", cvg, ")", sep = "")
    } else {
        cvg <- NA
    }
    # other features from strings: cvg constraints (efficiency, annealing), mismatches, binding
    #######
    # binding range:
    #####
    if (any(primer.df$Relative_Forward_Binding_Position_Start_fw != "")) {
        binding.pos.fw <- quantile(c(as.numeric(unlist(strsplit(primer.df$Relative_Forward_Binding_Position_Start_fw, split = ","))), 
                        as.numeric(unlist(strsplit(primer.df$Relative_Forward_Binding_Position_End_fw, split = ",")))))
        binding.pos.fw <- c(binding.pos.fw[2], binding.pos.fw[4]) # IQR
        binding.pos.fw <- paste0("[", paste0(ifelse(binding.pos.fw <= 0, as.character(binding.pos.fw), paste0("+", binding.pos.fw)), collapse = ","), "]")
    } else {
        binding.pos.fw <- NA
    }
    if (any(primer.df$Relative_Reverse_Binding_Position_Start_rev != "")) {
        binding.pos.rev <- quantile(as.numeric(unlist(strsplit(primer.df$Relative_Reverse_Binding_Position_End_rev, split = ","))), 
                         as.numeric(unlist(strsplit(primer.df$Relative_Reverse_Binding_Position_Start_rev, split = ","))))
        binding.pos.rev <- c(binding.pos.rev[2], binding.pos.rev[4]) # IQR
        binding.pos.rev <- paste0("[", paste0(ifelse(binding.pos.rev <= 0, as.character(binding.pos.rev), paste0("+", binding.pos.rev)), collapse = ","), "]")
    } else {
        binding.pos.rev <- NA
    }
    binding.df <- data.frame(Binding_Range_fw = binding.pos.fw, Binding_Range_rev = binding.pos.rev,
                             stringsAsFactors = FALSE)
                            
    string.constraints <- c("primer_efficiency", "annealing_DeltaG", "coverage_model")
    if (mode.directionality == "both") {
        string.constraints <- c("Nbr_of_mismatches_fw", "Nbr_of_mismatches_rev", string.constraints)
    } else if (mode.directionality == "fw") {
        string.constraints <- c("Nbr_of_mismatches_fw", string.constraints)
    } else {
        string.constraints <- c("Nbr_of_mismatches_rev", string.constraints)
    }
    string.constraints <- string.constraints[string.constraints %in% colnames(primer.df)]
    string.res <- vector("list", length(string.constraints))
    for (i in seq_along(string.constraints)) {
        con <- string.constraints[i]
        entries <- string.to.IQR(paste(primer.df[, con], collapse = ","))
        string.res[[i]] <- entries
    }
    string.df <- do.call(cbind, string.res)
    colnames(string.df) <- string.constraints
    all.entries <- data.frame(cbind(Template_Coverage = cvg, 
                                  Nbr_Primers = nbr.primers, 
                                  binding.df,
                                  string.df, 
                                  entry.df, stringsAsFactors=FALSE), 
                    stringsAsFactors = FALSE)
    # select only non-NA columns
    sel.cols <- apply(all.entries, 2, function(x) any(!is.na(x)))
    all.entries  <- all.entries[, sel.cols]
    #print("set stats:")
    #print(all.entries)
    return(all.entries)
}
#' Conversion of Comma-Separated String to IQR String
#'
#' @param string.values A vector of comma-separated numeric strings.
#' @return The IQR corresponding to the input string.
#' @keywords internal
string.to.IQR <- function(string.values) {
    quant <- lapply(strsplit(string.values, split = ","), function(x) quantile(as.numeric(x), na.rm = TRUE))
    IQR <- lapply(quant, function(x) c(x[2], x[4]))
    entries <- unlist(lapply(IQR, function(x) if (all(is.na(x))) {""} else {paste0("[", paste0(round(x, 2), collapse = ","), "]")}))
    return(entries)
}
#' Plot Dimer DeltaG
#'
#' Plot the distribution of dimerization free energies.
#'
#' @param dimer.data Data frame with dimerization information.
#' @param deltaG.cutoff Free energy cutoff for dimerization.
#' @return A plot of dimerization free energies.
#' @keywords internal
plot.dimer.dist <- function(dimer.data, deltaG.cutoff) {
    if (length(dimer.data) == 0) {
        return(NULL)
    }
    dimer.df <- dimer.data
    colors <- brewer.pal(8, "Pastel1")[c(1, 2)]
    dimer.df$Decision <- "No Dimer"
    dimer.df$Decision[dimer.data$DeltaG < deltaG.cutoff] <- "Dimer"
    dimer.df$Decision <- factor(dimer.df$Decision)
    if (all(dimer.df$Decision == "Dimer")) {
        colors <- colors[1]
    } else if (all(dimer.df$Decision == "No Dimer")) {
        colors <- colors[2]
    }
    p <- ggplot(dimer.df, aes_string(x = "DeltaG")) + geom_histogram(aes_string(fill = "Decision"), 
        colour = "grey") + geom_vline(xintercept = deltaG.cutoff, colour = "red") + 
        ggtitle(expression(paste("Dimerization ", Delta, " G values", 
            sep = ""))) + xlab(expression(paste(Delta, " G", sep = ""))) + scale_fill_manual(values = colors) + 
        theme(axis.text.x = element_text(hjust = 1, vjust = 0.5))
    return(p) 
}
#' Dimerization Table.
#'
#' Summarizes how often individual primers dimerize according to the \code{deltaG.cutoff}.
#'
#' @param dimer.data Data frame with dimerization data.
#' @param deltaG.cutoff Free energy cutoff for dimerization.
#' @param dimer.type String identifying whether \code{dimer.data} refers to cross-dimers or self-dimers?
#' @return Data frame with dimer counts.
#' @keywords internal
dimerization.table <- function(dimer.data, deltaG.cutoff, 
            dimer.type = c("Self-Dimerization", "Cross-Dimerization")) {
    # no need to compute if there's no data and no need to plot for self-dimerization as
    # well (maximum count is 1)
    if (length(dimer.type) == 0) {
        stop("Please provide a dimerization type.")
    }
    dimer.type <- match.arg(dimer.type)
    if (length(dimer.data) == 0 || dimer.type == "Self-Dimerization") {
        return(NULL)
    }
    ID.col <- c("Primer_1", "Primer_2") # only for cross-dimers ..
    dimer.idx <- which(dimer.data$DeltaG < deltaG.cutoff)
    dimer.IDs <- unlist(dimer.data[dimer.idx, ID.col])
    dimer.counts <- data.frame(table(factor(dimer.IDs)))
    if (nrow(dimer.counts) == 0) {
        return(data.frame())
    }
    colnames(dimer.counts)[1] <- "ID"
    dimer.mapping <- lapply(dimer.counts[, 1], function(x) unlist(lapply(ID.col, 
        function(y) which(dimer.data[, y] == x))))
    energy <- sapply(seq_along(dimer.mapping), function(x) mean(dimer.data[dimer.mapping[[x]], "DeltaG"]))
    dimer.counts$Mean_Delta_G <- energy
    o <- order(dimer.counts[, "Freq"], decreasing = TRUE)
    dimer.counts <- dimer.counts[o, ]
    colnames(dimer.counts) <- c("Primer", "Dimer Count", "Mean Free Energy")
    return(dimer.counts)
}

#' Plot of Constraint Values.
#'
#' Shows the distribution of the primer properties.
#' The current constraint settings are indicated with dashed lines in the plot.
#'
#' @param primers Either an evaluated object of class \code{Primers}
#' or a list of \code{Primers} objects.
#' @param settings A \code{DesignSettings} object containing the
#' settings for the constraints to be plotted. 
#' @param active.constraints Identifiers of constraints to be plotted.
#' If \code{active.constraints} is not provided, the plotting method
#' automatically plots all constraints defined in \code{settings}
#' that are annotated in \code{primers}.
#' @param ... \code{highlight.set} (a character vector identifying the set
#' that is to be highlighted when \code{primers} is a list).
#' @return A plot showing the distribution of primer properties.
#' @export
#' @include primers.R settings.R
#' @family constraint visualizations
#' @examples
#' # Plot histogram of constraints for a single primer set
#' data(Ippolito)
#' plot_constraint(primer.df, settings)
#' # Compare constraints across multiple primer sets
#' data(Comparison)
#' plot_constraint(primer.data, settings)
setGeneric("plot_constraint", 
    function(primers, settings, active.constraints = names(constraints(settings)), ...) {
        if (is(settings, "DesignSettings")) {
            # use only the constraints that are defined in the settings
            active.constraints <- active.constraints[active.constraints %in% names(constraints(settings))]
        } else {
            # no settings available -> trust that the input is ok for now
        }
        standardGeneric("plot_constraint")
})
#' Histogram of Constraint Values.
#'
#' Plots a histogram of constraint values.
#'
#' @param primers An evaluated object of class \code{Primers}.
#' @param settings A \code{DesignSettings} object containing the
#' settings for the constraints to be plotted. 
#' @param active.constraints Identifiers of constraints to be plotted.
#' provided settings are used to visualize the desired ranges of constraints.
#' If \code{active.constraints} is not provided, the plotting method
#' will automatically try to plot all constraints defined in \code{settings}.
#' @return A histogram of constraint values for the properties specified by
#' \code{constraints}.
#' @keywords internal
setMethod("plot_constraint", 
    methods::signature(primers = "Primers"), 
    function(primers, settings, active.constraints) {

    if (length(primers) == 0 || nrow(primers) == 0 || length(active.constraints) == 0) {
        return(NULL)
    }
    if (!is(primers, "Primers")) {
        stop("Please input a 'Primers' object for 'primers'.")
    }
    con.idx <- get.constraint.value.idx(active.constraints, primers)
    con.cols <- lapply(con.idx, function(x) names(primers)[x])
    if (length(settings) != 0 && is(settings, "DesignSettings")) {
        boundaries <- constraints(settings)[active.constraints]
    } else {
        boundaries <- NULL
    }
    constraints <- constraints_to_unit(active.constraints)
    p <- plot_constraint.histogram(primers, con.cols, constraints, boundaries)
    return(p)
})

#' Boxplot for Comparing Constraints.
#'
#' Creates a boxplot visualizing the physicochemical properties
#' of multiple primer sets. 
#'
#' @param primers List with evaluated objects of class \code{Primers}. 
#' Each list element corresponds to a single primer set.
#' @param settings A \code{DesignSettings} object containing
#' the constraints to be plotted.
#' @param active.constraints The names of the constraints to be plotted.
#' @param constraint.settings List with settings for each constraint.
#' @param highlight.set Identifiers of primer sets to be highlighted.
#' @return Boxplot comparing the values of the properties specified 
#' by \code{constraints}.
#' @keywords internal
setMethod("plot_constraint", 
    methods::signature(primers = "list"), 
    function(primers, settings, active.constraints,
             highlight.set = NULL) {
    if (length(primers) == 0 || length(active.constraints) == 0) {
        return(NULL)
    }
    # check type of primers
    primer.classes <- sapply(primers, function(x) class(x))
    if (any(primer.classes != "Primers")) {
        stop("Please check whether all supplied primers have the right type.")
    }
    if (length(settings) != 0 && is(settings, "DesignSettings")) {
        boundaries <- constraints(settings)[active.constraints]
    } else {
        boundaries <- NULL
    }
    # consider all constraints that are contained at least once in any of the primer data
    con.cols <- lapply(seq_along(primers), function(x) {
        idx <- get.constraint.value.idx(active.constraints, primers[[x]])
        lapply(idx, function(i) colnames(primers[[x]])[i])
    })
    con.cols <- unlist(con.cols, recursive = FALSE)
    con.cols <- con.cols[!duplicated(names(con.cols))]
    constraints <- constraints_to_unit(active.constraints)
    p <- plot_primer.comparison.box(primers, constraints, 
            con.cols, boundaries, highlight.set = highlight.set)
    return(p)
})
#' Constraint Fulfillment Plot.
#'
#' Visualizes which which primers pass the constraint settings and which
#' primers break the constraints.
#'
#' @param primers Either an object of class \code{Primers} or a list of such objects.
#' @param settings A \code{DesignSettings} object containing
#' the constraints to be evaluated.
#' @param active.constraints The identifiers of constraints to be plotted for fulfillment. By default \code{active.constraints} is set according to 
#' all active constarints defined in \code{settings}.
#' @param plot.p.vals An optional logical argument indicating whether
#' p-values computed via \code{\link{primer_significance}} should be annotated in the plot. The default is \code{FALSE}.
#' @param ... The optional arguments \code{ncol} (a numeric indicating the number of facet columns if \code{primers} is a list),
#' \code{highlight.set} (the identifier of the primer set to be highlighted if \code{primers} is a list)
#' @return A plot indicating the constraints that fulfilled by the input primers.
#' @export
#' @include primers.R settings.R
#' @family constraint visualizations
#' @examples
#' # Plot fulfillment for a single primer set:
#' data(Ippolito)
#' plot_constraint_fulfillment(primer.df, settings)
#' # Plot fulfillment for multiple primer sets:
#' data(Comparison)
#' plot_constraint_fulfillment(primer.data[1:5], settings)
setGeneric("plot_constraint_fulfillment", 
    function(primers, settings, active.constraints = names(constraints(settings)), 
             plot.p.vals = FALSE, ...) {
        standardGeneric("plot_constraint_fulfillment")
})
#' Overview of Constraint Fulfillment.
#'
#' Plots an overview of which primers passed the filtering constraints and which
#' primers did not.
#'
#' @param primers A \code{Primers} object.
#' @param settings A \code{DesignSettings} object.
#' @param active.constraints The identifiers of constarints to be plotted for fulfillment.
#' @param plot.p.vals Show p-value from Fisher's exact test 
#' for the significance
#' of primer constraint fulfillment in comparison to reference
#' primer sets.
#' @return A data frame with statistics on fulfilled constraints.
#' @keywords internal
setMethod("plot_constraint_fulfillment", 
    methods::signature(primers = "Primers"),
    function(primers, settings, active.constraints, plot.p.vals = TRUE) {
    if (length(primers) == 0 || nrow(primers) == 0) {
        return(NULL)
    }
    if (!is(primers, "Primers")) {
        stop("Please input a valid primer data frame.")
    }
    if (!is(settings, "DesignSettings")) {
        stop("Please provide a 'DesignSettings' object for 'settings'.")
    }
    constraint.settings <- constraints(settings)[names(constraints(settings)) %in% active.constraints]
    constraint.df <- asS3(primers)
    # use available EVAL column data
    eval.cols <- get.eval.cols(list(constraint.df), constraint.settings)
    if (length(eval.cols) == 0) {
        # nothing to plot
        return(NULL)
    }
    active.constraints <- names(constraint.settings)
    # re-evaluate primer.data and use constraint.settings constraints
    mode.directionality <- get.analysis.mode(constraint.df)
    eval.df <- eval.constraints(constraint.df, 
                        constraint.settings, active.constraints,
                        mode.directionality, constraint.df)
    constraint.df <- update.constraint.values(constraint.df, eval.df)
    title <- "Constraint fulfillment"
    if (plot.p.vals) {
        # compute p-vals
        p.data <- primer_significance(constraint.df, active.constraints = active.constraints)
        if (length(p.data) != 0) {
            p.val <- as.numeric(p.data)
            title <- paste0(title, " (p-value: ", format(p.val, scientific = TRUE, digits = 2), ")")
        } 
    }
    color.set <- c("#3e6ebc", "#c62f1d")
    names(color.set) <- c("Passed", "Failed")
    eval.df <- data.frame(constraint.df[, eval.cols, drop = FALSE])
    if (all(do.call(c, eval.df), na.rm = TRUE)) {
        # all TRUE?
        colors <- color.set[1]
    } else if (!any(do.call(c, eval.df), na.rm = TRUE)) {
        # all FALSE?
        colors <- color.set[2]
    } else {
        colors <- color.set
    }
    passed.name <- "Passed"
    failed.name <- "Failed"
    eval.df[eval.df == TRUE] <- passed.name  # blue
    eval.df[eval.df == FALSE] <- failed.name  # grey
    eval.df <- cbind(Identifier = constraint.df$Identifier, ID = constraint.df$ID, 
                    eval.df, stringsAsFactors = TRUE)
    colnames(eval.df) <- gsub("EVAL_", "", colnames(eval.df))
    eval.m <- reshape2::melt(eval.df, c("ID", "Identifier"), variable.name = "Constraint", 
        value.name = "Outcome")
    # sort columns by name
    levels <- unique(eval.m$Constraint)[order(as.character(unique(eval.m$Constraint)))]
    eval.m$Constraint <- factor(eval.m$Constraint, 
                         levels = levels)
    constraint.names <- unlist(constraints_to_unit(levels(eval.m$Constraint), FALSE))
    xlab <- "Constraint"
    ylab <- "Primer"
    nbr.constraints <- length(unique(eval.m$Constraint))
    levels(eval.m$ID) <- abbreviate(levels(eval.m$ID), getOption("openPrimeR.plot_abbrev"))
    p <- ggplot(eval.m, aes_string(x = "Constraint", y = "ID")) + geom_tile(aes_string(fill = "Outcome"), 
        colour = "black") + scale_fill_manual(values = colors) + 
        theme(axis.text.x = element_text(angle = 60, 
                hjust = 1)) +
        ggtitle(title) + xlab(xlab) + ylab(ylab) + 
        theme(legend.title = element_text(face = "bold")) + 
        scale_x_discrete(labels = constraint.names) + 
        scale_y_discrete(limits = rev(levels(eval.m$ID))) # top primer should be the first entry in the plot
    return(p)
})
#' Comparison of Evaluation Results.
#'
#' Plots the percentage of primers fulfilling the specified
#' constraints for multiple primer sets.
#'
#' @param primers A list of \code{Primers} objects.
#' @param settings A \code{DesignSettings} object.
#' @param active.constraints The identifiers of constarints to be plotted for fulfillment.
#' @param plot.p.vals Whether p-values from Fisher's exact test
#' should be annotated for every primer set.
#' @param ncol The number of columns for facet wrap.
#' @param highlight.set Identifiers of primer sets to be highlighted.
#' @return Plot indicating the ratio of primers
#' fulfilling the constraints specified in \code{constraint.settings}
#' for each primer set in \code{primers}. 
#' @keywords internal
setMethod("plot_constraint_fulfillment", 
    methods::signature(primers = "list"),
    function(primers, settings, active.constraints, 
             plot.p.vals = FALSE, ncol = 3, highlight.set = NULL) {

    constraint.settings <- constraints(settings)[names(constraints(settings)) %in% active.constraints]
    new.df <- prepare.constraint.plot(primers, constraint.settings, plot.p.vals) 
    if (length(new.df) == 0) {
        # no constraint to plot
        return(NULL)
    }
    if (length(highlight.set) != 0) {
        # check whether highlight set is specified correctly
        m <- match(highlight.set, new.df$Run)
        na.idx <- which(is.na(m))
        if (length(na.idx) != 0) {
            msg <- paste("Highlight set not available in data:",
                paste(highlight.set[na.idx], collapse = ","))
            warning(msg)
            highlight.set <- highlight.set[!is.na(m)]
        }
    }
    constraint.names <- unlist(constraints_to_unit(levels(new.df$Constraint), FALSE))
    pal <- getOption("openPrimeR.plot_colors")["Constraint"] # the RColorBrewer palette to use
    colors <- colorRampPalette(brewer.pal(8, pal))(length(unique(new.df[, "Constraint"])))
    # plot
    p <- ggplot(new.df) + 
            geom_bar(mapping = aes_string(x = "Constraint", y = "Ratio", 
            fill = "Constraint"), position = "dodge", stat = "identity") +
            scale_y_continuous(labels = scales::percent) + ylab("Fulfillment ratio") + 
            xlab("Constraint") + ggtitle("Evaluation of constraints") + 
            theme(
            axis.text.x = element_text(angle = 90, hjust = 1)) +
            facet_wrap(~Run, ncol = ncol) + 
            guides(fill = FALSE) + # remove legend: is redundant here
            scale_x_discrete(labels = constraint.names) + 
           scale_fill_manual(labels = constraint.names, values = colors)

    if (length(highlight.set) != 0) {
        # highlight selected sets (facet part can't be colored!!)
        highlights <- data.frame(Run = highlight.set)
        p <- p + geom_rect(data=highlights,aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf), fill='red', alpha=0.1)
    }
    return(p)
})
#' Plot of Coverage Constraints.
#'
#' Plots the distribution of the coverage constraint values.
#'
#' @param primers A \code{Primers} object or a list with objects of class \code{Primers}. 
#' @param settings A \code{DesignSettings} object.
#' @param active.constraints Names of coverage constraints to be plotted.
#' By default, all active coverage constraints in \code{settings} are plotted.
#' @param ... \code{highlight.set} (a character vector identifying the set
#' that is to be highlighted when \code{primers} is a list).
#' @return A plot showing the distribution of the coverage constraint values.
#' @export
#' @include primers.R settings.R
#' @examples
#' # Plot coverage constraints of a single primer set
#' data(Ippolito)
#' plot_cvg_constraints(primer.df, settings)
#' # Plot coverage constraints for mulitple primer sets
#' data(Comparison)
#' plot_cvg_constraints(primer.data, settings)
setGeneric("plot_cvg_constraints", 
    function(primers, settings, active.constraints = names(cvg_constraints(settings)), ...) {
        if (!is(settings, "DesignSettings")) {
            stop("'settings' should be of class 'DesignSettings'.")
        }
        standardGeneric("plot_cvg_constraints")
})

#' Histogram of Coverage Constraints.
#'
#' Plots a histogram of coverage constraint values.
#'
#' @param primers A \code{Primers} object.
#' @param settings A \code{DesignSettings} object.
#' @param active.constraints Names of coverage constraints to be plotted.
#' @return A plot of coverage constraints.
#' @keywords internal
setMethod("plot_cvg_constraints", 
    methods::signature(primers = "Primers"),
    function(primers, settings, active.constraints) {

    if (!is(primers, "Primers")) {
        stop("'primers' should be of class 'Primers'.")
    }
    if (nrow(primers) == 0) {
        # nothing to plot
        return(NULL)
    }
    cvg.constraints <- cvg_constraints(settings)
    # select only the constraints from the settings that should be plotted
    con.columns <- names(cvg.constraints)[names(cvg.constraints) %in% active.constraints]
    # select only the available constraints
    con.columns <- con.columns[con.columns %in% colnames(primers)]
    if (length(con.columns) == 0) {
        warning("No coverage constraint data available for plotting.")
        return(NULL)
    }
    con.identifier <- constraints_to_unit(con.columns, use.unit = FALSE)
    boundaries <- cvg.constraints[con.columns]
    # prepare every cvg constraint for plotting:
    dfs <- vector("list", length(con.columns))
    for (i in seq_along(con.columns)) {
        con.cols <- con.columns[i]
        res <- strsplit(primers[, con.cols] , split = ",")
        res <- data.frame(as.numeric(unlist(res)))
        colnames(res) <- con.cols
        dfs[[i]] <- res
    }
    # output: data frame with number of rows <=> nbr of total cvg events
    out.df <- asS3(do.call(rbind, lapply(seq_len(nrow(primers)), function(x) primers[rep(x, primers$primer_coverage[x])])))
    con.df <- do.call(cbind, dfs)
    out.df[, con.columns] <- con.df
    out.con.columns <- as.list(con.columns)
    names(out.con.columns) <- con.columns
    p <- plot_constraint.histogram(out.df, out.con.columns, con.identifier, boundaries)
    return(p)
})
#' Plot for Comparing Primer Coverage Constraints.
#'
#' @param primers List with objects of class \code{Primers}. 
#' @param settings A \code{DesignSettings} object.
#' @param active.constraints Names of the coverage constraints to be plotted.
#' @param highlight.set Primer sets to highlight in the plot.
#' @return Plot of primer coverage constraints for multiple sets.
#' @keywords internal
setMethod("plot_cvg_constraints", 
    methods::signature(primers = "list"),
    function(primers, settings, active.constraints, highlight.set = NULL) {

    cvg.constraints <- cvg_constraints(settings)
    # select only the active constraints
    con.columns <- names(cvg.constraints)[names(cvg.constraints) %in% active.constraints]
    if (length(con.columns) == 0) {
        warning("No coverage constraint data available for plotting.")
        return(NULL)
    }
    con.identifier <- constraints_to_unit(con.columns, use.unit = FALSE)
    boundaries <- cvg.constraints[con.columns]
    new.data <- primers
    # transform string columns to numeric by replication:
    for (i in seq_along(primers)) {
        primer.df <- primers[[i]]
        sel.cols <- con.columns[con.columns %in% colnames(primer.df)]
        missing.cols <- con.columns[!con.columns %in% colnames(primer.df)]
        df <- do.call(rbind, lapply(seq_len(nrow(primer.df)), function(x) asS3(primer.df[rep(x, primer.df$primer_coverage[x]), ])))
        if (length(df) != 0 && nrow(df) != 0) {
            # constraint is present in data set
            for (j in (seq_along(sel.cols))) {
                col <- sel.cols[j]
                val <- as.numeric(unlist(strsplit(as.character(primer.df[, col]), split = ",")))
                df[, col] <- val
            }
        } else {
            # constraint is not present
            df <- data.frame()
        }
        df[, missing.cols] <- NA
        new.data[[i]] <- df
    }
    out.con.columns <- as.list(con.columns)
    names(out.con.columns) <- con.columns
    p <- plot_primer.comparison.box(new.data, con.identifier, out.con.columns, 
        boundaries, show.points = FALSE, highlight.set = highlight.set)
    return(p)
})

#' Plot of Constraint Deviations.
#'
#' Shows the deviation of primer properties from
#' the target ranges.
#'
#' Deviations are computed in the following way. Let the
#' minimum and maximum allowed constraint values be given by
#' the interval \eqn{[s, e]} and the observed value be \eqn{p}. Then,
#' if \eqn{p < s}, we output \eqn{-p/|s|}, if \eqn{p > e} we output \eqn{p/|e|},
#' and otherwise, i.e. if \eqn{s <= p <= e}, we output 0.
#'
#' @param primer.data An evaluated object of class \code{Primers} or a list with \code{Primers} objects.
#' @param settings A \code{DesignSettings} object
#' containing the target ranges for the primer properties.
#' @param active.constraints Constraint identifiers to be plotted. By default,
#' all constraints found in \code{settings} are plotted.
#' @param ... \code{deviation.per.primer} (a boolean indicating whether 
#' the deviations should be plotted per primer rather than per constraint
#' if \code{primer.data} is a list)
#' @return A plot showing the deviations of the primer properties from the targets.
#' @export
#' @include primers.R settings.R
#' @family constraint visualizations
#' @examples
#' # Deviations for a single primer set
#' data(Ippolito)
#' plot_constraint_deviation(primer.df, settings)
#' # Deviations for multiple primer sets
#' data(Comparison)
#' plot_constraint_deviation(primer.data, settings)
setGeneric("plot_constraint_deviation", 
    function(primer.data, settings, active.constraints = names(constraints(settings)), ...) {
        if (is(settings, "DesignSettings")) {
            # use only the constraints that are defined in the settings
            active.constraints <- active.constraints[active.constraints %in% names(constraints(settings))]
        } else {
            # no settings available -> trust that the input is ok for now
        }
        standardGeneric("plot_constraint_deviation")
})

#' Plot of Constraint Deviations for a Single Primer Set.
#'
#' Plots a box plot of deviations of 
#' primer properties from the target ranges.
#'
#' @param primer.data An evaluated object of class \code{Primers}.
#' @param settings A \code{DesignSettings} object
#' containing the target ranges for the primer properties.
#' @param active.constraints Constraint identifiers to be plotted.
#' @return A boxplot of deviations
#' @keywords internal
setMethod("plot_constraint_deviation", 
    methods::signature(primer.data = "Primers"),
    function(primer.data, settings, active.constraints) {
        constraint.settings <- constraints(settings)
        constraint.settings <- constraint.settings[names(constraint.settings) %in% active.constraints]
        df <- get_constraint_deviation_data(primer.data, constraint.settings)
        if (length(df) == 0) {
            return(NULL)
        }
        # modify constraint names 
        constraint.names <- unlist(constraints_to_unit(levels(df$Constraint), TRUE))
        title <- "Constraint deviations"
        abs.total.deviation <- mean(abs(df$Deviation)) * 100
        abs.total.deviation <- paste0("mean |deviation| = ", 
                                round(abs.total.deviation, 0), "%")
        title <- paste0(title, " (", abs.total.deviation, ")")
        # set up colors
        pal <- getOption("openPrimeR.plot_colors")["Constraint"] # the RColorBrewer palette to use
        colors <- colorRampPalette(brewer.pal(8, pal))(length(unique(df$Constraint)))
        p <- ggplot(df, aes_string(x = "Constraint", y = "Deviation", colour = "Constraint")) + 
            theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
            ggtitle(title) +
            ylab("Deviation from target") +
            geom_boxplot(outlier.shape = NA) +
            geom_point(alpha = 0.4, position = position_jitter(width = 0.25, 
                height = 0)) + 
            scale_y_continuous(labels = scales::percent) +
            scale_x_discrete(labels = constraint.names) +
            guides(colour = FALSE) + 
            scale_colour_manual(values = colors)
        return(p)
})
#' Plot of Constraint Deviations for Multiple Primer Sets.
#'
#' Plots a box plot of the absolute mean deviation of each primer 
#' for comparing multiple primer sets.
#'
#' @param constraint.df An evaluated object of class \code{Primers}.
#' @param settings A \code{DesignSettings} object
#' containing the target ranges for the primer properties.
#' @param active.constraints Constraint identifiers to be plotted.
#' @param deviation.per.primer Whether to show the deviation per primer or per constraint.
#' @return A boxplot of deviations
#' @keywords internal
setMethod("plot_constraint_deviation", 
    methods::signature(primer.data = "list"),
    function(primer.data, settings, active.constraints, deviation.per.primer = FALSE) {
        # check class of primer.data entries
        c <- sapply(primer.data, class)
        if (!all(c == "Primers")) {
            stop("Please ensure that all 'primer.data' objects are of class 'Primers'.")
        }
        # select active.constraints
        constraint.settings <- constraints(settings)
        constraint.settings <- constraint.settings[names(constraint.settings) %in% active.constraints]
        if (length(constraint.settings) == 0) {
            warning("No active constraint settings available.")
            return(NULL)
        }
        plot.data <- vector("list", length(primer.data))
        for (i in seq_along(primer.data)) {
            primer.df <- primer.data[[i]]
            df <- get_constraint_deviation_data(primer.df, constraint.settings)
            plot.data[[i]] <- df
        }
        plot.data <- do.call(rbind, plot.data)
        if (length(plot.data) == 0) {
            warning("No constraint data available for plotting.") 
            return(NULL)
        }
        # transform to mean absolute deviation 
        ylab <- "Absolute deviation"
        if (deviation.per.primer) {
            # show deviation per primer across all constraints
            plot.data <- plyr::ddply(plot.data, c("ID", "Run"), plyr::summarize, Mean_Deviation = mean(abs(substitute(Deviation))))
            col.id <- "Run"
            title <- "Deviations per primer and set"
        } else {
            # show deviation per considered constraint across all primers
            plot.data <- plyr::ddply(plot.data, c("Constraint", "Run"), plyr::summarize, Mean_Deviation = abs(substitute(Deviation)))
            col.id <- "Constraint"
            title <- "Deviations per constraint and set"
        }
        # set up colors
        pal <- getOption("openPrimeR.plot_colors")[col.id] # the RColorBrewer palette to use
        colors <- colorRampPalette(brewer.pal(8, pal))(length(unique(df[, col.id])))
        p <- ggplot(plot.data, aes_string(x = "Run", y = "Mean_Deviation", colour = col.id)) + 
            theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
            ggtitle(title) +
            ylab(ylab) + 
            # don't show individual points too large...
            geom_boxplot(outlier.size = 0.75) + 
            scale_y_continuous(labels = scales::percent)

        if (deviation.per.primer) {
            # show deviation per primer -> don't show legend for each primer
            p <- p +
                guides(colour = FALSE) 
        } else {
            constraint.names <- constraints_to_unit(levels(plot.data$Constraint), FALSE)
            # change names of constraints
            p <- p + 
                scale_colour_manual(labels = constraint.names, values = colors)
        }
        return(p)
})


