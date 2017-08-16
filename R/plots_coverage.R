#' Plot of Primer Coverage Ratio vs Set Size.
#'
#' Plots the coverage ratios of the input primer sets 
#' against the size of the sets.
#'
#' @param primer.data List with objects of class \code{Primers} containing
#' the primer sets that are to be compared.
#' @param template.data List with objects of class \code{Templates} containing
#' the templates corresponding to \code{primer.data}.
#' @param show.labels Whether the identifiers of the primer sets
#' should be annotated in the plot. The default is \code{TRUE}.
#' @param highlight.set A character vector providing the identifiers
#' of primer sets to highlight. By default, \code{highlight.set} is \code{NULL}
#' such that no highlighting takes place.
#' @return A plot of coverage vs set size.
#' @family comparison visualizations
#' @export
#' @examples
#' data(Comparison)
#' plot_cvg_vs_set_size(primer.data, template.data)
plot_cvg_vs_set_size <- function(primer.data, template.data, show.labels = TRUE, highlight.set = NULL) {
    if (length(primer.data) == 0 || length(template.data) == 0) {
        return(NULL)
    }
    # check type of primer and template data
    template.classes <- sapply(template.data, function(x) class(x))
    primer.classes <- sapply(primer.data, function(x) class(x))
    if (any(template.classes != "Templates") || any(primer.classes != "Primers")) {
        stop("Check types of primers/templates.")
    }
    cvg <- unlist(lapply(seq_along(primer.data), function(x) get_cvg_ratio(primer.data[[x]], 
                                                      template.data[[x]])))
    set.size <- unlist(lapply(primer.data, function(primer.df) nrow(primer.df)))
    set.names <- get.run.names(primer.data)
    plot.df <- data.frame("Run" = set.names, "Coverage" = cvg, "Set_Size" = set.size)
    title <- "Coverage vs set size"
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
    pal <- getOption("openPrimeR.plot_colors")["Run"] # the RColorBrewer palette to use
    colors <- colorRampPalette(brewer.pal(8, pal))(length(unique(plot.df[, "Run"])))
    # determine rate of constraint fulfillment:
    fulfilled.rate <- unlist(lapply(primer.data, function(x) {
        fulfilled.counts <- create_fulfilled_counts(x)
        if (nrow(fulfilled.counts) != 0) {
            fulfilled.counts <- fulfilled.counts[!grepl("Run", colnames(fulfilled.counts))]
            ok.names <- which(!grepl("_failure", colnames(fulfilled.counts)))
            rate <- sum(fulfilled.counts[ok.names]) / (nrow(x) * length(ok.names))
        } else {
            # empty primer set
            rate <- 0
        }
    }))
    #print(fulfilled.rate)
    plot.df$Constraint_Fulfillment_Rate <- fulfilled.rate
    p <- ggplot(plot.df, aes_string(x = "Set_Size", y = "Coverage", size = "Constraint_Fulfillment_Rate"), key = substitute(Run)) + 
    # trans = scales::probability_trans("exp")
        scale_radius(limits = c(0,1), range = c(1,10),
            labels = scales::percent, 
            name = "Constraint Fulfillment") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
        xlab("Number of primers") +
        ggtitle(title) +
        # only show integer values
        scale_x_continuous(breaks = function(x) unique(as.integer(pretty(x)))) +
        scale_colour_manual(values = colors)

    # decide whether to add some jitter
    overplotted <- sapply(seq_len(nrow(plot.df)), function(x) any(plot.df$Set_Size[x] - plot.df$Set_Size[-x] == 0 && abs(plot.df$Coverage[x] - plot.df$Coverage[-x]) < 0.2))
    if (any(overplotted)) {
        # there's overplotting of points -> add some jitter
        p <- p + geom_point(aes_string(fill = "Run", color = NULL), alpha = 0.65, pch=21, position = position_jitter(width = 0.35, height = 0.0))
    } else {
        p <- p + geom_point(aes_string(fill = "Run", color = NULL), alpha = 0.65, pch=21)
    }
    if (show.labels) {
        p <- p + geom_text(aes_string(label = "Run",
                  hjust = 0, vjust = 0), check_overlap = TRUE, size = 3) 
    }
    if (length(highlight.set) != 0) {
        p <- p + geom_point(data = plot.df[plot.df$Run %in% highlight.set,], colour="black", shape=1, size=5, stroke = 3, alpha = 0.6)  
   } 
    return(p)
}
#' Plot of Primer Penalties vs Set Size.
#'
#' Plots the penalties of the input primer sets 
#' against the number of primers contained in each set.
#' The penalties are computed using \code{\link{score_primers}}
#' where more information is provided on how to set \code{alpha}.
#'
#' @param primer.data List with objects of class \code{Primers}.
#' @param settings An object of class \code{DesignSettings}.
#' @param active.constraints A character vector with constraint identifiers
#' to be considered for generating the plot.
#' @param alpha A numeric in the range [0,1] defining the trade-off between
#' the maximal deviation of a constraint (large code{alpha}) and
#' all constraint deviations (large \code{alpha}).
#' By default, \code{alpha} is set to 0 such that the absolute
#' deviation across all constraints is considered.
#' @return A plot showing the association between 
#' primer penalties and the size of the primer sets.
#' @family comparison visualizations
#' @export
#' @examples
#' data(Comparison)
#' plot_deviation_vs_set_size(primer.data, settings)
plot_penalty_vs_set_size <- function(primer.data, settings, 
        active.constraints = names(constraints(settings)), alpha = 0) {
    if (length(primer.data) == 0) {
        return(NULL)
    }
    # check types
    primer.classes <- sapply(primer.data, function(x) class(x))
    if (any(primer.classes != "Primers")) {
        stop("Ensure that 'primer.data' is a list of 'Primers' objects.")
    }
    set.size <- unlist(lapply(primer.data, function(primer.df) rep(nrow(primer.df), nrow(primer.df))))
    penalty.data <- lapply(primer.data, function(primer.df) score_primers(primer.df, settings, active.constraints = active.constraints, alpha = alpha))
    penalties <- unlist(lapply(penalty.data, function(x) x$Penalty))
    run.names <- get.run.names(primer.data)
    run.names <- unlist(lapply(seq_along(run.names), function(x) rep(run.names[x], nrow(primer.data[[x]]))))
    plot.df <- data.frame("Run" = run.names, "Penalty" = penalties, "Set_Size" = set.size)
    title <- "Penalty vs set size"
    pal <- getOption("openPrimeR.plot_colors")["Run"] # the RColorBrewer palette to use
    colors <- colorRampPalette(brewer.pal(8, pal))(length(unique(plot.df[, "Run"])))
    # determine rate of constraint fulfillment:
    p <- ggplot(plot.df, aes_string(x = "Set_Size", y = "Penalty", fill = "Run", group = "Run")) + 
        geom_boxplot() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        xlab("Number of primers") +
        ggtitle(title) +
        # only show integer values
        scale_x_continuous(breaks = function(x) unique(as.integer(pretty(x)))) +
        scale_fill_manual(values = colors)
   return(p)
}

cbind.fill <- function(...){
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

#' Primer Binding Region Data
#'
#' Collects all data concerning primer binding regions.
#'
#' @param primer.df Primer data frame.
#' @param template.df Template data frame.
#' @param direction Primer direction
#' @param group The groups for which binding data shall be retrieved.
#' @param relation Binding region data relative to forward/reverse binding region?
#' @return Data frame with primer binding data.
#' @keywords internal
primer.binding.regions.data <- function(primer.df, template.df, 
    direction = c("both", "fw", "rev"), group = NULL, relation = c("fw", "rev")) {
    if (length(relation) == 0) {
        stop("Please provide the 'relation' arg.")
    }
    relation <- match.arg(relation)
    if (length(direction) == 0) {
        stop("Please provide the 'direction' arg.")
    }
    direction <- match.arg(direction)
    if (length(primer.df) == 0 || nrow(primer.df) == 0) {
        return(NULL)
    }
    if (!"primer_coverage" %in% colnames(primer.df)) {
        stop(primer.df$Run[1], ": Please compute the primer coverage first.")
    }
    if (!is.null(group) && !"all" %in% group) { # select template subset
        idx <- which(template.df$Group %in% group)
        template.df <- template.df[idx,]
    }
    use.location.cols <- NULL
    if (relation == "fw") {
        use.location.cols <- c("Relative_Forward_Binding_Position_Start_fw", "Relative_Forward_Binding_Position_End_fw", 
            "Relative_Forward_Binding_Position_Start_rev", "Relative_Forward_Binding_Position_End_rev")
    } else {
        use.location.cols <- c("Relative_Reverse_Binding_Position_Start_fw", "Relative_Reverse_Binding_Position_End_fw", 
            "Relative_Reverse_Binding_Position_Start_rev", "Relative_Reverse_Binding_Position_End_rev")
    }
    # columns for 'start of binding':
    # 3 prime of primer: start
    # beginning of the target region'
    if (direction == "fw") {
        location.cols <- use.location.cols[c(1,2)]
    } else if (direction == "rev") {
        location.cols <- use.location.cols[c(3,4)]
    } else {
        location.cols <- use.location.cols
    }
    # select primer data belonging to the right group
    idx <- covered.seqs.to.idx(primer.df$Covered_Seqs, template.df)
    dfs <- do.call(rbind, lapply(seq_along(idx), function(x) if (length(idx[[x]]) == 0) 
        NULL else 
            data.frame(ID = primer.df$ID[x], Group = template.df[idx[[x]], "Group"], 
        StringIdx = seq_along(idx[[x]]))))  # id of primer + gene group
    if (any(!location.cols %in% colnames(primer.df))) {
        warning("Primer set doesn't have primers of indicated direction for plotting.")
        return(NULL)
    }
    # remove NA group entries
    if (length(dfs) != 0) {
        dfs <- dfs[!is.na(dfs$Group), ]
    }
    locations <- apply(asS3(primer.df)[, location.cols, drop = FALSE], 2, function(x) strsplit(as.character(x), 
        split = ","))
    all.IDs <- primer.df$ID
    if (length(group) != 0 && !"all" %in% group) {
        all.IDs <- dfs$ID
        m <- match(all.IDs, primer.df$ID)
        ## identify binding positions for the used location cols (fw and rev binding)
        locations <- lapply(seq_along(locations), function(y) lapply(seq_along(m),
                        function(x) {
                            pos <- locations[[y]][[m[x]]]
                            if (length(pos) != 0) {
                                pos[dfs$StringIdx[x]]
                            } else {
                                NULL
                            }
                        }
                    ))
    }
    my.locations <- lapply(locations, function(x) as.numeric(unlist(x)))
    my.IDs <- NULL
    cvg.per.primer <- lapply(locations, function(y) unlist(lapply(y, length)))
    IDs <- unlist(lapply(seq_along(cvg.per.primer), 
        function(x) unlist(lapply(seq_along(cvg.per.primer[[x]]), function(y) rep(all.IDs[y], cvg.per.primer[[x]][y])))))

    if (length(location.cols) == 2) {
        # one direction
        my.IDs <- all.IDs[unlist(lapply(seq_along(locations[[1]]), function(x)
                         rep(x, length(unlist(locations[[1]][x])))))]
        directions <- unlist(lapply(seq_along(cvg.per.primer[[1]]),    
                        function(y) rep(primer.df$Direction[y], cvg.per.primer[[1]][y])))
        location.df <- data.frame(ID = my.IDs, Direction = directions, do.call(cbind, my.locations))
    } else {
        # both directions
        my.IDs.fw <- all.IDs[unlist(lapply(seq_along(locations[[1]]), function(x)
                         rep(x, length(unlist(locations[[1]][x])))))]
         my.IDs.rev <- all.IDs[unlist(lapply(seq_along(locations[[3]]), function(x)
                         rep(x, length(unlist(locations[[3]][x])))))]
        my.IDs <- c(my.IDs.fw, my.IDs.rev)
        directions.fw <- unlist(lapply(seq_along(cvg.per.primer[[1]]),    
                        function(y) rep(primer.df$Direction[y], cvg.per.primer[[1]][y])))
        directions.rev <- unlist(lapply(seq_along(cvg.per.primer[[3]]),    
                        function(y) rep(primer.df$Direction[y], cvg.per.primer[[3]][y])))
        location.df.fw <- data.frame(ID = my.IDs.fw, Direction = directions.fw, cbind(my.locations[[1]], my.locations[[2]]))
        location.df.rev <- data.frame(ID = my.IDs.rev, Direction = directions.rev, cbind(my.locations[[3]], my.locations[[4]]))
        location.df <- rbind(location.df.fw, location.df.rev)
    }
    colnames(location.df) <- c("ID", "Direction", "Start", "End")
    if (relation == "rev") {
        old.cnames <- colnames(location.df)
        colnames(location.df)[old.cnames == "Start"] <- "End"
        colnames(location.df)[old.cnames == "End"] <- "Start"
    }
    if (nrow(location.df) == 0) {
        return(NULL)
    }
    return(location.df)
}

#' Data Preparation for Mismatch Plot.
#' @param primer.df A \code{Primers} object.
#' @param template.df A \code{Templates} object.
#' @param mode Whether to compute for on-target or off-target events.
#' @return A data frame with binding information for every primer.
#' @keywords internal
prepare_mm_plot <- function(primer.df, template.df, 
                    mode = c("on_target", "off_target")) {
    if (length(mode) == 0) {
        mode <- "on_target"
    } else {
        mode <- match.arg(mode)
    }
    all.mm.cols <- c("Mismatch_pos_fw", "Mismatch_pos_rev")
    if (mode == "off_target") {
        all.mm.cols <- paste0("Off_", all.mm.cols)
    }
    all.cvg.cols <- c("Covered_Seqs")
    if (mode == "off_target") {
        all.cvg.cols <- paste0("Off_", all.cvg.cols)
    }
    cvg.criteria <- c("off_annealing_DeltaG", "annealing_DeltaG", "off_primer_efficiency", "primer_efficiency", "coverage_model", "off_coverage_model")
    # get constrained coverage data:
    cvg.defs <- c("constrained", "basic")
    dfs <- vector("list", length(cvg.defs))
    for (i in seq_along(cvg.defs)) {
        # define current columns depending on coverage definition
        cvg.def <- cvg.defs[i]
        if (cvg.def == "basic") { 
            cvg.cols <- paste0("Basic_", all.cvg.cols)
            mm.cols <- paste0("Basic_", all.mm.cols)
        } else {
            cvg.cols <- all.cvg.cols
            mm.cols <- all.mm.cols
        }
        if (any(!mm.cols %in% colnames(primer.df))) {
            # no data for selected coverage definition
            next
        }
        abs.mm.pos.fw <- mismatch.string.to.list(primer.df[, mm.cols[1]])
        abs.mm.pos.rev <- mismatch.string.to.list(primer.df[, mm.cols[2]])
        # determine mismatch positions for forward and reverse primers
        fw.idx <- which(primer.df$Direction == "fw")
        rev.idx <- which(primer.df$Direction == "rev")
        both.idx <- which(primer.df$Direction == "both")
        # for primer pairs: consider the primer with the higher number of mismatches
        both.sel <- sapply(seq_along(both.idx), function(x) if (length(abs.mm.pos.fw[[both.idx[x]]]) > abs.mm.pos.fw[[both.idx[x]]]) {"fw"} else {"rev"})
        # assign 'primer pairs either to fw' or 'rev'
        fw.idx <- c(fw.idx, both.idx[which(both.sel == "fw")])
        rev.idx <- c(rev.idx, both.idx[which(both.sel == "rev")])
        abs.mm.pos <- rep(NA, nrow(primer.df))
        abs.mm.pos[fw.idx] <- abs.mm.pos.fw[fw.idx]
        abs.mm.pos[rev.idx] <- abs.mm.pos.rev[rev.idx]
        pos.3prime <- rep(NA, length(abs.mm.pos))
        pos.3prime[fw.idx] <- lapply(seq_along(fw.idx), function(x) lapply(abs.mm.pos[[fw.idx[x]]], function(y) nchar(primer.df$Forward[x]) - y + 1))
        pos.3prime[rev.idx] <- lapply(seq_along(rev.idx), function(x) lapply(abs.mm.pos[[rev.idx[x]]], function(y) nchar(primer.df$Reverse[x]) - y + 1))
        pos.worst <- lapply(seq_along(abs.mm.pos), function(x) lapply(pos.3prime[[x]], function(y) rep(min(y), length(y))))
        primer.ids <- lapply(seq_along(pos.3prime), function(x) rep(as.character(primer.df$ID[x]), length(unlist(pos.3prime[[x]]))))
        directions <- lapply(seq_along(pos.3prime), function(x) rep(primer.df$Direction[x], length(unlist(pos.3prime[[x]]))))
        cvd.idx <- covered.seqs.to.idx(primer.df[, cvg.cols], template.df)
        if (length(unlist(cvd.idx)) == 0) {
            #warning("Nothing covered")
            next
        }
        template.ids <- lapply(seq_along(pos.3prime), function(x) lapply(seq_along(pos.3prime[[x]]), function(y) rep(template.df$ID[cvd.idx[[x]][y]], length(unlist(pos.3prime[[x]][y])))))
        # count number of mismatches per primer-template pair
        nbr.mismatches <- lapply(seq_along(pos.3prime), function(x) lapply(pos.3prime[[x]], function(y) rep(length(which(!is.na(y))), length(y))))
        # annotate with group of templates
        m <- match(unlist(template.ids), template.df$ID)
        template.groups <- template.df$Group[m]
		#print(primer.ids) # character on windows
        df <- data.frame(Primer = as.character(unlist(primer.ids)), 
                        Direction = unlist(directions),
                        Template = unlist(template.ids),
                        Group = template.groups,
                        Position_3prime = unlist(pos.3prime),
                        Position_3terminus = unlist(pos.worst),
                        Number_of_mismatches = unlist(nbr.mismatches),
                        Coverage_Type = rep(cvg.def, length(template.groups)))
		#print(df$Template)
        # if coverage criteria are available, also include these
        for (z in seq_along(cvg.criteria)) {
            crit <- cvg.criteria[z]
            if (crit %in% colnames(primer.df)) {
                val <- lapply(strsplit(primer.df[, crit], split = ","), as.numeric)
                val <- lapply(seq_along(pos.3prime), function(x) lapply(seq_along(pos.3prime[[x]]), function(y) rep(val[[x]][y], length(pos.3prime[[x]][[y]]))))
                df[, crit] <- unlist(val)
            }
        }
        ####
        # add new hexamer encoding:
        #######
        df$Position_3terminusLocal <- df$Position_3terminus
        df$Position_3terminusLocal[is.na(df$Position_3terminusLocal) | df$Position_3terminusLocal >= 7]  <- 7
        # ensure that position is increasing monotonically: 0 is no mismatch, 6 is mismatch at last 3' hexamer position
        df$Position_3terminusLocal <- abs(df$Position_3terminusLocal - 7)
        # if we have a subsetted template data frame, some primer matches might not be identified -> remove those that shouldn't be considered.
        df <- df[!is.na(df$Template ),]
        dfs[[i]] <- df
    }
    df <- do.call(rbind, dfs)
    # rename off_columns for later use for cvg model
    if (length(df) != 0 && mode == "off_target") {
        # note: this may 'overwrite' the regular columns ..
        colnames(df) <- gsub("^off_", "", colnames(df))
    }
    return(df)
}
#' Retrieval of Template Coverage Data.
#'
#' Determines the coverage of the templates for individual
#' allowed mismatch settings and coverage definitions.
#'
#' @param primer.df A \code{Primers} object.
#' @param template.df A \code{Templates} object.
#' @return Computes a data frame providing the coverage of the 
#' templates for the basic as well as expected (constrained) coverage.
#' @keywords internal
get_template_cvg_data <- function(primer.df, template.df) {
    mode.directionality <- get.analysis.mode(primer.df)
    full.df <- prepare_mm_plot(primer.df, template.df)
    # select only unique primer-template pairs with a certain number of mismatches
    unique.cvg.types <- unique(full.df$Coverage_Type)
    plot.data <- vector("list", length(unique.cvg.types))
    for (z in seq_along(unique.cvg.types)) {
        df <- full.df[full.df$Coverage_Type == unique.cvg.types[z],]    
        # select one event per primer-template pair (worst-case terminal mismatch)
        ddf <- plyr::ddply(df, c("Primer", "Direction", "Template", "Group"), plyr::summarize,
                            Position = unique(substitute(Position_3terminus)), 
                            Number_of_mismatches = unique(substitute(Number_of_mismatches)))
        # select best-case binding mode for each template (lowest number of positions, mismatch furthest from 3 prime)
        dff <- plyr::ddply(ddf, c("Direction", "Template"), function(x) plyr::arrange(x, 
                                            substitute(Number_of_mismatches), -substitute(Position))[1, ])
        # if mode.directionality is "both" -> retain only templates where we have coverage from a fw & rev primer!
        if (mode.directionality == "both") {
            dir.count <- plyr::ddply(ddf, "Template", plyr::summarize, DirectionCount = length(unique(substitute(Direction))))
            # define template coverage events to remove:
            rm.template.id <- dir.count$Template[dir.count$DirectionCount <= 1]
            m <- match(rm.template.id, dff$Template)
            if (length(m) != 0) {
                dff <- dff[-m,]
            }
        }
        # if there's no mismatches, plot these primers with mismatch position @ position 0
        dff[is.na(dff$Position), "Position"] <- 0
        dff$Status <- unique.cvg.types[z]
        plot.data[[z]] <- dff
    }
    plot.df <- do.call(rbind, plot.data)
    return(plot.df)
}
#' Statistics on the Number of Coverage Events per Primer.
#'
#' Creates a table summarizing the coverage events
#' of each primer according to the number of mismatches
#' between primers and templates.
#'
#' Entries in numeric table columns
#' indicate the percentage of coverage events occurring
#' with a certain number of mismatches. For example
#' column \emph{3} provides all coverage events
#' with exactly three mismatches between primers and templates.
#' The column \emph{Group_Coverage} provides a listing
#' of the percentage of covered templates per group.
#'
#' @param primer.df A \code{Primers} object with evaluated coverage providing
#' the set of primers for which the coverage statistics shall be computed.
#' @param template.df A \code{Templates} object providing
#' the template sequences for which the primer coverage
#' has been computed.
#' @param cvg.definition If \code{cvg.definition} is set to
#' "constrained", the statistics for the expected
#' coverage (after applying the coverage constraints) are retrieved.
#' If \code{cvg.definition} is set to "basic", the coverage is determined 
#' solely by string matching (i.e. without applying the coverage constraints).
#' By default, \code{cvg.definition} is set to "constrained".
#' @return A data frame listing the number of binding events
#' broken down according to the number of expected mismatches between
#' primers and templates.
#' @export
#' @examples
#' data(Ippolito)
#' primer.cvg.stats <- get_cvg_stats_primer(primer.df, template.df)
get_cvg_stats_primer <- function(primer.df, template.df,
                                cvg.definition = c("constrained", "basic")) {

    cvg.definition <- match.arg(cvg.definition)
    full.df <- prepare_mm_plot(primer.df, template.df)
    # select only constrained cvg events
    full.df <- full.df[full.df$Coverage_Type == cvg.definition, ]
    df <- plyr::ddply(full.df, c("Primer", "Template", "Group"), plyr::summarize,
                            Position = unique(substitute(Position_3terminus)), 
                            Number_of_mismatches = unique(substitute(Number_of_mismatches)))
    if (length(df) == 0 || nrow(df) == 0) {
        # no stats available
        return(NULL)
    }
    mm.range <- seq(0, max(df$Number_of_mismatches))
    # create a matrix containing the coverage of every primer for every mismatch setting
    mm.stats <- matrix(rep(0, length(mm.range) * nrow(primer.df)), nrow = nrow(primer.df), ncol = length(mm.range))
    for (i in seq_along(mm.range)) {
        mm <- mm.range[i]
        sub.df <- df[df$Number_of_mismatches == mm,]
        count.df <- plyr::ddply(sub.df, c("Primer"), plyr::summarize,
                            Coverage = length(unique(substitute(Template))))
        m <- match(count.df$Primer, primer.df$ID)
        mm.stats[m,i] <- count.df$Coverage
    }
    ####
    # annotate coverage counts per mismatch number with percentages
    ####
    mm.stats.x <- mm.stats / rowSums(mm.stats)
    mm.stats.x[is.na(mm.stats.x)] <- 0
    for (i in seq_len(nrow(mm.stats))) {
        mm.stats[i, ] <- paste0(mm.stats[i,], " (", round(mm.stats.x[i,], 3) * 100, "%)")
    }
    mm.stats <- data.frame(mm.stats)
    colnames(mm.stats) <- mm.range
    ########
    # enrich table with other information
    ########
    # what is the main group of covered templates per primer?
    count.df <- plyr::ddply(df, c("Primer", "Group"), plyr::summarize,
                            Coverage = length(unique(substitute(Template))))
    # order by largest number of coverage events per group
    count.df <- count.df[order(count.df$Coverage, decreasing = TRUE),]
    count.df$Group_Coverage <- sapply(seq_len(nrow(count.df)), function(x) count.df$Coverage[x] / length(which(template.df$Group == count.df$Group[x])))
    group.cvg <- unlist(lapply(primer.df$ID, function(x) {
                        idx <- which(as.character(count.df$Primer) == as.character(x))
                        strings <- sapply(idx, function(y) paste0(count.df[y, "Group"], " (",
                                                       paste0(round(count.df[y, "Group_Coverage"], 3) * 100, "%", 
                                                       ")", sep = "")))
                        res <- paste0(strings, collapse = ",")
                        }))
    #mm.stats <- cbind("Primer" = rownames(mm.stats), "Forward" = primer.df$Forward, "Reverse" = primer.df$Reverse, "Group_Coverage" = group.cvg, mm.stats)
    mm.stats <- cbind("Primer" = primer.df$ID, "Group_Coverage" = group.cvg, mm.stats)
    return(mm.stats)
}
#' Plot of Primer Subset Coverage.
#'
#' Visualizes the coverage of optimized primer subsets.
#'
#' The input for the \code{primer.subsets} argument can be computed using
#' \code{\link{subset_primer_set}}. 
#' The line plot indicates the ratio of covered templates when considering
#' all primers in a primer set of a given size. The bar plots indicate
#' the coverage ratios of individual primers in a set. The target coverage
#' ratio is indicated by a horizontal line. Bars exceeding the target ratio
#' possibly indicate the existence of redundant coverage events.
#'
#' @param primer.subsets A list with optimal primer subsets, each of class \code{Primers}. The \emph{k}-th list entry should correspond to an object of class \code{Primers}
#' representing the primer subset of size \emph{k} whose coverage ratio
#' is the largest among all possible subsets of size \emph{k}.
#' @param template.df An object of class \code{Templates} containing the
#' template sequences corresponding to the primers specified in \code{primer.subsets}.
#' @param required.cvg The required coverage ratio.
#' The default is 100\%; this value is plotted as a horizontal line.
#' @return Plot of the coverages of the primer subsets in \code{primer.subsets}.
#' @export
#' @family coverage visualizations
#' @examples
#' data(Ippolito)
#' primer.subsets <- subset_primer_set(primer.df, template.df)
#' plot_primer_subsets(primer.subsets, template.df)
plot_primer_subsets <- function(primer.subsets, template.df, required.cvg = 1) {
                            
    if (length(primer.subsets) == 0) {
        return(NULL)
    }
    if (!all(sapply(primer.subsets, function(x) is(x, "Primers")))) {
        stop("All subsets should represent primer sets.")
    }
    if (!is(template.df, "Templates")) {
        stop("Please provide a valid template data frame.")
    }
    # create plot.df: need cvg and size of set
    N <- sapply(primer.subsets, nrow)
    cvg <- sapply(seq_along(primer.subsets), function(x) get_cvg_ratio(primer.subsets[[x]], 
        template.df))
    # retrieve the IDs of the primer(s) that are different in the next subset size
    IDs <- lapply(primer.subsets, function(x) x$ID)
    # collect all set members?
    collect.all <- TRUE
    if (!collect.all) {
        # only consider primers that are different from previous ones
    IDs <- lapply(seq_along(IDs), function(x) {
        if (x == 1) {
            out <- IDs[[x]]
        } else {
            out <- setdiff(IDs[[x]], IDs[[x-1]])
        }
        })
    }
    primer.df <- tail(primer.subsets, n = 1)[[1]] # the complete primer set
    individual.cvg <- lapply(IDs, function(x) sapply(x, function(y) get_cvg_ratio(primer.df[match(y, primer.df$ID),], template.df))) # cvg of added primers for this subset size
    #IDs <- unlist(lapply(IDs, function(x) paste(x, collapse = "\n")))
    subset.df <- data.frame(Set_Size = N, Coverage = cvg, Type = "Overall")
    single.df <- data.frame(Set_Size = unlist(lapply(seq_along(N), function(x) rep(N[x], length(individual.cvg[[x]])))),
                            Single_Coverage = unlist(individual.cvg),
                            ID = unlist(IDs), Type = "Individual")
    set.labels <- seq_along(primer.subsets)
    max.cvg <- sum(single.df[single.df$Set_Size == tail(N, n = 1), "Single_Coverage"])
    # order primers by time of addition to set
    levels <- rev(unique(unlist(lapply(primer.subsets, function(x) x$ID)))) # first primers at the bottom
    single.df$ID <- factor(single.df$ID, levels = levels, labels = abbreviate(levels, getOption("openPrimeR.plot_abbrev")))
    pal <- getOption("openPrimeR.plot_colors")["Primer"] # the RColorBrewer palette to use
    colors <- colorRampPalette(brewer.pal(8, pal))(length(unique(single.df[, "ID"])))
    indi.color <- ifelse(length(colors) >= 1, colors[1], NA)
    p <- ggplot() + 
    # add bars for individual cvg
        geom_bar(data = single.df, stat = 'identity', 
                aes_string(x = "Set_Size", y = "Single_Coverage", 
                           fill = "ID", colour = "Type"), 
                alpha = 0.80) +
        scale_color_manual(values=c("Individual"= NA, "Overall" = "grey30"),
            guide = ggplot2::guide_legend(override.aes = list(linetype = c(0, 1), 
                                fill = c(indi.color, NA), 
                                colour = c(NA, "grey20")),
                                title = "Coverage type")) +
        scale_fill_manual(values = colors) +
        geom_point(data = subset.df, size = 1.5, aes_string(x = "Set_Size", y = "Coverage", colour = "Type"))  +
        geom_line(data = subset.df, size = 1, aes_string(x = "Set_Size", 
                  y = "Coverage", colour = "Type")) +
        xlab("Subset size") + ylab("Coverage of templates") + 
        ggtitle("Primer subset coverage") + 
        scale_y_continuous(labels = scales::percent, 
                            limits = c(0, max(max.cvg, 1))) + 
        scale_x_continuous(breaks = seq_along(primer.subsets),
                           labels = set.labels) +
        geom_hline(yintercept = required.cvg, linetype = 2, colour = "red") + # indicate the desired cvg
        theme(axis.text.x = element_text(angle = 60, hjust = 1))
    if (length(unique(single.df$ID)) > 30) {
        # don't show individual primer legend if nbr of primers exceeds 30
        p <- p + guides(fill = FALSE)
    }
    return(p)
}

#' Data for Primer Plot.
#'
#' Constructs a data frame containing information
#' about primer binding events.
#'
#' @param primer.df An object of class \code{Primers} containing
#' primers with evaluated primer coverage.
#' @param template.df An object of class \code{Templates} with template sequences
#' corresponding to \code{primer.df}.
#' @param identifier Identifiers of primers that are to be considered.
#' If \code{identifier} is set to \code{NULL} (the default), all primers are considered.
#' @param relation Compute binding positions relative to forward (\code{fw}) or reverse (\code{rev}) binding regions.
#' The default is "fw".
#' @return Data frame with primer binding data.
#' @keywords internal
get_plot_primer_data <- function(primer.df, template.df, identifier = NULL, relation = c("fw", "rev")) {

    if (!is(primer.df, "Primers")) {
        stop("Please input a valid primer data frame.")
    }
    if (!is(template.df, "Templates")) {
        stop("Please provide a valid template data frame.")
    }
    if (!"primer_coverage" %in% colnames(primer.df)) {
        stop(primer.df$Run[1], ": Please compute the primer coverage first.")
    }
    if (is.null(identifier)) { # plot all primers
        identifier <- primer.df$Identifier
    }
    if (length(relation) == 0) {
        stop("Please provide the 'relation' arg.")
    }
    relation <- match.arg(relation)
    m <- match(identifier, primer.df$Identifier)  # determine primers to be plotted
    if (any(is.na(m))) {
        stop("Could not find specified identifiers.")
    }
    if (length(m) == 0 || is.na(identifier) || identifier == "") {
        # nothing to plot
        return(NULL)
    }
    primer.info <- primer.df[m, ]  # info about selected primers
    s.idx <- lapply(covered.seqs.to.idx(primer.info$Covered_Seqs, template.df), function(x) if (length(x) != 0) {which(!is.na(x))} else {NULL})  # idx of of covered seqs per primer
    sel.idx <- which(!is.na(unlist(covered.seqs.to.idx(primer.info$Covered_Seqs, 
        template.df))))  # unlist idx of all primers covering something
    # need to find the idx of fw covering and rev covering seqs
    fw.idx <- which(primer.info$Forward != "")
    rev.idx <- which(primer.info$Reverse != "")
    if (any(primer.info$Forward != "")) {
        sel.idx.fw <- which(!is.na(unlist(covered.seqs.to.idx(primer.info$Covered_Seqs[fw.idx], 
            template.df))))
        s.idx.fw <- sapply(covered.seqs.to.idx(primer.info$Covered_Seqs[fw.idx], 
            template.df), function(x) if (length(x) != 0) {which(!is.na(x))} else {NULL})
    } else {
        sel.idx.fw <- NULL
        s.idx.fw <- NULL
    }
    if (any(primer.info$Reverse != "")) {
        sel.idx.rev <- which(!is.na(unlist(covered.seqs.to.idx(primer.info$Covered_Seqs[rev.idx], 
            template.df))))
        s.idx.rev <- sapply(covered.seqs.to.idx(primer.info$Covered_Seqs[rev.idx], 
            template.df), function(x) if (length(x) != 0) {which(!is.na(x))} else {NULL})
    } else {
        sel.idx.rev <- NULL
        s.idx.rev <- NULL
    }
    # retrieve coverage info:
    covered.seqs <- lapply(seq_along(primer.info$Covered_Seqs), function(x) as.numeric(unlist(strsplit(primer.info$Covered_Seqs[x], 
        split = ",")[s.idx[[x]]])))
    # determine relative primer starts and ends depending on whether we want to
    # analyze regarding fw allowed region or rev allowed region
    if (relation == "fw") {
        primer.start <- lapply(1:nrow(primer.info), function(x) as.numeric(unlist(strsplit(primer.info$Relative_Forward_Binding_Position_Start_fw[x], 
            split = ","))[s.idx[[x]]]))
        primer.end <- lapply(1:nrow(primer.info), function(x) as.numeric(unlist(strsplit(primer.info$Relative_Forward_Binding_Position_End_fw[x], 
            split = ","))[s.idx[[x]]]))
        primer.start.rev <- lapply(1:nrow(primer.info), function(x) as.numeric(unlist(strsplit(primer.info$Relative_Forward_Binding_Position_End_rev[x], 
            split = ","))[s.idx[[x]]]))
        primer.end.rev <- lapply(1:nrow(primer.info), function(x) as.numeric(unlist(strsplit(primer.info$Relative_Forward_Binding_Position_Start_rev[x], 
            split = ","))[s.idx[[x]]]))
    } else {
        primer.start <- lapply(1:nrow(primer.info), function(x) as.numeric(unlist(strsplit(primer.info$Relative_Reverse_Binding_Position_Start_fw[x], 
            split = ","))[s.idx[[x]]]))
        primer.end <- lapply(1:nrow(primer.info), function(x) as.numeric(unlist(strsplit(primer.info$Relative_Reverse_Binding_Position_End_fw[x], 
            split = ","))[s.idx[[x]]]))
        #print("primer posis:")
        #print(primer.start)
        #print(primer.end)
        primer.start.rev <- lapply(1:nrow(primer.info), function(x) as.numeric(unlist(strsplit(primer.info$Relative_Reverse_Binding_Position_End_rev[x], 
            split = ","))[s.idx[[x]]]))
        primer.end.rev <- lapply(1:nrow(primer.info), function(x) as.numeric(unlist(strsplit(primer.info$Relative_Reverse_Binding_Position_Start_rev[x], 
            split = ","))[s.idx[[x]]]))
        
    }
    covered.seq.identifiers <- as.numeric(unlist(strsplit(primer.info$Covered_Seqs, 
        split = ",")))[sel.idx]
    # determine nbr covered for given template.df
    nbr.covered <- sapply(seq_along(1:nrow(primer.info)), function(x) length(strsplit(primer.info$Covered_Seqs[x], 
        split = ",")[[1]][s.idx[[x]]]))  # nbr of covered seqs per primer
    m <- match(covered.seq.identifiers, template.df$Identifier)
    not.na.mapping <- which(!is.na(m))
    m <- m[not.na.mapping]
    covered.seqs <- template.df[m, ]
    not.m <- setdiff(1:nrow(template.df), m)
    uncovered.seqs <- template.df[not.m, ]
    repeat.idx <- unlist(sapply(seq_along(nbr.covered), function(x) rep(x, each = nbr.covered[x])))  # assignment of index in primer.info to all covered seqs from each primers
    primer.info.m <- primer.info[repeat.idx, ]
    total.covered <- sum(nbr.covered)
    total.y <- total.covered + nrow(template.df)  # dimension of plot: nbr of primers + nbr of seqs
    plot.df <- data.frame(Type = rep(NA, total.y), y = rep(NA, total.y), 
                          Map = rep(NA, total.y))
    y.idx <- 1
    mode.directionality <- get.analysis.mode(primer.info)
    for (i in 1:nrow(template.df)) {
        id <- template.df$Identifier[i]
        # determine nbr of occurrences of this seq in covered.seq.identifiers vector
        count.idx <- which(covered.seq.identifiers == id)  # how often was this sequence covered?
        p.idx <- repeat.idx[count.idx] # index in primer info
        covered <- FALSE
        if (mode.directionality == "both") {
            # require at least 1 fw and 1 rev primer
            fw.count <- length(which(primer.info[p.idx,]$Forward != ""))
            rev.count <- length(which(primer.info[p.idx,]$Reverse != ""))
            count <- min(c(fw.count, rev.count))
        } else {
            # require any primer to cover a sequence
            count <- length(count.idx)
        }
        if (count != 0) { # 
            covered <- TRUE
        }
        seq.idx <- y.idx
        primer.idx <- (seq.idx + 1):(seq.idx + 1 + length(count.idx) - 1)
        plot.df$Type[seq.idx] <- ifelse(covered, "Covered Sequence", "Not-covered Sequence")
        plot.df$y[seq.idx] <- seq.idx
        plot.df$Map[seq.idx] <- i
        if (length(count.idx) != 0) {
            # determine fw and rev primers
            plot.df$Type[primer.idx] <- "Primer"
            plot.df$y[primer.idx] <- primer.idx
            plot.df$Map[primer.idx] <- p.idx  # gives the index of primer in primer.info
        } 
        y.idx <- y.idx + 1 + length(count.idx)
    }
    # enrich with info
    primer.data.idx <- plot.df$Map[which(plot.df$Type == "Primer")]
    pos.idx <- get.unlist.idx(primer.start, primer.data.idx)
    covered.seq.identifiers <- sapply(seq_along(primer.info$Covered_Seqs), function(x) as.numeric(strsplit(primer.info$Covered_Seqs[x], 
        split = ",")[[1]][s.idx[[x]]]))
    primer.cov.idx <- sapply(covered.seq.identifiers, function(x) match(x, template.df$Identifier))
    primer.gene.groups <- sapply(primer.cov.idx, function(x) template.df[x, "Group"])
    fw.primer.name <- "fw primer"
    rev.primer.name <- "rev primer"
    primer.names <- c(fw.primer.name, rev.primer.name)
    primer.data <- data.frame(Identifier = as.character(primer.info$Identifier[primer.data.idx]), 
                    ID = as.character(primer.info$ID[primer.data.idx]),
                    stringsAsFactors = FALSE) 
    primer.data <- cbind(primer.data, x_start = unlist(primer.start)[pos.idx][not.na.mapping], 
        x_end = unlist(primer.end)[pos.idx][not.na.mapping], x_start_rev = unlist(primer.start.rev)[pos.idx][not.na.mapping], 
        x_end_rev = unlist(primer.end.rev)[pos.idx][not.na.mapping], Group = unlist(primer.gene.groups)[pos.idx][not.na.mapping])
    seq.data.idx <- plot.df$Map[which(plot.df$Type != "Primer")]
    #print("primer data:")
    #print(primer.data)
    # n.b.: not using initial here anymore -> plot relates to the current binding region
    seq.data <- data.frame(Identifier = template.df$Identifier[seq.data.idx], 
        ID = template.df$ID[seq.data.idx], 
        x_start = -template.df$Allowed_End_fw[seq.data.idx], 
        x_end = nchar(template.df$Sequence)[seq.data.idx] - 
                template.df$Allowed_End_fw[seq.data.idx] - 1, 
        x_start_rev = NA,  # only relevant for primers
        x_end_rev = NA, 
        Group = template.df$Group[seq.data.idx],
        stringsAsFactors = FALSE)
    if (relation == "fw") {
        seq.data$x_start <- -template.df$Allowed_End_fw[seq.data.idx]
        seq.data$x_end <- nchar(template.df$Sequence)[seq.data.idx] - 
                            template.df$Allowed_End_fw[seq.data.idx] - 1 
    } else {
        seq.data$x_start <- template.df$Allowed_Start_rev[seq.data.idx] - nchar(template.df$Sequence)[seq.data.idx]
        seq.data$x_end <- nchar(template.df$Sequence[seq.data.idx]) - seq.data$x_start - 1
    }
    add.df <- rbind(primer.data, seq.data)[FALSE, ]  # primer data + seq data should make up the plot.df
    add.df[which(plot.df$Type == "Primer"), ] <- primer.data
    add.df[which(plot.df$Type != "Primer"), ] <- seq.data
    rownames(add.df) <- NULL
    rownames(plot.df) <- NULL
    plot.data <- cbind(plot.df, add.df)
    plot.data$Type <- factor(plot.data$Type, levels = c("Primer", "Covered Sequence", 
        "Not-covered Sequence"))
    ########## determine limits of plot
    #seq.idx <- which(plot.data$Type != "Primer")
    #s.min <- min(plot.data[seq.idx, "x_start"], na.rm = TRUE)
    #seq.max <- max(plot.data[seq.idx, "x_end"], na.rm = TRUE)
    #primer.idx <- which(plot.data$Type == "Primer")
    #x.lim <- c(s.min, seq.max)  # max was p.max before
    #print("before plot lim adjustment")
    #print(plot.data)
    #plot.data$x_start <- sapply(plot.data$x_start, function(x) max(x, x.lim[1]))
    #plot.data$x_end <- sapply(plot.data$x_end, function(x) min(x, x.lim[2]))
    #plot.data$x_start_rev <- sapply(plot.data$x_start_rev, function(x) max(x, x.lim[1]))
    #plot.data$x_end_rev <- sapply(plot.data$x_end_rev, function(x) min(x, x.lim[2]))
    ########## 
    cur.plot.df <- plot.data
    cur.plot.df$ID <- as.character(cur.plot.df$ID)
    # plot only the templates
    rel.cols <- c("ID", "y", "Type", "x_start", "x_end", "x_start_rev", "x_end_rev")
    m <- melt(cur.plot.df[, rel.cols], id = c("y", "Type", "ID"))
    # remove NA entries
    na.idx <- which(is.na(m[, "value"]))
    if (length(na.idx) != 0) {
        m <- m[-na.idx, ]
    }
    # update types
    types <- ifelse(m$variable == "x_start_rev" | m$variable == "x_end_rev", "Reverse Primer", 
        ifelse(m$Type == "Primer", "Forward Primer", as.character(m$Type)))
    m$Type <- types
    m$variable <- ifelse(grepl("start", m$variable), "x_start", "x_end")
    d <- dcast(m, ID + y + Type ~ variable)
    return(d)
}
#' Primer View Plot.
#' 
#' Visualizes the binding positions of every primer relative to
#' the target binding region in the corresponding template sequences.
#'
#' @param primer.df An object of class \code{Primers} containing
#' primers with evaluated primer coverage.
#' @param template.df An object of class \code{Templates} with template sequences
#' corresponding to \code{primer.df}.
#' @param identifier Identifiers of primers that are to be considered.
#' If \code{identifier} is set to \code{NULL} (the default), all primers are considered.
#' @param relation Compute binding positions relative to forward ("fw") or reverse ("rev") binding regions.
#' The default is "fw".
#' @param region.names Character vector of length 2 providing the names
#' of the binding and amplification region.
#'
#' @return A plot of primer binding sites in the templates.
#' @export
#' @family coverage visualizations
#' @examples
#' data(Ippolito)
#' plot_primer(primer.df[1,], template.df[1:30,])
plot_primer <- function(primer.df, template.df, identifier = NULL, 
                        relation = c("fw", "rev"), 
                        region.names = c("Binding region", "Amplification region")) {

    if (length(region.names) != 2) {
        stop("Need 2 region names.")
    }
    relation <- match.arg(relation)
    d <- get_plot_primer_data(primer.df, template.df, identifier = identifier, relation = relation)
    # region annotation:
    y.val <- 0.05 * nrow(d)
    y.ext <- min(10, y.val)
    ymin <- min(-y.ext, -3)
    ymax <- min(-1, -0.01 * nrow(d))
    region.df <- create_region_boxes(list(primer.df), list(template.df), relation, region.names, ymin, ymax, max(d$x_end, na.rm = TRUE))

    if (length(d) == 0 || nrow(d) == 0) {
        # nothing to plot
        return(NULL)
    }
    col <- c("Forward Primer" = brewer.pal(8, "Accent")[5], 
            "Reverse Primer" = "#0B4D46", 
            "Covered Sequence" = "black", 
            "Not-covered Sequence" = "grey70")
    m <- match(unique(d$Type), names(col))
    col <- col[m[which(!is.na(m))]]
    label.idx <- which(d$Type %in% c("Covered Sequence", "Not-covered Sequence"))
    label.idx.primer <- which(d$Type %in% c("Forward Primer", "Reverse Primer") & 
        !duplicated(d$y))
    label.idx <- c(label.idx, label.idx.primer)
    label.order <- order(d$y[label.idx])
    labels <- d$ID[label.idx][label.order]
    template.data <- d[d$Type %in% c("Covered Sequence", "Not-covered Sequence"), ]
    primer.data <- d[d$Type %in% c("Forward Primer", "Reverse Primer"), ]
    if (length(primer.data) == 0 || nrow(primer.data) == 0) {
        # no primers to plot
        return(NULL)
    }
    break.step.size <- 10
    lim.min <- min(d$x_start, na.rm = TRUE)
    lim.max <- max(d$x_end, na.rm = TRUE)
    ticks <- c(seq(lim.min, lim.max, break.step.size))
    # determine coverage count
    N.total <- length(which(!d$Type %in% c("Forward Primer", "Reverse Primer")))
    N.covered <- length(which(d$Type == "Covered Sequence"))
    N.p <- N.covered/N.total
    stat.text <- paste(N.covered, " of ", N.total, " (", round(N.p, 3) * 100, "%) seqs covered", 
        sep = "")
    title <- stat.text
    segment.size <- 1.5
    arrow.size <- 2
    # arrow angle needs to be a function of the number of coverage events -> the more events, the larger the angle: from 5 to 30
    # number of primers
    min.arrow.angle <- 30 # was 5 before, not necessary anymore
    max.arrow.angle <- 30
    arrow.angle <- max(min(0.5 * (nrow(primer.data) + nrow(template.data)), max.arrow.angle), min.arrow.angle)
    # arrow length should be scaled according to the x-extent of the plot
    x.extent <- max(template.data$x_end) - min(template.data$x_start)
    p.extent <- max(primer.data$x_end) - min(primer.data$x_start)
    arrow.length <- 80 * (p.extent / x.extent^2)
    x.ticks <- c(seq(floor(lim.min/10) * 10, floor(lim.max/10) * 10, break.step.size))
    x.labels <- x.ticks
    x.labels[x.labels > 0] <- paste0("+", x.labels[x.labels > 0])
    r.colors <- c("#e5f4ff", "#ffefe5")
    ggplot() + 
        # need show.legend = FALSE for the segment, otherwise arrow shows
        geom_segment(data = primer.data, show.legend = FALSE,
            lineend = "butt", size = arrow.size,  
            aes_string(x = "x_start", xend = "x_end", y = "y", 
                       yend = "y", colour = "Type"),
            arrow = arrow(angle = arrow.angle, ends = "last", type = "closed")) + 
             geom_segment(data = template.data, 
                     aes_string(x = "x_start", xend = "x_end", 
                                y = "y", yend = "y", colour = "Type"), 
                    lineend = "square", size = segment.size) + 

        # x-axis rectangles to annotate binding/amplification region:
        geom_rect(data = region.df, 
            mapping = aes_string(xmin="xmin", xmax="xmax", 
                        ymin="ymin", ymax="ymax"), 
            fill = r.colors, alpha = 0.5,
            colour = "#3d3835", 
            size = segment.size) +
        # text for region annotation
        geom_text(data=region.df, 
                  aes_string(x = "xmin+(xmax-xmin)/2", 
                             y = "ymin+(ymax-ymin)/2", 
                             label = "Region"), 
                  size = 4) +
        ggtitle(title) + 
        xlab("Relative primer position") + 
        ylab("Sequence") + 
        scale_x_continuous(limits = c(lim.min, lim.max), 
                            breaks = x.ticks,
                            labels = x.labels) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        geom_vline(data = region.df, aes_string(xintercept = "RelStartPosition"), colour = "red") + # start of binding region
        geom_vline(xintercept = -1, colour = "red") + 
        #scale_y_continuous(limits = c(1, max(d$y, na.rm = TRUE)), 
        scale_y_continuous(
                breaks = 1:length(labels), labels = abbreviate(labels, getOption("openPrimeR.plot_abbrev"))) + 
        scale_colour_manual(values = col) +
        theme(legend.title = element_text(face = "bold"),
              legend.position = "top")
        # arrow arg doesn't work for linetype ...
        #guides(colour = guide_legend(override.aes = 
                #list(size = 1, linetype = "solid", arrow = arrow(length = unit(1, "cm")))))
       
}

#' Plot of Excluded Primers
#' 
#' Plots histogram of excluded primers.
#'
#' @param excluded.df Data frame with excluded primers.
#' @param filtered.stats Data frame with statistics of the filtering procedure.
#' @param template.df Template data frame.
#' @return A plot of excluded primers.
#' @keywords internal
plot.excluded.hist <- function(excluded.df, filtered.stats, template.df) {
    # histogram of excluded sequences
    if (length(excluded.df) == 0 || nrow(excluded.df) == 0) {
        # nothing to plot
        return(NULL)
    }
    if (!"Coverage_Ratio" %in% excluded.df || all(is.na(excluded.df$Coverage_Ratio))) {
        # no cvg to plot
        return(NULL)
    }
    ylab <- "Coverage"
    xlab <- "Constraint"
    tit <- "Excluded primer coverage"
    plot.df <- excluded.df
    plot.df$Filter_Reason <- factor(plot.df$Exclusion_Reason, levels = unique(filtered.stats$Constraint))
    ggplot(plot.df, aes_string(x = "Filter_Reason", y = "Coverage_Ratio", colour = "Coverage_Ratio")) + 
        ylab(ylab) + xlab(xlab) + ggtitle(tit) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        scale_colour_gradient(name = "Coverage") + 
        geom_boxplot(outlier.colour = NA) + 
        geom_point(position = position_jitter(width = 0.5, height = 0)) +
        scale_y_continuous(labels = scales::percent) +
        facet_wrap(~Direction)
}
#' Bar Plot of Template Coverage.
#'
#' Creates a bar plot showing the coverage for every group of template sequences.
#'
#' @param primers Either a \code{Primers} object with evaluated primer coverage
#' or a list containing \code{Primers} objects.
#' @param templates If \code{primers} is a \code{Primers} object, \code{templates} should be a \code{Templates} object.
#' If \code{primers} is a list, then \code{templates} should be a list of \code{Templates} objects.
#' @param per.mismatch A logical specifying whether the visualization should be stratified
#' according to the allowed number of mismatches. By default,
#' \code{per.mismatch} is set to \code{FALSE} such that the overall coverage
#' is plotted.
#' @param ... Optional arguments \code{groups} (a character vector of groups to be plotted when \code{primers} is a single primer set), \code{highlight.set} (the identifier of a primer set to be highlighted when \code{primers} is a list)
#'
#' @return A plot showing the number of covered template sequences.
#' @family templates
#' @export
#' @include primers.R templates.R
#' @family coverage visualizations
#' @examples
#' # Visualize the template coverage of a single primer set
#' data(Ippolito)
#' plot_template_cvg(primer.df, template.df)
#' # Stratify by allowed mismatches:
#' plot_template_cvg(primer.df, template.df, per.mismatch = TRUE)
#' # Compare the coverage of multiple primer sets
#' data(Comparison)
#' plot_template_cvg(primer.data[1:3], template.data[1:3])
#' # Stratify by allowed mismatches:
#' plot_template_cvg(primer.data[1:3], template.data[1:3], per.mismatch = TRUE)
setGeneric("plot_template_cvg", 
    function(primers, templates, per.mismatch = FALSE, ...) {
        standardGeneric("plot_template_cvg")
})

#' Bar Plot of Template Coverage.
#'
#' Creates a bar plot showing the coverage for every group of template sequences.
#'
#' @param primers A \code{Primers} object with evaluated primer coverage.
#' @param templates A \code{Templates} object containing the template sequences.
#' @param per.mismatch Whether to stratify by mismatches.
#' @param groups Identifiers of template groups for which plot should be created. By default, \code{groups} is set to \code{NULL} such that all 
#' templates are considered.
#' according to the number of mismatches between primer-template pairs.
#' @return A plot showing the number of covered template sequences.
#' @keywords internal
setMethod("plot_template_cvg", 
    methods::signature(primers = "Primers", templates = "Templates"),
    function(primers, templates, per.mismatch, groups = NULL) {
    if (length(primers) == 0 || length(templates) == 0) {
        return(NULL)
    } 
    if (!is(primers, "Primers")) {
        stop("Please input a valid primer data frame.")
    }
    if (!is(templates, "Templates")) {
        stop("Please input a valid template data frame.")
    }
    if (per.mismatch) {
        # stratify by mismatches
        p <- plot_template_cvg_mismatches(primers, templates)
    } else {
        p <- plot_template_cvg_unstratified(primers, templates)
    }
    return(p)
})

#' Bar Plot of Template Coverage.
#'
#' Creates a bar plot showing the coverage for every group of template sequences.
#'
#' @param primers A \code{Primers} object with evaluated primer coverage.
#' @param templates A \code{Templates} object containing the template sequences.
#' @param groups Identifiers of template groups for which plot should be created. By default, \code{groups} is set to \code{NULL} such that all 
#' templates are considered.
#' according to the number of mismatches between primer-template pairs.
#' @return A plot showing the number of covered template sequences.
#' @keywords internal
plot_template_cvg_unstratified <- function(primers, templates, groups = NULL) {
    # get statistics on expected coverage
    stats <- get_cvg_stats(primers, templates)
    if (length(stats) == 0) {
        warning("No coverage statistics available for the input primers.")
        return(NULL)
    }
    cvg.con <- paste0(round(stats[stats$Group == "Total", "Coverage_Ratio"] * 100, 0), "%")
    # compute coverage stats according to text identity:
    stats.txt <- get_cvg_stats(primers, templates, allowed.mismatches = 0, cvg.definition = "basic")
    cvg.txt <- paste0(round(stats.txt[stats.txt$Group == "Total", "Coverage_Ratio"] * 100, 0), "%")
    # need to melt
    # change names of columns
    vars <- c("Group", "primer_coverage", "N")
    p.df.con <- reshape2::melt(stats[, vars], c("Group"), variable.name = "Status", value.name = "Count")
    levels(p.df.con$Status) <- c("Expected Coverage (E)", "Available Templates")
    vars <- c("Group", "primer_coverage")
    p.df.txt <- reshape2::melt(stats.txt[, vars], c("Group"), variable.name = "Status", value.name = "Count")
    p.df.txt$Status <- "Identity Coverage (I)"
    # integrate both coverage results
    plot.df <- rbind(p.df.txt, p.df.con)
    plot.df$Status <- factor(plot.df$Status, levels = unique(plot.df$Status))
    # remove total cvg column
    plot.df <- plot.df[plot.df$Group != "Total",]
    # select groups
    if (length(groups) != 0 && !"all" %in% groups) {
        plot.df <- plot.df[plot.df$Group %in% groups, ]
    }
    colnames(plot.df)[colnames(plot.df) == "N"] <- "Available Templates"
    xlab <- "Group"
    ylab <- "Number of templates"
    title <- paste0("Covered templates (I = ", cvg.txt, ", E = ", cvg.con, ")")
    cols.1 <- brewer.pal(8, "Accent")[c(5, 8)]
    cols.2 <- brewer.pal(8, "Blues")[3]
    colors <- c(cols.2, cols.1)
    names(colors) <- levels(plot.df$Status)
    p <- ggplot(plot.df) + geom_bar(aes_string(x = "Group", y = "Count", fill = "Status"), stat = "identity", 
        position = "dodge") + xlab(xlab) + ggtitle(title) + 
        ylab(ylab) + 
        theme(axis.text.x = element_text(
                angle = 90, 
                hjust = 1, vjust = 0.5)) +
        scale_fill_manual(values = colors)
    return(p)
}
#' Preparation of Data for Plotting Mismatch Template Coverage.
#'
#' Creates a data frame for plotting a bar plot for the covered templates per allowed mismatches.
#'
#' @param primer.df A \code{Primers} object.
#' @param template.df A \code{Templates} object.
#' @param allowed.mismatches An optional numeric specifying the number of mismatches to be considered at most for plotting.
#' If not provided, the maximal number of mismatches found for the input primer 
#' set is used.
#' @return A data frame for creating a plot.
#' @keywords internal
prepare_template_cvg_mm_data <- function(primer.df, template.df, allowed.mismatches = NULL) {
    plot.df <- get_template_cvg_data(primer.df, template.df)
    # get only the unique events for every template (for both: fw & primer are found here)
    plot.df <- plyr::ddply(plot.df, c("Template", "Status"), function(x) plyr::arrange(x, substitute(Number_of_mismatches))[1,])
    # get data for available templates for every mismatch category:
    t.df <- data.frame("Template" = unique(template.df$ID), 
                       "Group" = NA, "Position" = NA, "Number_of_mismatches" = NA)
    m <- match(t.df$Template, template.df$ID)
    t.df$Group <- template.df$Group[m]
    if (length(allowed.mismatches) == 0) {
        allowed.mismatches <- max(plot.df$Number_of_mismatches)
    }
    mm.settings <- seq(0, allowed.mismatches)
    mm.settings <- mm.settings[order(mm.settings)]
    # identify available number of template sequences
    available.data <- vector("list", length(mm.settings))
    for (i in seq_along(mm.settings)) {
        df.a <- t.df
        df.a$Number_of_mismatches <- mm.settings[i]
        available.data[[i]] <- df.a
    }
    available.df <- do.call(rbind, available.data)
    available.df$Status <- "Available"
    p.df <- rbind(plot.df[, !colnames(plot.df) %in% c("Primer", "Direction")], available.df)
    cvg.groups <- c("Coverage", "Basic Coverage")
    all.groups <- c(cvg.groups, "Available Templates")
    levels(p.df$Status) <- all.groups
    # do not consider the exact nbr of mismatches, but the cumulative number (cvg with >= x mismatches rather than cvg with == x mismatches)
    cum.data <- p.df[p.df$Status == "Available Templates", ] # start from available data
    cum.data$Maximal_mismatches <- cum.data$Number_of_mismatches
    p.df$Number_of_mismatches[p.df$Status == "Basic Coverage"]
    for (i in seq_along(mm.settings)) {
        mm <- mm.settings[i]
        for (j in seq_along(cvg.groups)) {
            cvg.status <- cvg.groups[j]
            event.idx <- which(p.df$Status == cvg.status & p.df$Number_of_mismatches <= mm)
            if (length(event.idx) != 0) {
                cur.df <- p.df[event.idx,]
                cur.df$Maximal_mismatches <- mm
                cum.data <- rbind(cum.data, cur.df) 
            }
        }
    }
    o <- order(unique(as.character(cum.data$Group)))
    cum.data$Group <- factor(cum.data$Group, levels = unique(as.character(cum.data$Group))[o])
    # need to plot using identity as stat in order to ensure that bars have all the same width
    plot.df <- plyr::ddply(cum.data, c("Group", "Maximal_mismatches", "Status"), plyr::summarize, Count = length(substitute(Group)))
    additional.df <- expand.grid(Maximal_mismatches = mm.settings, Group = unique(template.df$Group), Status = unique(plot.df$Status), Count = 0)
    plot.df <- merge(plot.df, additional.df, all = TRUE)
    # for the added, duplicate events, select the 'real events'
    plot.df <- plyr::ddply(plot.df, c("Group", "Maximal_mismatches", "Status"), plyr::summarise,
                           Count = max(substitute(Count)))
    total.percentages <- TRUE
    if (!total.percentages) {
        idx <- which(plot.df$Status == "Available Templates")
        idx <- idx[!duplicated(plot.df$Group[idx])]
        available.per.group <- plot.df$Count[idx]
        names(available.per.group) <- plot.df$Group[idx]
        m <- match(plot.df$Group, names(available.per.group))
        plot.df$Coverage_Ratio <- plot.df$Count / available.per.group[m]
    } else {
        plot.df$Coverage_Ratio <- plot.df$Count / nrow(template.df)
    }
    return(plot.df)
}
#' Bar Plot of Template Coverage for Mismatches.
#'
#' Creates a bar plot showing the coverage for every group of template sequences.
#'
#' @param primer.df A \code{Primers} object.
#' @param template.df A \code{Templates} object.
#' @return A plot showing the number of covered template sequences.
#' @keywords internal
plot_template_cvg_mismatches <- function(primer.df, template.df) {
    plot.df <- prepare_template_cvg_mm_data(primer.df, template.df)
    # compute cvg ratio per mismatch setting to show in facet labels:
    cvg.per.mm <- plyr::ddply(plot.df, c("Maximal_mismatches", "Status"), plyr::here(plyr::summarise),
                           Coverage_Ratio = sum(substitute(Count)) / nrow(template.df)) 
    cvg.per.mm <- cvg.per.mm[cvg.per.mm$Status == "Coverage",]
    cvg.info <- paste0("(", round(cvg.per.mm$Coverage_Ratio * 100, 0), "% coverage)")
    plot.df$Maximal_mismatches <- paste0(plot.df$Maximal_mismatches, " ", cvg.info[match(plot.df$Maximal_mismatches, cvg.per.mm$Maximal_mismatches)])
    xlab <- "Group"
    ylab <- "Number of templates"
    title <- "Covered templates: mismatches"
    cols.1 <- brewer.pal(8, "Accent")[c(5, 8)]
    cols.2 <- brewer.pal(8, "Blues")[3]
    colors <- c(cols.1[1], cols.2, cols.1[2])
    names(colors) <- levels(plot.df$Status)
    p <- ggplot(plot.df) + 
        geom_bar(aes_string(x = "Group", y = "Count", fill = "Status"),
                 position = "dodge", stat = "identity") + 
        xlab(xlab) + ggtitle(title) + 
        ylab(ylab) + 
        theme(axis.text.x = element_text(
                angle = 90, 
                hjust = 1, vjust = 0.5)) +
        scale_fill_manual(values = colors) + 
        facet_wrap(~Maximal_mismatches, ncol = 2,
        labeller = label_bquote("Mismatches"<=.(substitute(Maximal_mismatches))))
    return(p)
}
#' Templates Coverage for Multiple Primer Sets.
#'
#' Plots the coverage of multiple primer sets.
#'
#' @param primers List with primer data frames.
#' @param templates List with template data frames.
#' @param colors Color for every primer set.
#' @param highlight.set Primer sets to be highlighted.
#' @return A plot for comparing primer coverage.
#' @keywords internal
setMethod("plot_template_cvg", 
    methods::signature(primers = "list", templates = "list"),
    function(primers, templates, per.mismatch, highlight.set = NULL) {
    if (length(primers) == 0 || length(templates) == 0) {
        return(NULL)
    }
    if (per.mismatch) {
        p <- plot_template_cvg_comparison_mismatch(primers, templates, highlight.set)
    } else {
        p <- plot_template_cvg_comparison_unstratified(primers, templates, highlight.set)
    }
    return(p)
})
#' Getter for Run Names.
#'
#' Retrieves the run names of the input data.
#'
#' @param primer.data A list with \code{Primers} or \code{Templates}.
#' @return A vector with identifiers for every set.
#' @keywords internal
get.run.names <- function(primer.data) {
    # if primer.data is a named list, use the names of the list
    if (length(primer.data) == 0) {
        return(NULL)
    }
    if (length(names(primer.data) != 0)) {
        run.names <- names(primer.data)
    } else {
        # primer.data is not a named list -> use the 'Run' identifier instead
        #run.names <- unname(unlist(lapply(primer.data, function(x) unique(x$Run))))
        run.names <- lapply(primer.data, function(x) unique(x$Run))
        print(lapply(primer.data, function(x) nrow(x) == 0))
        empty.idx <- which(unlist(lapply(primer.data, function(x) nrow(x) == 0)))
        run.names[empty.idx] <- "Unknown"
        run.names <- unname(unlist(run.names))
    }
    return(run.names)
}
#' Templates Coverage for Multiple Primer Sets.
#'
#' Plots the coverage of multiple primer sets.
#'
#' @param primers List with primer data frames.
#' @param templates List with template data frames.
#' @param highlight.set Primer sets to be highlighted.
#' @return A plot for comparing primer coverage.
#' @keywords internal
plot_template_cvg_comparison_unstratified <- function(primers, templates, highlight.set = NULL) {

    if (length(primers) == 0 || length(templates) == 0) {
        return(NULL)
    }
    run.names <- get.run.names(primers)
    cvg <- lapply(seq_along(primers), function(x) {
        stats <-  openPrimeR::get_cvg_stats(primers[[x]], 
                    templates[[x]],
                    total.percentages = TRUE)
        if (length(stats) != 0) {
            stats <- cbind("Run" = run.names[x], stats)
        }
    })
    plot.df <- do.call(plyr::rbind.fill, cvg)
    if (length(plot.df) == 0 || nrow(plot.df) == 0) {
        return(NULL)
    }
    # remove the total group
    plot.df <- plot.df[plot.df$Group != "Total", ]
    plot.df$Group <- factor(plot.df$Group, levels = unique(plot.df$Group)[order(unique(plot.df$Group))])
    title <- "Coverage of the template sequences"
    plot.df$Run <- abbreviate(plot.df$Run, getOption("openPrimeR.plot_abbrev"))
    plot.df$Run <- factor(plot.df$Run, levels = unique(plot.df$Run)[order(as.character(unique(plot.df$Run)))])
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
    pal <- getOption("openPrimeR.plot_colors")["Group"] # the RColorBrewer palette to use
    colors <- colorRampPalette(brewer.pal(8, pal))(length(unique(plot.df$Group)))
    # add label for top of the bars (the first group) 
    labels <- rep(NA, nrow(plot.df))
    for (i in seq_along(levels(plot.df$Run))) {
        run <- levels(plot.df$Run)[i]
        idx <- which(plot.df$Run == run)
        cur.data <- plot.df[idx,]
        res <- sum(cur.data$Coverage_Ratio) # the coverage to show for this bar
        out.idx <- which(plot.df$Run == run & plot.df$Group == unique(cur.data$Group)[1]) # row index to show the coverage
        labels[out.idx] <- res
    }
    plot.df$Label <- ifelse(is.na(labels), "", paste0(round(labels * 100, 1), "%"))

    p <- ggplot(plot.df, aes_string(x = "Run", y = "Coverage_Ratio", fill = "Group")) + 
        geom_bar(stat = "identity") + 
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
        scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
        ggtitle(title) +
        ylab("Coverage") +
        scale_fill_manual(values = colors) + 
        # add coverage above bars:
        geom_text(aes(label = substitute(Label)),
                    size = 5, 
                    position = position_stack(vjust = 0.5),
                    check_overlap = FALSE)
    if (length(highlight.set) != 0) {
        # highlight selected sets
        sel <- levels(plot.df$Run) %in% highlight.set
        p <- p + theme(
                axis.text.x = element_text(
                    face=ifelse(sel, "bold","plain"),
                    colour = ifelse(sel, 
                            "grey20", "grey30")))
    }
    return(p)
}
#' Templates Coverage for Multiple Primer Sets.
#'
#' Plots the coverage of multiple primer sets.
#'
#' @param primers List with primer data frames.
#' @param templates List with template data frames.
#' @param highlight.set Primer sets to be highlighted.
#' @param show.percentages Whether to show the total coverage percentage above the bars.
#' @return A plot for comparing primer coverage.
#' @keywords internal
plot_template_cvg_comparison_mismatch <- function(primers, templates, 
            highlight.set = NULL, show.percentages = FALSE) {

    if (length(primers) == 0 || length(templates) == 0) {
        return(NULL)
    }
    run.names <- get.run.names(primers)
    max.mm <- max(unlist(lapply(primers, function(x) as.numeric(strsplit(c(x$Nbr_of_mismatches_fw, x$Nbr_of_mismatches_rev), split = ",")[[1]]))))
    cvg <- parallel::mclapply(seq_along(primers), function(x) {
        # supply allowed mismatch arg to show the same number of mismatches for all sets in the plot's facets
        plot.df <- prepare_template_cvg_mm_data(primers[[x]], templates[[x]],
                    allowed.mismatches = max.mm)
        if (length(plot.df) != 0) {
            # only select the expected coverage events
            plot.df <- plot.df[plot.df$Status == "Coverage",]
            # annotate with Run
            plot.df <- cbind("Run" = run.names[x], plot.df)
        }
    })
    plot.df <- do.call(plyr::rbind.fill, cvg)
    if (length(plot.df) == 0 || nrow(plot.df) == 0) {
        return(NULL)
    }
    title <- "Coverage of the template sequences"
    plot.df$Run <- abbreviate(plot.df$Run, getOption("openPrimeR.plot_abbrev"))
    plot.df$Run <- factor(plot.df$Run, levels = unique(plot.df$Run)[order(as.character(unique(plot.df$Run)))])
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
    pal <- getOption("openPrimeR.plot_colors")["Group"] # the RColorBrewer palette to use
    colors <- colorRampPalette(brewer.pal(8, pal))(length(unique(plot.df$Group)))
    # add label for top of the bars (the first group) 
    labels <- rep(NA, nrow(plot.df))
    for (i in seq_along(levels(plot.df$Run))) {
        run <- levels(plot.df$Run)[i]
        idx <- which(plot.df$Run == run)
        cur.data <- plot.df[idx,]
        res <- plyr::ddply(cur.data, c("Maximal_mismatches"), plyr::summarize, Label = sum(substitute(Coverage_Ratio)))$Label
        out.idx <- which(plot.df$Run == run & plot.df$Group == levels(plot.df$Group)[1])
        labels[out.idx] <- res
    }
    plot.df$Label <- ifelse(is.na(labels), "", paste0(round(labels * 100, 1), "%"))
    p <- ggplot(plot.df, aes_string(x = "Run", y = "Coverage_Ratio", fill = "Group", group = "Group")) + 
        geom_bar(stat = "identity") + 
        facet_wrap(~Maximal_mismatches,
            labeller = label_bquote("Mismatches"<=.(substitute(Maximal_mismatches)))
        ) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
        ggtitle(title) +
        ylab("Coverage") +
        scale_fill_manual(values = colors)
        # add coverage above bars:
    if (show.percentages) {
        p <- p + geom_text(aes(label = substitute(Label)),
                    size = 3, 
                    position = position_stack(vjust = 0.5),
                    check_overlap = FALSE)
    }
    if (length(highlight.set) != 0) {
        # highlight selected sets
        sel <- levels(plot.df$Run) %in% highlight.set
        p <- p + theme(
                axis.text.x = element_text(
                    face=ifelse(sel, "bold","plain"),
                    colour = ifelse(sel, "grey20", "grey30")))
    }
    return(p)
}
#' Coverage Ratios per Group of Templates.
#'
#' Retrieve statistics on covered templates, either for a single primer set
#' or for multiple primer sets.
#'
#' @param primers To retrieve coverage statistics for a single primer set, 
#' please provide an object of class \code{Primers} containing primers with evaluated coverage.
#' To retrieve coverage statistics for multiple primer sets, pelase provide
#' a list with evaluated \code{Primers} objects.
#' @param templates If \code{primers} is an object of class \code{Primers},
#' please provide an object of class \code{Templates} containing the
#' template sequences targeted by \code{primers}. If \code{primers} is a list,
#' \code{templates} should be a list of \code{Template} objects.
#' @param for.viewing Whether the table should be formatted
#' to be human-readable. By default, \code{for.viewing} is \code{FALSE}.
#' @param total.percentages Whether group coverage percentages
#' should be computed in relation to the total number of template sequences
#' or in relation to the number of templates belonging to a specific group.
#' By default, \code{total.percentages} is \code{FALSE} suc that the
#' percentages are group-specific.
#' @param allowed.mismatches The maximal allowed number of mismatches.
#' By default, \code{allowed.mismatches} is set to \code{Inf} such that the number of mismatches is not restricted additionally.
#' @param cvg.definition If \code{cvg.definition} is set to
#' "constrained", the statistics for the expected
#' coverage (after applying the coverage constraints) are retrieved.
#' If \code{cvg.definition} is set to "basic", the coverage is determined 
#' solely by string matching (i.e. without applying the coverage constraints).
#' By default, \code{cvg.definition} is set to "constrained".
#' @return Data frame whose entries provide the coverage of templates belonging to a specific group.
#' @export
#' @include primers.R templates.R
#' @examples
#' # Coverage statistics for a single primer set
#' data(Ippolito)
#' cvg.stats <- get_cvg_stats(primer.df, template.df)
#' # Coverage statistics for multiple primer sets
#' data(Comparison)
#' cvg.stats.comp <- get_cvg_stats(primer.data, template.data)
setGeneric("get_cvg_stats", 
    function(primers, templates, for.viewing = FALSE, total.percentages = FALSE, 
            allowed.mismatches = Inf, cvg.definition = c("constrained", "basic")) {
        standardGeneric("get_cvg_stats")
})
#' Coverage Statistics of a Primer Set.
#'
#' Retrieve statistics on the templates that are covered by a primer set.
#'
#' @param primer.df An object of class \code{Primers} containing
#' primers with evaluated coverage.
#' @param template.df An object of class \code{Templates} containing
#' templates with evaluated coverage.
#' @param for.viewing Whether the table should be formatted
#' for viewing rather than processing.
#' @param total.percentages Whether group coverage percentages
#' should relate to all template sequences or just those templates
#' belonging to a specific group.
#' @param allowed.mismatches The maximal allowed number of mismatches.
#' By default, the number of mismatches is not restricted. 
#' @param cvg.definition If \code{cvg.definition} is set to
#' "constrained", the statistics for the expected
#' coverage (after applying the coverage constraints) are retrieved.
#' If \code{cvg.definition} is set to "basic", the coverage is determined 
#' solely by string matching (i.e. without applying the coverage constraints).
#' By default, \code{cvg.definition} is set to "constrained".
#' @return Data frame with coverage statistics.
#' @keywords internal
setMethod("get_cvg_stats", methods::signature(primers = "Primers"), 
    function(primers, templates, for.viewing = FALSE, total.percentages = FALSE, 
             allowed.mismatches = Inf, cvg.definition = c("constrained", "basic")) {

    # compute some stats on the coverage table structure: group | nbr templates |
    # [nbr covered | covered percent |]
    if (!is(primers, "Primers")) {
        stop("Please input a valid primer data frame.")
    }
    if (!is(templates, "Templates")) {
        stop("Please input a valid template data frame.")
    }
    if (!"Covered_Seqs" %in% colnames(primers) || !"primer_coverage" %in% colnames(primers)) {
        warning("get_cvg_stats: Primer/template coverage not available.")
        return(NULL)  # not possible
    }
    cvg.definition <- match.arg(cvg.definition)
    # update primer & template coverage according to the input conditions:
    primers <- update_primer_cvg(primers, templates, allowed.mismatches = allowed.mismatches, cvg.definition = cvg.definition)
    templates <- update_template_cvg(templates, primers)
    if (total.percentages) {
        # normalize coverage by total number of templates
        df <- ddply(templates, "Group", plyr::here(summarize), N = nrow(templates),
            N_primer = length(unique(strsplit(as.character(substitute(Covered_By_Primers)), split = ",")[[1]])), 
            primer_coverage = length(which(substitute(primer_coverage) > 0)), 
            Coverage_Ratio = sum(substitute(primer_coverage))/nrow(templates),
            N_primer_fw = length(unique(strsplit(as.character(substitute(Covered_By_Primers_fw)), split = ",")[[1]])), 
            primer_coverage_fw = length(which(substitute(primer_coverage_fw) > 0)), 
            Coverage_Ratio_fw = sum(substitute(primer_coverage_fw))/nrow(templates),
            N_primer_rev = length(unique(strsplit(as.character(substitute(Covered_By_Primers_rev)), split = ",")[[1]])), 
            primer_coverage_rev = length(which(substitute(primer_coverage_rev) > 0)), 
            Coverage_Ratio_rev = sum(substitute(primer_coverage_rev))/nrow(templates))
    } else {
        # normalize coverage by number of templates per group
        df <- ddply(templates, "Group", plyr::here(summarize), N = length(substitute(Group)), N_primer = length(unique(strsplit(as.character(substitute(Covered_By_Primers)), 
            split = ",")[[1]])), primer_coverage = length(which(substitute(primer_coverage) > 0)), 
            Coverage_Ratio = sum(substitute(primer_coverage))/length(substitute(Group)), N_primer_fw = length(unique(strsplit(as.character(substitute(Covered_By_Primers_fw)), 
                split = ",")[[1]])), primer_coverage_fw = length(which(substitute(primer_coverage_fw) > 
                0)), Coverage_Ratio_fw = sum(substitute(primer_coverage_fw))/length(substitute(Group)), 
            N_primer_rev = length(unique(strsplit(as.character(substitute(Covered_By_Primers_rev)), 
                split = ",")[[1]])), primer_coverage_rev = length(which(substitute(primer_coverage_rev) > 
                0)), Coverage_Ratio_rev = sum(substitute(primer_coverage_rev))/length(substitute(Group)))
    }
    # order by groups
    o <- order(df$Group)
    df <- df[o, ]
    # summary of all groups in 'Total' row
    df <- rbind(data.frame(Group = "Total", N = nrow(templates), N_primer = nrow(primers), 
        primer_coverage = sum(df$primer_coverage), Coverage_Ratio = sum(df$primer_coverage)/nrow(templates), 
        N_primer_fw = length(which(primers$Forward != "")), primer_coverage_fw = sum(df$primer_coverage_fw), 
        Coverage_Ratio_fw = sum(df$primer_coverage_fw)/nrow(templates), N_primer_rev = length(which(primers$Reverse != 
            "")), primer_coverage_rev = sum(df$primer_coverage_rev), Coverage_Ratio_rev = sum(df$primer_coverage_rev)/nrow(templates), 
        stringsAsFactors = FALSE), data.frame(df, stringsAsFactors = FALSE))
    df$Coverage <- paste(df$primer_coverage, " of ", df$N, " (", round(df$Coverage_Ratio * 
        100, 2), "%)", sep = "")
    df$Coverage_fw <- paste(df$primer_coverage_fw, " of ", df$N, " (", round(df$Coverage_Ratio_fw * 
        100, 2), "%)", sep = "")
    df$Coverage_rev <- paste(df$primer_coverage_rev, " of ", df$N, " (", round(df$Coverage_Ratio_rev * 
        100, 2), "%)", sep = "")
    out.df <- df
    if (for.viewing) {
        cols.both <- c("Group", "Coverage", 
                "Coverage_fw", "Coverage_rev")
        out.names <-  c("Group", "Coverage", "Coverage (fw)", "Coverage (rev)")
        mode <- "both"
        if (sum(out.df$N_primer_rev) == 0 || sum(out.df$N_primer_fw) == 0) {
            mode <- "single"
        }
        if (mode == "single") {
            out.df <- out.df[, cols.both[1:2]]
            colnames(out.df) <- out.names[1:2]
        } else {
            out.df <- out.df[, cols.both]
            colnames(out.df) <- out.names
        }
    }
    return(out.df)
})
#' Coverage Statistics for Multiple Primer Sets.
#'
#' Retrieve statistics on covered templates for multiple primer sets.
#'
#' @param primers A list with objects of class \code{Primers} containing
#' primers with evaluated coverage.
#' @param templates A list with objects of class \code{Templates} containing
#' templates with evaluated coverage.
#' @param for.viewing Whether the table should be formatted
#' for viewing rather than processing.
#' @param total.percentages Whether group coverage percentages
#' should relate to all template sequences or just those templates
#' belonging to a specific group.
#' @param allowed.mismatches The maximal allowed number of mismatches.
#' By default, the number of mismatches is not restricted. 
#' @param cvg.definition If \code{cvg.definition} is set to
#' "constrained", the statistics for the expected
#' coverage (after applying the coverage constraints) are retrieved.
#' If \code{cvg.definition} is set to "basic", the coverage is determined 
#' solely by string matching (i.e. without applying the coverage constraints).
#' By default, \code{cvg.definition} is set to "constrained".
#' @return Data frame with coverage statistics.
#' @keywords internal
setMethod("get_cvg_stats", methods::signature(primers = "list"),
    function(primers, templates, for.viewing = FALSE, total.percentages = FALSE,
             allowed.mismatches = Inf, cvg.definition = c("constrained", "basic")) {


    if (length(primers) == 0 || length(templates) == 0) {
        return(NULL)
    }
    template.classes <- sapply(templates, function(x) class(x))
    primer.classes <- sapply(primers, function(x) class(x))
    if (any(template.classes != "Templates") || any(primer.classes != "Primers")) {
        stop("Check types of primers/templates.")
    }
    runs <- get.run.names(primers)
    stats <- lapply(seq_along(primers), function(x)  {
                df <- get_cvg_stats(primers[[x]], templates[[x]], FALSE, total.percentages) 
                groups <- as.character(df$Group)
                if (for.viewing) {
                    # only select the coverage entry
                    cvg <- df[, which(colnames(df) == "Coverage_Ratio") ]
                    # format cvg
                    cvg <- paste0(round(cvg * 100, 1), "%")
                    df <- data.frame(Run = runs[x], Coverage = t(cvg), stringsAsFactors = FALSE)
                } else {
                    df <- data.frame(cbind(Run = runs[x], t(df)), stringsAsFactors = FALSE)
                }
                colnames(df) <- c("Run", groups)
                return(df)
            })
    stat.df <- do.call(plyr::rbind.fill, stats)
    if (ncol(stat.df) >= 3) {
        # don't change Run & Total columns
        idx <- seq(3, ncol(stat.df))
        group.order <- colnames(stat.df)[idx]
        group.order <- group.order[order(group.order)]
        colnames(stat.df)[idx] <- group.order
    }
    # order rows by 'Run':
    stat.df <- stat.df[order(stat.df$Run), ]
    return(stat.df)
})
#' Plot of Primer Coverage.
#'
#' Shows which groups of templates are covered by individual primers.
#' 
#' @param primers An object of class \code{Primers} or a list with with objects
#' of class \code{Primers}.
#' @param templates If \code{primers} is an object of class \code{Primers},
#' please supply a \code{Templates} object. If \code{primers} is a list,
#' please supply a corresponding list with \code{Templates} objects.
#' @param per.mismatch A logical identifiying whether the coverage should
#' be plotted for individual settings of allowed mismatches. By default
#' \code{per.mismatch} is set to \code{FALSE} such that the overall
#' coverage is plotted.
#' @param ... \code{groups} (the identifiers of template groups
#' to be excluded from the plot if \code{primers} is a single primer set)
#' @return A plot showing the coverage of individual primers.
#' @export
#' @family coverage visualizations
#' @include primers.R templates.R
#' @examples
#' # Plot expected coverage per primer
#' data(Ippolito)
#' plot_primer_cvg(primer.df, template.df)
#' # Plot coverage stratified by allowed mismatches:
#' plot_primer_cvg(primer.df, template.df, per.mismatch = TRUE)
#' # Plot coverage of multiple primer sets
#' data(Comparison)
#' plot_primer_cvg(primer.data[1:3], template.data[1:3])
setGeneric("plot_primer_cvg", 
    function(primers, templates, per.mismatch = FALSE, ...) {
        standardGeneric("plot_primer_cvg")
})
#' Plot Individual Primer Coverage.
#'
#' Shows which templates are covered by individual primers.
#'
#' @param p.df Primer data frame.
#' @param template.df Template data frame.
#' @param per.mismatch Whether to stratify by allowed mismatches.
#' @param excluded.seqs Identifiers of templates that should not be considered.
#' @param per.mismatch Whether the coverage should
#' be broken down for individual settings of allowed mismatches.
#' @return A bar plot showing the coverage of individual primers.
#' @keywords internal
setMethod("plot_primer_cvg", 
    methods::signature(primers = "Primers", templates = "Templates"),
    function(primers, templates, per.mismatch = FALSE, groups = NULL) {
    
    if (length(primers) == 0 || nrow(primers) == 0 || length(templates) == 0) {
        return(NULL)
    }
    if (!is(primers, "Primers")) {
        stop("Please input a valid primer data frame.")
    }
    if (!is(templates, "Templates")) {
        stop("Please input a valid template data frame.")
    }
    if (per.mismatch) {
        p <- plot_primer_cvg_mismatches(primers, templates)
    } else {
        p <- plot_primer_cvg_unstratified(primers, templates, groups = groups)
    }
    return(p)
})
#' Plot Multiple Primer Coverages.
#'
#' Plots the coverage of individual primers for multiple sets.
#'
#' @param primers List with \code{Primers} objects.
#' @param templates List with \code{Templates} objects.
#' @param per.mismatch Whether the coverage should
#' be broken down for individual settings of allowed mismatches.
#' @return A bar plot showing the coverage of individual primers.
#' @keywords internal
setMethod("plot_primer_cvg", 
    methods::signature(primers = "list", templates = "list"),
    function(primers, templates, per.mismatch = FALSE) {

    if (length(primers) == 0 || length(templates) == 0) {
        return(NULL)
    }
    if (per.mismatch) {
        warning("Not supported yet.")
        return(NULL)
    } else {
        # alias for existing exported function
        p <- plot_constraint(primers, NULL, "primer_coverage")
    }
    return(p)
})
#' Plot Individual Primer Coverage.
#'
#' Plots the coverage of individual primers.
#'
#' @param p.df Primer data frame.
#' @param template.df Template data frame.
#' @param groups Optional identifiers of template groups to be considered.
#' If not provided, all template groups are considered.
#' @return A bar plot showing the coverage of individual primers.
#' @keywords internal
plot_primer_cvg_unstratified <- function(p.df, template.df, groups = NULL) {
    # Select only the selected groups of templates:
     if (!is.null(groups) && !"all" %in% groups) { # select subset
        idx <- which(template.df$Group %in% groups)
        # select relevant templates
        template.df <- template.df[idx,]
        # set excluded seqs
        excluded.seqs <- setdiff(template.df$Identifier[seq_len(nrow(template.df))], template.df$Identifier[idx])
        # re-evalaute coverage with the new templates
        p.df <- evaluate.diff.primer.cvg(p.df, excluded.seqs, template.df)
    }
    title <- "Template coverage per primer"
    xlab <- "Primer Identifier"
    title <- "Template coverage per primer and group"
    # annotate p.df with strat column
    m <- covered.seqs.to.idx(p.df$Covered_Seqs, template.df)
    plot.df <- NULL
    strat <- "Group"
    # TODO: maybe modify the code inline with prepare_mm_plot?
    for (i in seq_along(m)) {
        missing.groups <- unique(template.df$Group)
        new.data <- NULL
        if (length(m[[i]]) != 0) {
            hit.groups <- unique(template.df[m[[i]], strat])
            new.data <- data.frame(ID = rep(p.df[i, "ID"], length(hit.groups)))
            new.data[, strat] <- hit.groups
            cvd.idx <- m[[i]]
            cvg.idx <- lapply(seq_along(hit.groups), function(x) {
              s <- hit.groups[x]
              sel <- which(template.df[cvd.idx, strat] == s)
              cvd.idx[sel]
            })
            cvg.counts <- sapply(cvg.idx, length)
            cvd.seqs <- sapply(cvg.idx, function(x) paste(template.df$Identifier[x], 
              collapse = ","))
            available.seqs.per.group <- sapply(hit.groups, function(x) length(which(template.df$Group == x)))
            new.data$primer_coverage <- cvg.counts/available.seqs.per.group # percentage of group templates covered
            new.data$Covered_Seqs <- cvd.seqs
            missing.groups <- setdiff(unique(template.df$Group), new.data$Group)
        }
        # also add all combinations for missing groups -> ensures bar width stays constant throughout the plot observations
        missing.df <- data.frame(ID = rep(p.df$ID[i], length(missing.groups)), Group = missing.groups, 
                                primer_coverage = rep(0, length(missing.groups)), Covered_Seqs = rep("", length(missing.groups)))
        new.data <- rbind(new.data, missing.df)
        plot.df <- rbind(plot.df, new.data)
    }
    ###
    ylab <- "Covered Templates"
    pal <- getOption("openPrimeR.plot_colors")[strat] # the RColorBrewer palette to use
    group.colors <- colorRampPalette(brewer.pal(8, pal))(length(unique(plot.df[, strat])))
    plot.df$ID <- factor(plot.df$ID)
    unique.cvd.idx <- compute.unique.covered.idx(plot.df, template.df)
    unique.cvg <- sapply(unique.cvd.idx, length)
    unique.cvg.ratio <- unique.cvg/nrow(template.df)
    p.df.unique <- plot.df
    p.df.unique$Coverage_Type <- "Unique Coverage"
    plot.df$Coverage_Type <- "Overall Coverage"
    p.df.unique$primer_coverage <- unique.cvg.ratio
    plot.df <- rbind(plot.df, p.df.unique)
    levels(plot.df$ID) <- abbreviate(levels(plot.df$ID), getOption("openPrimeR.plot_abbrev"))
    plot.df <- plot.df[plot.df$Coverage_Type == "Overall Coverage",] # only plot overall cvg
    bar.width <- 0.75
    dist.buffer <- (1- 0.75) / 2
    # add shading for every primer region in the plot to differentiate them better
    base.pos <- seq_along(unique(plot.df$ID))
    starts <- base.pos - 0.5 + dist.buffer
    ends <- base.pos + 0.5 - dist.buffer
    rects <- data.frame(xstart = starts, xend = ends,
                        ymin = 0, ymax = 1)
    rects$col <- factor(rep(c(0,1), nrow(rects))[seq_len(nrow(rects))])
    # only retain rectangles alternately
    rects <- rects[rects$col == 0, ]
    p <- ggplot() +
        geom_bar(data = plot.df, stat = "identity", position = "dodge", 
                    aes_string(x = "ID", y = "primer_coverage",
                    fill = strat), width = bar.width, show.legend = TRUE) + 
        geom_rect(data = rects, aes_string(xmin = "xstart", xmax = "xend", ymin = "ymin", ymax = "ymax"), 
                  fill = "grey60", alpha = 0.15, inherit.aes = FALSE, show.legend = FALSE) +
        xlab(xlab) + ggtitle(title) + ylab(ylab) + 
        theme(axis.text.x = element_text(angle = 60, 
                                hjust = 1)) +
        scale_fill_manual(values = group.colors) + 
        scale_y_continuous(limits = c(0, 1), labels = scales::percent)
    return(p)
}

#' Data for Mismatch Primer Coverage Plot.
#'
#' Ensures that there's an entry for every possible mismatch setting.
#'
#' @param primer.df A \code{Primers} object.
#' @param template.df A \code{Templates} object.
#' @return A data frame for plotting mismatch primer coverage.
#' @keywords internal
get_primer_cvg_mm_plot_df <- function(primer.df, template.df) {
    full.df <- prepare_mm_plot(primer.df, template.df)
    full.df <- full.df[full.df$Coverage_Type == "constrained",]
    df <- plyr::ddply(full.df, c("Primer", "Template", "Group"), plyr::summarize,
                            Position = unique(substitute(Position_3terminus)), 
                            Number_of_mismatches = unique(substitute(Number_of_mismatches)))
    # for every primer and group, ensure that there's an entry for every number of mismatches#
    # otherwise there would be missing data for individual facets (individual nbr of mismatches)
    count.df <- plyr::ddply(df, c("Primer", "Group", "Number_of_mismatches"), plyr::summarize,
                            Coverage = length(unique(substitute(Template))))
    if (nrow(df) == 0) {
        # just plot a single facet to show nothing is covered
        max.mm <- 0 
    } else {
        max.mm <- max(df$Number_of_mismatches)
    }
    result.data <- vector("list", max.mm + 1)
    for (i in seq(0, max.mm)) {
        additional.df <- data.frame(Primer = unlist(lapply(seq_len(nrow(primer.df)), function(x) rep(primer.df$ID[x], length(unique(template.df$Group))))),
                                    Group = unlist(lapply(seq_len(nrow(primer.df)), function(x) unique(template.df$Group))),
                                    Number_of_mismatches = i, Coverage = 0)
        count.df <- merge(count.df, additional.df, all = TRUE)
    }
    # ensure that for duplicate entries (0 coverage entries), the one with the maximal coverage is chosen:
    count.df <- plyr::ddply(count.df, c("Primer", "Group", "Number_of_mismatches"), plyr::summarise,
                            Coverage = max(substitute(Coverage)))
    # make the counts cumulative
    count.df$Cumulative_Coverage <- stats::ave(count.df$Coverage, count.df$Primer, count.df$Group, FUN = cumsum)
    # determine ratio of covered seqs per group
    group.counts <- table(template.df$Group) # match to count.df$Group for division
    m <- match(count.df$Group, names(group.counts))
    count.df$Coverage_Ratio <- count.df$Cumulative_Coverage / as.vector(group.counts[m])
    o <- order(unique(as.character(count.df$Group)))
    count.df$Group <- factor(count.df$Group, levels = unique(as.character(count.df$Group))[o])
    colnames(count.df)[colnames(count.df) == "Number_of_mismatches"] <- "Maximal_mismatches" 
    count.df$Primer <- factor(abbreviate(count.df$Primer, getOption("openPrimeR.plot_abbrev")), levels = abbreviate(levels(count.df$Primer), getOption("openPrimeR.plot_abbrev"))) # shorter identifiers
    return(count.df)
}
#' Plot of Individual Primer Coverage and Mismatches.
#'
#' Plots the coverage of individual primers for different mismatch settings.
#'
#' @param primer.df A \code{Primers} object.
#' @param template.df A \code{Templates} object.
#' @return A bar plot showing the coverage of individual primers for different mismatch settings.
#' @keywords internal
plot_primer_cvg_mismatches <- function(primer.df, template.df) {
    # retrieve cvg stats of each primer:
    count.df <- get_primer_cvg_mm_plot_df(primer.df, template.df)
    pal <- getOption("openPrimeR.plot_colors")["Group"] # the RColorBrewer palette to use
    group.colors <- colorRampPalette(brewer.pal(8, pal))(length(unique(count.df[, "Group"])))
    p <- ggplot(count.df) + 
        geom_bar(width = 0.75, position = "stack", stat = "identity",
            aes_string("Primer", 
                "Cumulative_Coverage",
                fill = "Group")) + 
        xlab("Primer") + ggtitle("Mismatch primer coverage") + 
        ylab("Number of covered templates") + 
        theme(axis.text.x = element_text(
                angle = 90,  # angle @ 90 to prevent overplotting when facetting
                hjust = 1)) +
        scale_fill_manual(values = group.colors) + 
        facet_wrap(~Maximal_mismatches,
            labeller = label_bquote("Mismatches"<=.(substitute(Maximal_mismatches))))
    return(p)
}

