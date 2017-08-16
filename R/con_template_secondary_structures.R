########
# Secondary structure of templates
########

#' Template Secondary Structures
#' 
#' Computes template secondary structures.
#'
#' @param template.df Template data frame.
#' @param annealing.temperature Temperature [C] at which 
#' to compute secondary structures.
#' @param regions List containing the positional intervals for which the template
#' secondary structure should be computed.
#' @param constraints String giving secondary structure constraints.
#' For example \code{xxxxxx...} would forbid folding in the first 6 bases
#' of a template with length 9 and allow folding in its last 3 bases.
#' @return Data frame with info on template secondary structures.
#' @keywords internal
compute.template.secondary.structures <- function(template.df, annealing.temperature, 
    regions = NULL, constraints = NULL) {
    seqs <- template.df$Sequence
    if (length(regions) != 0) {
        # compute local secondary structures
        if (length(regions) != length(seqs)) {
            stop("Length of regions does not agree with length of seqs for template secondary structure analysis.")
        }
        seqs <- sapply(seq_along(regions), function(x) substring(seqs[x], regions[[x]][1], 
            regions[[x]][2]))  # extract region substrings
    }
    struct <- compute.structure.vienna(seqs, annealing.temperature, constraints)
    return(struct)
}

#' Optimization of Template Binding
#'
#' Optimizes template binding regions according to secondary structures.
#'
#' @param template.df Template data frame.
#' @param annealing.temperature Temperature at which to compute secondary structures.
#' @param primer.lengths Target length of primers that are to be used.
#' @param mode.directionality Direction of primers.
#' @return List with new binding intervals for every template.
#' @keywords internal
optimize.template.binding.regions.single <- function(template.df, annealing.temperature, 
    primer.lengths, mode.directionality = c("fw", "rev")) {
    # determine local template regions for secondary structure analysis for fw/rev
    # direction
    
    # annealing.temperature: annealing temp for which to compute deltaG
    # primer.length: interval of minimal and maximal primer length
    # mode.directionality: adjust fw or rev binding region?
    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' argument.")
    }
    mode.directionality <- match.arg(mode.directionality)
    if (length(annealing.temperature) == 0) {
        # no annealing temp available yet choice of temperature shouldn't matter too much
        # here, since we compute Delta Delta G only
        annealing.temperature <- 25  # set to standard thermodynamic temperature. 
    }
    region.slack <- 40  # how many additional positions (outside the target region) to consider after the target region has ended. ensures 'correct folding'.
    if (mode.directionality == "fw") {
        sub.s.base <- template.df$Allowed_Start_fw
        sub.e.base <- template.df$Allowed_End_fw
        sub.e <- template.df$Allowed_End_fw + region.slack
        sub.e <- sapply(seq_along(sub.e), function(x) ifelse(sub.e[x] > nchar(template.df$Sequence)[x], 
            nchar(template.df$Sequence)[x], sub.e[x]))  # make sure that we don't overextend
        sub.s <- sub.s.base
    } else {
        sub.s.base <- template.df$Allowed_Start_rev
        sub.e.base <- template.df$Allowed_End_rev
        sub.s <- template.df$Allowed_Start_rev - region.slack
        sub.s <- sapply(seq_along(sub.s), function(x) ifelse(sub.s[x] < 1, 1, sub.s[x]))  # make sure that we don't underextend
        sub.e <- sub.e.base
    }
    if (all(sub.s == 0) || all(sub.e == 0)) {
        # nothing to optimize
        return(NULL)
    }
    regions <- lapply(seq_along(sub.s), function(x) c(sub.s[x], sub.e[x]))  # regions for secondary structure computation (extended region)
    regions.base <- lapply(seq_along(sub.s.base), function(x) c(sub.s.base[x], sub.e.base[x]))  # regions for constraining (only the target region)
    best.structs <- compute.template.secondary.structures(template.df, annealing.temperature, 
        regions, NULL)  # deltaG of unconstrained structures
    best.structs$ID <- template.df$ID
    # determine constrained structures in windows 1. create windows determine windows
    # for constrained folding
    pol.size <- 0  # is it necessary to add some space for polymerase binding without secondary structure in the way?
    win.size <- min(primer.lengths) + pol.size  # window should at least span the minimal primer size
    # very important: assure that all windows have the same size, otherwise there's a
    # biasA
    # a) individual windows -> problem for plotting
    #windows <- lapply(regions.base, function(x) Hmisc::cut2(x[1]:x[2], m = win.size, 
        #onlycuts = TRUE))
    # b) uniform windows
    windows <- replicate(length(sub.s.base), seq(min(sub.s.base), max(sub.e.base), win.size), simplify = FALSE)
    # 2. create constraint strings for selected region
    constraints <- vector("list", length(windows))
    for (i in seq_along(windows)) {
        win <- windows[[i]]
        win.pos <- win - min(win) + 1
        base.condition <- rep(".", max(regions[[i]]) - min(regions[[i]]) + 1)  # dot indicates no constraints on binding, base.condition is for the extended target region
        # message(length(base.condition))
        cons <- rep(NA, length(win) - 1)  # constraints for each window into current seq
        for (j in seq_along(win[-length(win)])) {
            # don't create windows after the target region
            condition <- base.condition
            s <- win.pos[j]
            max.pos <- sub.e.base[i] - sub.s.base[i] + 1
            if (s > length(condition) || s > max.pos) {
                next
            }
            if (win[j+1] > sub.e.base[i]) {
                # always extend the last window until the end
                #e <- win.pos[j+1]
                e <- sub.e.base[i] - sub.s.base[i] + 1
            } else {
                e <- win.pos[j+1] - 1
            }
            # adjust end to the end of the target region if it exceeds it 
            if (e > max.pos) {
                e <- max.pos
            }
            condition[s:e] <- "x"  # x indicates constraints on binding
            cons[j] <- paste(condition, collapse = "")
        }
        if (length(cons) == 0) {
            # window is a single position
            base.condition[win[1]] <- "x"
            cons <- paste(base.condition, collapse = "")
        }
        constraints[[i]] <- cons
    }
    # 3. compute constrained structures deltaG
    constrained.foldings <- NULL
    max.nbr.intervals <- max(sapply(windows, function(x) length(x) - 1))
    if (max.nbr.intervals == 0) {
        return(NULL)
    }
    for (i in seq_len(max.nbr.intervals)) {
        cur.constraints <- sapply(constraints, function(x) x[i])  # constraints for all seqs
        sel <- which(!is.na(cur.constraints))
        if (length(sel) == 0) {
            next
        }
        con.structs <- compute.template.secondary.structures(template.df[sel, ], annealing.temperature, 
            regions[sel], cur.constraints[sel])  # deltaG of constrained structures
        con.structs$Constraint <- cur.constraints[sel]
        con.structs$Constraint_Idx <- i  # the i-th window of constraints on the secondary structure
        con.structs$Constraint_Region_Start <- sapply(windows[sel], function(x) x[i])
        # determine the window ends
        if (i == max.nbr.intervals) {
            con.structs$Constraint_Region_End <- sapply(windows[sel], function(x) x[i + 1])
        } else {
            con.structs$Constraint_Region_End <- sapply(windows[sel], function(x) x[i + 1] - 1)
        }
        # compute effective region length (the constrained region)
        con.structs$Length <- unlist(lapply(strsplit(cur.constraints[sel], split = ""), function(x) length(which(x == "x"))))
        con.structs$Constraint_Region_Actual_End <- sapply(seq_along(sel), function(x) windows[sel][[x]][i] + con.structs$Length[x] - 1)
        con.structs$ID <- template.df[sel, "ID"]
        constrained.foldings <- rbind(constrained.foldings, con.structs)
    }
    colnames(constrained.foldings)[colnames(constrained.foldings) == "deltaG"] <- "DeltaG_constrained"
    colnames(constrained.foldings)[colnames(constrained.foldings) == "Structure"] <- "Structure_constrained"
    constrained.foldings$Structure_best <- best.structs$Structure[match(constrained.foldings$ID, 
        best.structs$ID)]
    
    constrained.foldings$DeltaG_best <- best.structs$deltaG[match(constrained.foldings$ID, 
        best.structs$ID)]
    # Delta Delta G: the difference between the constrained and unconstrained folding: the smaller Delta Delta G, the better (in this case there is little benefit of folding with secondary structure in the current region)
    constrained.foldings$Delta_Delta_G <- constrained.foldings$DeltaG_constrained - 
        constrained.foldings$DeltaG_best
    # annotate with geneGroups
    constrained.foldings$Group <- template.df$Group[match(constrained.foldings$ID, 
        template.df$ID)]
    # plot.Delta.DeltaG(constrained.foldings, stratify = TRUE)
    # my_ggsave(file.path(results.loc,
    # 'template_secondary_structure_delta_delta_G.pdf')) 4. select windows: a)
    # determine standard deviation of Delta Delta G for each template
    deviations <- sapply(seq_len(nrow(template.df)), function(x) 
                    sd(constrained.foldings$Delta_Delta_G[constrained.foldings$ID == template.df[x, "ID"]]))
    myfun <- function(len, deltaG, win.size) {
        return(deltaG[which(len >= win.size)])
    }
    # don't select windows that are too short for defining the minDeltaG
    minG.per.template <- ddply(constrained.foldings, c("ID", "Group"), plyr::here(plyr::summarize),
        MinG = myfun(substitute(Length), substitute(Delta_Delta_G), win.size))
    m <- match(template.df$ID, minG.per.template$ID)
    minG.per.template <- minG.per.template[m,]
    # determine the maximal DeltaDeltaG of the optimized regions:
    maxG.per.template <- minG.per.template$MinG + deviations
    # b) retrieve windows in each template with DeltaDeltaG <= maxG.per.template
    sel.windows <- lapply(seq_along(template.df$ID), function(x) {
        temp.data <- constrained.foldings[constrained.foldings$ID == template.df$ID[x], ]
        allowed.deltaG <- maxG.per.template[match(template.df$ID[x], minG.per.template$ID)]
        sel <- temp.data$Delta_Delta_G <= allowed.deltaG
        con <- temp.data$Constraint_Idx[sel]
        # c) select longest consecutive windows for implementation reasons
        if (length(con) > 1) {
            # multiple best windows found
            result <- rle(diff(con))
            idx <- which(result$values == 1)
            if (length(idx) == 0) {
                # no consecutive runs -> select window with smallest DeltaDeltaG that has allowed effective length
                data <- temp.data[sel, ]
                len.sel <- which(data$Length >= win.size)
                con <- con[len.sel[which.min(data$Delta_Delta_G[len.sel])]]
            } else {
                longest.runs <- idx[which.max(result$lengths[idx])]
                if (longest.runs == 1) {
                  sel.start <- 1
                } else {
                  sel.start <- sum(result$length[1:(longest.runs - 1)]) + 1
                }
                sel.end <- sum(result$length[1:longest.runs]) + 1
                con <- con[sel.start:sel.end]
            }
        } else {
            # only one best window found
            con
        }
    })
    new.regions <- lapply(seq_along(template.df$ID), function(x) {
        sel <- which(constrained.foldings$ID == template.df$ID[x])  # retrieve start/end of regions
        sel <- sel[sel.windows[[x]]]
        s <- min(constrained.foldings$Constraint_Region_Start[sel])  # retrieve selected consecutive regions
        e <- max((constrained.foldings$Constraint_Region_Actual_End[sel])) # use the actual end of the region here
        return(c(s, e))
    }) 
    constrained.foldings$Direction <- mode.directionality
    # don't output incomplete regions
    constrained.foldings <- constrained.foldings[constrained.foldings$Length >= win.size,]
    result <- list("Intervals" = new.regions, "Foldings" = constrained.foldings)
    return(result)
}
#' Optimization of Binding Regions
#' 
#' Optimizes the template binding regions.
#'
#' @param template.df Template data frame.
#' @param annealing.temperature Temperature at which to compute secondary structures.
#' @param primer.lengths Target length of primers that are to be used.
#' @param mode.directionality Direction of primers.
#' @return List with intervals indicating improved primer binding regions.
#' @keywords internal
optimize.template.binding.regions.dir <- function(template.df, annealing.temperature = NULL, 
    primer.lengths, mode.directionality = c("fw", "rev", "both")) {
    # compute uniform result for all possible direction inputs
    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' argument.")
    }
    mode.directionality <- match.arg(mode.directionality)
    if (mode.directionality == "both") {
        binding.regions.fw <- optimize.template.binding.regions.single(template.df, annealing.temperature, 
            primer.lengths, "fw")
        binding.regions.rev <- optimize.template.binding.regions.single(template.df, annealing.temperature, 
            primer.lengths, "rev")
        intervals <- list(fw = binding.regions.fw$Intervals, rev = binding.regions.rev$Intervals)
        foldings <- plyr::rbind.fill(binding.regions.fw$Foldings, binding.regions.rev$Foldings)
    } else {
        binding.regions <- optimize.template.binding.regions.single(template.df, annealing.temperature, 
            primer.lengths, mode.directionality)
        intervals <- list(binding.regions$Intervals)
        names(intervals) <- mode.directionality
        foldings <- binding.regions$Foldings
    }
    result <- list("Intervals" = intervals, "Foldings" = foldings)
    return(result)
}
#' Delta DeltaG Plot
#'
#' Plots the difference between the free energy of constrained 
#' and unconstrained foldings.
#'
#' @param constrained.foldings Data frame with info from constrained foldings.
#' @param stratify Stratify according to template groups?
#' @return Plot of Delta DeltaG.
#' @keywords internal
plot.Delta.DeltaG <- function(constrained.foldings, stratify = FALSE) {
    # plot line plot of DeltaDeltaG (constrained - MFE structure) for all templates
    # at once, stratify by class?
    plot.df <- constrained.foldings
    groups <- max(constrained.foldings$Constraint_Idx)
    plot.df$Binding_Region <- Hmisc::cut2(c(1, constrained.foldings$Constraint_Region_End), 
        g = groups)[-1]
    
    if (!stratify) {
        p <- ggplot(plot.df, aes_string(x = "Binding_Region", y = "Delta_Delta_G")) + geom_bar(stat = "summary")
    } else {
        p <- ggplot(plot.df, aes_string(x = "Binding_Region", y = "Delta_Delta_G")) + geom_boxplot() + 
            facet_wrap(~Group)
    }
    x.names <- unique(paste(plot.df$Constraint_Region_Start, "-", plot.df$Constraint_Region_End, 
        sep = ""))
    title <- "Tendency for secondary structure formation in template target regions"
    # with geom_bar
    p <- p + xlab("Region constrained to be free of base pairing") + ylab("Delta Delta G") + 
        ggtitle(title)
    return(p)
}


