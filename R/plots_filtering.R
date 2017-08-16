#' Plot Filtering Runtimes
#'
#' Plots the runtimes of individual evaluation steps in the filtering procedure.
#'
#' @param filered.stats Stats from filtering.
#' @return A plot of the runtime for each filtering step.
#' @keywords internal
plot.filtering.runtime <- function(filtered.stats) {
    plot.stats <- filtered.stats[which(filtered.stats$Constraint != ""), ]  # make sure that we don't have duplicates when we form levels if filtering stopped early
    plot.stats$Constraint <- factor(plot.stats$Constraint, levels = rev(unique(plot.stats$Constraint)))
    constraint.names <- unlist(constraints_to_unit(levels(plot.stats$Constraint), FALSE))
    if (nrow(filtered.stats) == 0) {
        return(NULL)
    }
    ggplot(plot.stats, aes_string(y = "Time", x = "Constraint", fill = "Time")) + geom_bar(stat = "identity") + 
        coord_flip() + xlab("Applied constraint") + ylab("Required computation time [s]") + 
        ggtitle("Runtime for each filtering step") + 
        theme(axis.text.x = element_text(hjust = 1)) +
        scale_fill_gradient() + 
        facet_wrap(~Direction) +
        scale_x_discrete(labels = constraint.names)
}
#' Plot of Filtering Stats for Coverage.
#'
#' Plots the remaining coverage after each filtering step.
#'
#' @param stats Statistics of the filtering procedure.
#' @param stats.relax Statistic of the relaxation procedure.
#' @return A plot showing the possible coverage after each filtering step.
#' @keywords internal
plot.filtering.stats.cvg <- function(stats, stats.relax = NULL) {
    plot.stats <- stats_plot_data(stats, stats.relax)
    if (length(plot.stats) == 0 || nrow(plot.stats) == 0) {
        return(NULL)
    }
    if (all(is.na(plot.stats$Current_Coverage))) { # no coverage available
        return(NULL)
    }
    yl <- "Primer coverage after filtering"
    ttl <- "Filtering coverage"
    filter.color <- RColorBrewer::brewer.pal(8, "Accent")[5]
    relax.color <- RColorBrewer::brewer.pal(8, "Blues")[4]
    color <- c("Filtering" = filter.color, "Relaxation" = relax.color)

    constraint.names <- unlist(constraints_to_unit(levels(plot.stats$Constraint), FALSE))
    ggplot(plot.stats, aes_string(y = "Current_Coverage", x = "Constraint", fill = "Type")) + 
        geom_bar(stat = "identity") + coord_flip() + xlab("Applied constraint") + 
        ylab(yl) + ggtitle(ttl) + 
        scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
        theme(axis.text.x = element_text(hjust = 1)) +
        scale_fill_manual(values = color) + 
        facet_wrap(~Direction) + 
        scale_x_discrete(labels = constraint.names)
}
#' Combination of Filtering Stats.
#'
#' Summarizes filtering/relaxation statistics for plotting.
#'
#' @param stats Statistics of the filtering procedure.
#' @param stats.relax Statistic of the relaxation procedure.
#' @return A data frame combinin filtering/relaxation stats.
#' @keywords internal
stats_plot_data <- function(stats, stats.relax) {
    if (length(stats) == 0 && nrow(stats) == 0) {
        return(NULL)
    }
    plot.stats <- stats
    plot.stats$Constraint <- factor(plot.stats$Constraint, levels = rev(unique(plot.stats$Constraint)))
    ratios <- plyr::ddply(plot.stats, "Direction", plyr::summarize, 
                Ratio = substitute(Remaining) / (substitute(Remaining)[1] + substitute(Excluded)[1]))$Ratio
    plot.stats$Ratio <- ratios
    plot.stats$Type <- "Filtering"
    if (length(stats.relax) != 0 && nrow(stats.relax) != 0) {
        # augment with relaxation info
        plot.relax <- stats.relax
        plot.relax$Type <- "Relaxation"
        # don't take first element for reverse .. split here
        fw.idx <- plot.stats$Direction == "fw"
        fw.count <- plot.stats$Remaining[fw.idx][1] + plot.stats$Excluded[fw.idx][1]
        rev.idx <- plot.stats$Direction == "rev"
        rev.count <- plot.stats$Remaining[rev.idx][1] + plot.stats$Excluded[rev.idx][1]
        fw.ratio <- plot.relax$Remaining[plot.relax$Direction == "fw"] / fw.count
        rev.ratio <- plot.relax$Remaining[plot.relax$Direction == "rev"] / fw.count
        plot.relax$Ratio <- c(fw.ratio, rev.ratio)
        # adjust actual coverage using the relaxation info
        plot.relax$Current_Coverage[is.na(plot.relax$Current_Coverage)] <- 0
        cum.sums <- plyr::ddply(plot.relax, "Direction", plyr::summarize, 
                Cumulative_Coverage = cumsum(substitute(Current_Coverage)))
        plot.relax$Current_Coverage <- cum.sums$Cumulative_Coverage
        idx <- which(!is.na(plot.relax$Current_Coverage))
        if (length(idx) != 0) {
            plot.stats$Current_Coverage[idx] <- sapply(seq_along(idx), function(x) max(0,plot.stats$Current_Coverage[idx[x]] - plot.relax$Current_Coverage[idx[x]]))
        }
        plot.stats <- rbind(plot.stats, plot.relax)
        plot.stats$Type <- factor(plot.stats$Type, levels = c("Relaxation", "Filtering"))
    } 
    return(plot.stats)
}
#' Plot of Overall Filtering Stats.
#' 
#' Plots the number of primers remaining after each filtering step.
#'
#' @param stats Statistics on the filtering procedure
#' @param stats.relax Statistic on the filtering procedure after relaxation.
#' @return A plot for the number of primers after filtering. 
#' @keywords internal
plot.filtering.stats <- function(stats, stats.relax = NULL) {
    plot.stats <- stats_plot_data(stats, stats.relax)
    if (length(plot.stats) == 0 || nrow(plot.stats) == 0) {
        return(NULL)
    }
    constraint.names <- unlist(constraints_to_unit(levels(plot.stats$Constraint), FALSE))
    filter.color <- RColorBrewer::brewer.pal(8, "Accent")[5]
    relax.color <- RColorBrewer::brewer.pal(8, "Blues")[4]
    color <- c("Filtering" = filter.color, "Relaxation" = relax.color)
    ggplot(plot.stats, aes_string(y = "Ratio", x = "Constraint", fill = "Type")) + 
        geom_bar(stat = "identity") + coord_flip() + xlab("Filtering constraint") + 
        ylab("Remaining primers") + 
        ggtitle("Filtering overview") + 
        scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
        scale_fill_manual(values = color) + 
        theme(axis.text.x = element_text(hjust = 1)) +
        facet_wrap(~Direction) + 
        scale_x_discrete(labels = constraint.names)
}
