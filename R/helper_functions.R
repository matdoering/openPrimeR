#' Identification of File extension.
#'
#' Identifies the file extension of \code{x}.
#' @param x A string for a filename.
#'
#' @return The extension of \code{x}.
#' @keywords internal
get.extension <- function(x) {
    pos <- regexpr("\\.([[:alnum:]]+)$", x)
    return(ifelse(pos > -1L, substring(x, pos + 1L), ""))
}
################ 
#' Copy Directories.
#'
#' Copies a directory to another location.
#' 
#' @param src.dir The directory to be copied.
#' @param dest.dir The target directory.
#' @param overwrite Overwrite existing files in \code{dest.dir}.
#' @return TRUE if copying was successful, FALSE otherwise.
#' @keywords internal
dir.copy <- function(src.dir, dest.dir, overwrite) {
    file.names <- dir(src.dir)
    suppressWarnings(dir.create(dest.dir))
    return(file.copy(from = file.path(src.dir, file.names), to = dest.dir, overwrite = overwrite, 
        recursive = TRUE))
}
#' Index for Unlisting.
#'
#' Determines indices for unlisting.
#'
#' @param primer.start Numeric vector.
#' @param primer.data.idx Selection indices.
#'
#' @return Indices.
#' @keywords internal
get.unlist.idx <- function(primer.start, primer.data.idx) {
    un <- unique(primer.data.idx)
    t <- rep(1, length(un))
    cat.count <- table(primer.data.idx)
    cat <- as.numeric(names(cat.count))
    names(t) <- un
    result <- rep(NA, length(primer.data.idx))
    prev.count <- 0
    for (i in seq_along(primer.data.idx)) {
        idx <- t[as.character(primer.data.idx[i])]
        cat.c <- sum(cat.count[which(cat < primer.data.idx[i])])  # nbr of list entries appearing before
        result[i] <- idx + cat.c
        t[as.character(primer.data.idx[i])] <- t[as.character(primer.data.idx[i])] + 
            1
    }
    return(result)
}

#' Rbind for Primer Data Frames.
#'
#' Merges all primer data frames in \code{primer.data} into one data frame.
#'
#' @param primer.data List with primer data frames.
#' @return A data frame containing all data in \code{primer.data}.
#' @keywords internal
rbind.primer.data <- function(primer.data) {
    if (length(primer.data) == 0 || all(sapply(primer.data, function(x) nrow(x) == 
        0))) {
        return(NULL)
    }
    # select non-zero primer sets
    idx <- sapply(primer.data, function(x) nrow(x) != 0)
    run.names <- get.run.names(primer.data)[idx]
    unique.names <- make.unique(run.names)
    idx <- which(idx)
    names(primer.data)[idx] <- unique.names  # list names need to be unique for my_rbind!
    data <- my_rbind(primer.data[idx])
    return(data)
}
#' Format Strings
#'
#' Changes the representation of the comma-separated string input.
#'
#' @param values A comma-separated string with values.
#' @return A percentage-formatted representation of the input string.
#' @keywords internal
string.list.format.total <- function(values) {
    # format a list of values to percentages for the whole data set
    na.idx <- which(is.na(values))
    if (length(na.idx) != 0) {
        values[na.idx] <- ""
    }
    v <- unlist(strsplit(values, split = ","))
    t <- table(v)
    t <- t[order(t, decreasing = TRUE)]
    res <- paste(names(t), "(", round((t/sum(t)) * 100, 2), "%)", collapse = ",", 
        sep = "")
    return(res)
}
#############
# helper functions for plotting
#############
#' Plot Extent
#'
#' Returns the extent of a plot.
#'
#' @param N Number of observations to plot.
#' @param px.per.n Pixels required per observations.
#' @param min.size Minimal extent of plot in pixels.
#' @param max.size Maximal extent of plot in pixels.
#' @return The extent of the plot.
#' @keywords internal
get.plot.height <- function(N, px.per.n = 50, min.size = 300, 
    max.size = 1500000) {

    height <- min(N * px.per.n + min.size, max.size)
    return(height)
}
#' Smartbind preserving classes.
#'
#' Rbind allowing for column mismatch, retains the classes of the data frames.
#' Motivation: smartbind/rbind.fill only keep the data.frame class
#' but not additional classes.
#'
#' @param ... Data frames.
#' @return A data frame resulting from row binding of \code{...}.
#' @keywords internal
my_rbind <- function(...) {
    df <- plyr::rbind.fill(...)
    args <- list(...)
    classes <- sapply(args, class)
    if (length(args) == 0) {
        return(NULL)
    }
    if (any(classes == "Templates")) {
       df <- Templates(df) 
    } 
    if (any(classes == "Primers")) {
        df <- Primers(df)
    }
    return(df)
}
#' Conversion of Positions to Ranges.
#' 
#' Converts two numeric values to a range.
#'
#' @param pos1 The first value.
#' @param pos2 The second value.
#' @return A character vector range.
#' @keywords internal
pos.to.range <- function(pos1, pos2) {
    s1 <- strsplit(pos1, split = ",")
    s2 <- strsplit(pos2, split = ",")
    res <- unlist(lapply(seq_along(s1), function(x) if (length(s1[[x]]) != 0 && length(s2[[x]]) != 
        0) 
        paste(s1[[x]], " to ", s2[[x]], collapse = ",", sep = "") else ""))
    return(res)
}
