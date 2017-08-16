############
# Initialization of Primers
############

#' Primer Length Check.
#'
#' Checks whether it is possible to construct primers of the desired length.
#'
#' @param template.df Template data frame.
#' @param allowed.region.definition Definition of allowed binding regions.
#' @param primer.lengths The desired lengths of the priemrs.
#' @param mode.directionality The primer directionality.
#' @return TRUE, if primers of the desired length can be constructed, 
#' @keywords internal
#' FALSE otherwise.
check.init.primer.length <- function(template.df, 
    allowed.region.definition = c("within", "any"), primer.lengths, 
    mode.directionality = c("fw", "rev", "both")) {

    if (length(allowed.region.definition) == 0) {
        stop("Please supply the 'allowed.region.definition' arg.")
    }
    allowed.region.definition <- match.arg(allowed.region.definition)
    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' arg.")
    }
    mode.directionality <- match.arg(mode.directionality)
    if (length(template.df) == 0) {
        stop("No templates available.")
    }
    if (length(primer.lengths) == 0) {
        stop("No primer length specified.")
    }
    check <- NULL
    if (mode.directionality == "fw") {
        check <- check.init.primer.length.single(template.df$Allowed_fw, allowed.region.definition, 
            min(primer.lengths))
    } else if (mode.directionality == "rev") {
        check <- check.init.primer.length.single(template.df$Allowed_rev, allowed.region.definition, 
            min(primer.lengths))
    } else {
        check.fw <- check.init.primer.length.single(template.df$Allowed_fw, allowed.region.definition, 
            min(primer.lengths))
        check.rev <- check.init.primer.length.single(template.df$Allowed_rev, allowed.region.definition, 
            min(primer.lengths))
        check <- check.fw && check.rev
    }
    return(check)
}
#' Primer Length Check.
#'
#' Checks whether it is possible to construct primers of the desired length.
#'
#' @param allowed String containing the allowed binding sequence.
#' @param allowed.region.definition Definition of allowed binding regions.
#' @param min.len Minimal desired primer lengths.
#' @return TRUE if primers of the desired length can be constructed,
#' FALSE otherwise.
#' @keywords internal
check.init.primer.length.single <- function(allowed, 
    allowed.region.definition = c("within", "any"), min.len) {
    # allowed: allowed binding region min.len: minimal primer length required by user
    if (length(allowed.region.definition) == 0) {
        stop("Please supply the 'allowed.region.definition' arg.")
    }
    allowed.region.definition <- match.arg(allowed.region.definition)
    if (allowed.region.definition == "within" && any(nchar(allowed) < min.len)) {
        warning("Selected region is too short to create primers for all templates.")
        return(FALSE)
    } else {
        return(TRUE)
    }
}
#' Creation of an Initial Primer Set.
#'
#' Creates an initial set of candidate primers for primer design.
#'
#' @param template.df Template data frame.
#' @param primer.lengths Vector containing the permissible primer lengths.
#' @param mode.directionality Direction of primers to be created.
#' @param sample Name of the template sample.
#' @param allowed.region.definition Definition of the allowed binding region.
#' @param init.algo Algorithm for initializing primers.
#' @param max.degen Maximal allowed degeneration of created primers.
#' @param conservation Required conservation of primers.
#' The value of \code{conservation} should be in the range[0,1].
#' @param updateProgress Shiny progress object.
#' @return An initialized data frame of candidate primers.
#' @keywords internal
create.initial.primer.set <- function(template.df, primer.lengths, mode.directionality = c("fw", "rev"), 
    sample, allowed.region.definition = c("within", "any"), 
    init.algo = c("naive", "tree"), max.degen, conservation, updateProgress = NULL) {

    if (length(init.algo) == 0) {
        stop("Please supply the 'init.algo' arg.")
    }
	init.algo <- match.arg(init.algo)
    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' arg.")
    }
	mode.directionality <- match.arg(mode.directionality)
    if (length(allowed.region.definition) == 0) {
        stop("Please supply the 'allowed.region.definition' arg.")
    }
    allowed.region.definition <- match.arg(allowed.region.definition)
    seqs <- template.df$Sequence
    if (mode.directionality == "fw") {
        # Set the positions in the sense sequence
        l.s <- template.df$Allowed_Start_fw
        e.s <- template.df$Allowed_End_fw
    } else if (mode.directionality == "rev") {
        # Start primer design start/end in reverse complement sequence
        l.s <- nchar(template.df$Sequence) - template.df$Allowed_End_rev + 1 # smaller number (e.g. 1)
        e.s <- nchar(template.df$Sequence) - template.df$Allowed_Start_rev + 1 # larger number (e.g. 30)
        seqs <- rev.comp.sequence(seqs)
    }
    # if l.s/e.s is NA -> create primers for the full sequence
    l.s[is.na(l.s)] <- 1
    e.s[is.na(e.s)] <- nchar(seqs)[is.na(e.s)]
    # ensure that required primer lengths can be fulfilled
    if (allowed.region.definition == "within") {
        if (any(max(primer.lengths) > (e.s - l.s + 1))) {
            stop("The specified primer length cannot be obtained with the specified binding regions.")
        }
    } else {
        # any intersection with binding range: primer length may not exceed the sequence length
        if (any(max(primer.lengths) > nchar(template.df$Sequence))) {
            stop("Specified primer length exceeded the length of at least one sequence.")
        }
    }
    # change group IDs to be unique
    group.tags <- uniqtag::uniqtag(unique(template.df$Group), 3)
    seq.groups <- group.tags[match(template.df$Group, unique(template.df$Group))]
    seq.IDs <- uniqtag::uniqtag(template.df$ID, 5)
    if (init.algo == "naive") {
        primers <- create.primers.naive(seqs, seq.IDs, seq.groups, l.s, e.s, primer.lengths, 
            allowed.region.definition, max.degen, sample, mode.directionality, updateProgress)
    } else {
        primers <- create.primers.tree(seqs, seq.IDs, seq.groups, l.s, e.s, primer.lengths, 
            allowed.region.definition, max.degen, conservation, sample, mode.directionality, 
            updateProgress)
    }
    return(primers)
}
#' Primer Identifier Creation.
#' 
#' Creates identifiers for generated primers.
#'
#' @param sample Sample name of the templates.
#' @param seq.IDs Identifiers of the templates.
#' @param seq.identifiers The index of the seq.
#' @param all.starts Primer positions (start).
#' @param all.ends Primer positions (end).
#' @param identifier Direction keyword.
#' @param seq.primers The primer sequences as strings.
#' @return Identifiers for each primer.
#' @keywords internal
get.primer.identifier.string <- function(sample, seq.IDs, seq.identifier, all.starts, all.ends,
                                    identifier, seq.primers) {
    enum <- seq_along(seq.IDs)
    seq.names.short <- paste(paste0(seq.identifier, "-", enum), "|", seq.IDs, "|", 
                        all.starts, ":", all.ends, "|_",
                        identifier, sep = "")
    return(seq.names.short)
}
#' Naive Initialization of Primers.
#'
#' Initialize primers by extracting substrings from all templates.
#'
#' @param seqs The template sequence strings.
#' @param seq.IDs The identifiers of the templates.
#' @param seq.groups The group identifiers of the templates.
#' @param l.s The positions where the allowed region starts for each template.
#' @param e.s The positions where the allowed reigon ends for each template.
#' @param primer.lengths Vector of desired primer lengths.
#' @param allowed.region.definition Definition of the allowed region.
#' @param max.degen Maximum allowed degeneracy of primers.
#' @param sample Template sample identifier.
#' @param updatProgress Shiny progress object.
#' @return Data frame with initialized primer candidates.
#' @keywords internal
create.primers.naive <- function(seqs, seq.IDs, seq.groups, l.s, e.s, primer.lengths, allowed.region.definition, 
    max.degen, sample = "", identifier = "", updateProgress = NULL) {
    if (length(seqs) == 0) {
        return(NULL)
    }
    seq.idx <- NULL
    # use group identifier for primers?
    use.group.identifier <- TRUE
    if (length(unique(seq.groups)) <= 1) {
        use.group.identifier <- FALSE
    }
    all.primer.seqs <- foreach(seq.idx = seq_along(seqs), .combine = c) %dopar% {
        if (is.function(updateProgress)) {
            detail <- identifier
            updateProgress(1/length(seqs), detail, "inc")
        }
        # create possible primer start positions:
        if (allowed.region.definition == "within") {
            primer.starts <- l.s[seq.idx]:(max(1, e.s[seq.idx] - min(primer.lengths)))
        } else {
            # any overlap:
            primer.starts <- (max(1, l.s[seq.idx] - max(primer.lengths))):(e.s[seq.idx])
        }
        # create allowed primer intervals
        all.starts <- unlist(lapply(primer.starts, function(x) rep(x, length(primer.lengths))))
        all.ends <- unlist(lapply(seq_along(primer.starts), function(x) primer.starts[x] + primer.lengths - 1))
        if (allowed.region.definition == "within") {
            # remove primers that extend over the target region
            rm.idx <- which(all.ends > e.s[seq.idx])
        } else {
            # remove primers that would exceed the template length
            rm.idx <- which(all.ends > nchar(seqs[seq.idx]))
        }
        if (length(rm.idx) != 0) {
            all.starts <- all.starts[-rm.idx]
            all.ends <- all.ends[-rm.idx]
        }
        seq.primers <- substring(seqs[seq.idx], all.starts, all.ends)
        # identify primers: 
        if (use.group.identifier) {
            # use group identifiers 
            use.ids <- seq.groups[seq.idx]
        } else {
            # use the first sequence identifier
            use.ids <- seq.IDs[seq.idx]
        }
        names(seq.primers) <- get.primer.identifier.string(sample, abbreviate(use.ids, 5),  
                                    seq.idx, all.starts, all.ends, identifier, seq.primers)
        seq.primers
    }
    all.primer.seqs <- all.primer.seqs[!duplicated(all.primer.seqs)]
    degen <- score_degen(strsplit(all.primer.seqs, split = ""))
    all.primer.seqs <- all.primer.seqs[degen <= max.degen]
    return(all.primer.seqs)
}
#' Creation of Candidate Primers.
#' 
#' Creates a set of primer candidates based on the input template
#' sequences. This set of primers can be used to create
#' custom primer design algorithms. 
#'
#' @param sample Character vector providing an identifier for the templates.
#' @param template.df An object of class \code{Templates} providing
#' the template sequences for which an initial set of primers is constructed.
#' @param primer.lengths Numeric interval providing the 
#' minimal and maximal allowed lengths of generated primers.
#' @param mode.directionality Character vector giving the direction of primers to be created. Either "fw" to create forward primers 
#' or "rev" to create reverse primers. The default is "fw".
#' @param allowed.region.definition A character vector providing the definition
#' of region where primers are to be constructed.
#' If \code{allowed.region.definition} is "within", constructed primers lie within the allowed binding region.
#' If \code{allowed.region.definition} is "any", primers overlap with the allowed binding region.
#' The default is "within".
#' @param init.algo Algorithm for initializing primers.
#' If you set \code{init.algo} to "tree", then initial primers will be generated
#' by forming degenerate consensus sequences using a tree-based approach. 
#' This option requires MAFFT for multiple alignments (see notes). 
#' If \code{init.algo} is set to "naive", initial primers are
#' generated by extracting substrings from the template target regions.
#' @param max.degen A numeric providing the maximal allowed degeneration of created primers.
#' @param conservation The percentile of top-conserved template
#' regions to consider. The value of \code{conservation} should be 
#' in the range[0,1].
#' If the \code{conservation} is set to 1 (the default), all regions 
#' are considered.
#' @param updateProgress A Shiny progress object; by default
#' this is set to \code{NULL} such that no callback is used.
#' @return A data frame with candidate primers for optimization.
#' @note If you want to set \code{init.algo} to "tree", please install
#' MAFFT (http://mafft.cbrc.jp/alignment/software/) on your computer.
#' @export 
#' @family primer functions
#' @examples
#' data(Ippolito)
#' # Naive primer initialization
#' init.primers <- get_initial_primers("InitialPrimers", template.df, 
#'                          c(18,18), "fw", init.algo = "naive")
#' # Tree-based primer initialization (requires MAFFT)
#' \dontrun{
#' init.primers <- get_initial_primers("InitialPrimers", template.df, 
#'                          c(18,18), "fw", init.algo = "tree")
#' }
get_initial_primers <- function(sample, template.df, primer.lengths, 
    mode.directionality = c("fw", "rev"), allowed.region.definition = c("within", "any"), 
    init.algo = c("naive", "tree"), max.degen = 16, conservation = 1.0, updateProgress = NULL) {

    if (length(init.algo) == 0) {
        stop("Please supply the 'init.algo' arg.")
    }
	init.algo <- match.arg(init.algo)
    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' arg.")
    }
	mode.directionality <- match.arg(mode.directionality)
    if (length(allowed.region.definition) == 0) {
        stop("Please supply the 'allowed.region.definition' arg.")
    }
    allowed.region.definition <- match.arg(allowed.region.definition)
    if (length(primer.lengths) != 2) {
        stop("Primer lengths should be an interval")
    }
    if (primer.lengths[2] < primer.lengths[1]) {
        stop("Primer length doesn't specify an interval")
    }
    if (conservation < 0 || conservation > 1) {
        stop("The conservation percentile should be in the range [0,1]")
    }
    primer.lengths <- primer.lengths[1]:primer.lengths[2]
    if (any(primer.lengths <= 0)) {
        stop("All primer lengths have to be strictly positive.")
    }
    if (max.degen <= 0) {
        stop("The maximal primer degeneracy has to be strictly positive.")
    }
    if (is.function(updateProgress)) {
        detail <- "Creating primers"
        updateProgress(1/2, detail, "inc")
    }
    primer.set <- create.initial.primer.set(template.df, primer.lengths, mode.directionality, 
        sample, allowed.region.definition, init.algo, max.degen, conservation, updateProgress)
    if (is.function(updateProgress)) {
        detail <- "Structuring"
        updateProgress(1/2, detail, "inc")
    }
    primer.df <- read_primers.internal(primer.set, names(primer.set), "fw", "rev", 
        merge.ambig = "none", max.degen, sample)
    return(primer.df)
}
#' Creation of Initial Primers
#' 
#' Creates a set of candidate primers.
#'
#' @param template.df Template data frame.
#' @param sample.name Name of the template sample.
#' @param primer.lengths Interval of minimal and maximal desired primer length.
#' @param mode.directionality Direction of primers to be created.
#' @param init.algo Algorithm for initializing primers.
#' @param allowed.region.definition Definition of the allowed binding region.
#' @param max.degen Maximal allowed degeneration of created primers.
#' @param conservation Required conservation of primers.
#' The value of \code{conservation} should be in the range[0,1].
#' @param cur.results.loc Location for writing the primers as csv.
#' @return An initial primer data frame.
#' @keywords internal
initialize.primer.set <- function(template.df, sample.name, primer.lengths, allowed.region.definition, 
    mode.directionality, init.algo, max.degen, conservation, cur.results.loc) {
    out.name <- get.init.file.name(cur.results.loc, sample.name, primer.lengths, 
        mode.directionality, allowed.region.definition, init.algo, max.degen, conservation)


    primer.df <- get_initial_primers(sample.name, template.df, primer.lengths, mode.directionality, 
        allowed.region.definition, init.algo, max.degen, conservation)
    if (length(out.name) != 0) {
        dir.create(dirname(out.name), showWarnings = FALSE, recursive = TRUE)
        write.csv(primer.df, file = out.name, row.names = FALSE)
    }
    return(primer.df)
}
