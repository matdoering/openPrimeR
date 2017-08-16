##########################
# Primer functionalities
##########################
#' Creation of k-mers of a single sequence.
#' @param seq A character vector.
#' @param k The size of the k-mer.
#' @return A names character vector, where the names 
#' are the relative positions of the k-mers and the values
#' give the character vector of the k-mer.
#' @keywords internal
create.kmer <- function(seq, k) {
    # generates all unique substrings of length k contained in the single sequence 'seq' 
    if (length(seq) != 1) {
        stop("This function is only for single sequences")
    }
    if (k > nchar(seq)) {
        k <- nchar(seq)
        warning("Specified k-mer size was larger than sequence length; ",
                "Reduced k-mer size to: ", k, ".")
    }
    k.mers <- unique(substring(seq, 1:(nchar(seq) - k + 1), k:nchar(seq)))
    kmer.pos <- -seq(nchar(seq), k)
    names(k.mers) <- kmer.pos
    return(k.mers)
}
#' Creation of k-mers for multiple sequences.
#' @param seq A character vector.
#' @param k The size of the k-mer.
#' @return A list with named character vectors, containing the k-mers.
#' @keywords internal
create.k.mers <- function(seqs, k) {
    # generates all substrings of length k contained in the single sequence 'seq' 
    i <- NULL
    k.mers <- foreach (i = 1:length(seqs), .combine = c) %dopar% {
        k.mer <- list(create.kmer(seqs[i], k))
    }
    return(k.mers)
}
#' Estimation of Primer Coverage.
#'
#' Estimates the possible coverage of primers using
#' probes of size \code{k} and only considering perfect
#' matches without consideration of ambiguities.
#'
#' @param seqs A character vector of sequences to evaluate coverage for.
#' @param k A numeric giving the size of the primers.
#' @param id An optional identifier for the primers.
#' @return A data frame with binding information.
#' @keywords internal
estimate.cvg.dir <- function(seqs, k, id = "") {
    # 1. create all k-mers for input seqs
    kmers <- create.k.mers(seqs, k)
    # since we're creating before the target region, k mer pos is negative
    kmer.df <- do.call(rbind, lapply(seq_along(kmers), function(x) data.frame("Template_ID" = x, "Primer" = kmers[[x]], "Start" = as.numeric(names(kmers[[x]])), "End" = as.numeric(names(kmers[[x]])) + k)))
    # use ddply to determine binding positions / cvg
    binding.df <- plyr::ddply(kmer.df, "Primer", plyr::summarize, 
                    "Start" = paste(substitute(Start), collapse = ","), 
                    "End" = paste(substitute(End), collapse = ","),
                    "Covered_Seqs" = paste(substitute(Template_ID), collapse = ","),
                    "primer_coverage" = length(substitute(Template_ID)))
    binding.df$Coverage <- unlist(lapply(strsplit(binding.df$Start, split = ","), function(x) length(x)))
    binding.df$Coverage_Ratio <- binding.df$Coverage / length(seqs)
    return(binding.df)
}
#' Estimation of Primer Coverage.
#'
#' Estimates the possible coverage of primers using
#' probes of size \code{k} and only considering perfect
#' matches without consideration of ambiguities.
#'
#' @param seqs A character vector of sequences to evaluate coverage for.
#' @param k A numeric giving the size of the primers.
#' @param mode.directionality Estimation of coverage for forward/reverse/both?
#' @param sample An optional identifier for the sample.
#' @return A list with entries \code{fw} and \code{rev} giving
#' data frames for forward/reverse binding.
#' @keywords internal
estimate.cvg <- function(lex.df, k = 18, mode.directionality, sample = "") {
    cvg.fw <- NULL
    cvg.rev <- NULL
    if (mode.directionality == "fw") {
        # seqs <- lex.df$Allowed_fw
        cvg.fw <- estimate.cvg.dir(lex.df$Allowed_fw, k, id = paste0(sample, "_", "fw"))
    } else if (mode.directionality == "rev") {
        cvg.rev <- estimate.cvg.dir(lex.df$Allowed_rev, k, id = paste0(sample, "_", "rev"))
    } else {
        cvg.fw <- estimate.cvg.dir(lex.df$Allowed_fw, k, id = paste0(sample, "_", "fw"))
        cvg.rev <- estimate.cvg.dir(lex.df$Allowed_rev, k, id = paste0(sample, "_", "rev"))
    }
    return(list("fw" = cvg.fw, "rev" = cvg.rev))
}
#' Classification of the Difficulty of a Primer Design Task.
#'
#' Uses reference beta distributions of primer coverage ratios to
#' classify a primer design task into the groups ranging from \emph{easy} 
#' to \emph{hard}. For \emph{easy} tasks, it should
#' not be a problem to design a small primer set. For
#' \emph{hard} tasks, however, a small set of primers may not be
#' achievable.
#'
#' The difficulty of a primer design task is evaluated by
#' estimating the distribution of coverage ratios per primer
#' by performing exact string matching with 
#' primers of length \code{primer.length}, which are constructed
#' by extracting template subsequences. Next, a beta distribution
#' is fitted to the estimated coverage distribution, which is
#' then compare to reference distributions representing
#' primer design problems of different difficulties via the
#' total variance distance. The difficulty of the input primer design
#' problem is found by selecting the class of the 
#' reference distributions that has the smallest distance
#' to the estimated coverage distribution.
#' An estimate of the required number of primers to reach a given
#' \code{required.cvg} can be computed by setting
#' \code{primer.estimate} to \code{TRUE}. Since this estimate
#' is based solely on perfect matching primers, the number of
#' primers that would actually be required is typically less.
#'
#' @param template.df A \code{Templates} object providing
#' the template sequences for which the difficulty of designing
#' primers shall be estimated.
#' @param mode.directionality The directionality of the
#' primers that are to be designed. Either
#' \code{fw} for forward primers, \code{rev} for reverse primers,
#' or \code{both} for primers of both directions. By default,
#' both directions are considered.
#' @param primer.length A scalar numeric providing the 
#' target length of the designed primers. The default length 
#' of generated primers is set to \code{18}.
#' @param primer.estimate Whether the number of required primers shall be estimated. By default (\code{FALSE}), the number of required primers is not estimated.
#' @param required.cvg A scalar numeric in the range [0,1] providing the target coverage ratio for designing primers. 
#' The \code{required.cvg} is used only when \code{primer.estimate} is set to \code{TRUE} such that a solution to the set cover problem is required.
#' @return A list with the following fields:
#' \describe{
#' \item{\code{Classification}}{The estimated difficulty of the primer design task.}
#' \item{\code{Class-Distances}}{The total variance distance of the fitted
#' beta distribution to the reference distribution.}
#' \item{\code{Confidence}}{The confidence in the estimate of the
#' design tasks' difficulty as based on the class distances.}
#' \item{\code{Uncertain}}{Whether the classification is highly uncertain, that is
#' low-confidence.}
#' \item{\code{Nbr_primers_fw} and \code{Nbr_primers_rev}}{The respective number of 
#' required forward and reverse primers if \code{primer.estimate} was set to \code{TRUE}.}
#' }
#' @export
#' @examples
#' data(Ippolito)
#' design.estimate <- classify_design_problem(template.df)
#' # Estimate the number of required primers
#' design.estimate.nbr <- classify_design_problem(template.df, mode.directionality = "fw",
#'                          primer.length = 20, primer.estimate = TRUE)
classify_design_problem <- function(template.df, 
                                    mode.directionality = c("both", "fw", "rev"),
                                    primer.length = 18, 
                                    primer.estimate = FALSE,
                                    required.cvg = 1) {
    if (!is(template.df, "Templates")) {
        stop("Please supply a valid template data frame")
    }
    mode.directionality <- match.arg(mode.directionality)
    # 1. Estimate lower bound of possible primer coverage ratios
    cvg.df <- estimate.cvg(template.df, k = primer.length, mode.directionality)
    # 2. Fit a Beta distribution to the coverage ratios
    x <- c(cvg.df$fw$Coverage_Ratio, cvg.df$rev$Coverage_Ratio)
    fit.beta <- try(fitdistrplus::fitdist(x, "beta")) # beta has high error
    if (class(fit.beta) == "try-error") {
        # if there's too few unique x values, fitdistrplus won't be able to find a fit -> stop here!
        my.warning("ProblemEstimationProblem", "Could not estimate problem difficulty, problably because the estimated coverage distribution was too narrow.")
        return(NULL)
    }
    fit <- distr::Beta(fit.beta$estimate[1], fit.beta$estimate[2])
    #print("fit is: ")
    #print(fit.beta)
    #hist(rbeta(10000, shape1 = fit.beta$estimate[1], shape2 = fit.beta$estimate[2]))
    # 3. Compare the beta distribution to the reference distributions
    dists <- unlist(lapply(cvg.ref.dists, function(x) {
                    distrEx::TotalVarDist(x, fit)
    }))
    # 4. Classify
    best.idx <- which.min(dists)
    names(dists) <- names(cvg.ref.dists)
    classification <- names(cvg.ref.dists)[best.idx]
    # 5. Confidence of classification
    if (classification %in% c("very_easy", "easy")) {
        check.cols <- c("hard", "very_hard")
    } else if (classification %in% c("very_hard", "hard")) {
        check.cols <- c("easy", "very_easy")
    } else {
       check.cols <- c("easy", "very_easy", "hard", "very_hard") 
    }
    other.dist <- mean(dists[names(dists) %in% check.cols]) # mean of completely different distributions
    selected.dist <- dists[best.idx] # lowest distance
    confidence <- 1 - (selected.dist / other.dist)
    # refuse to give a classification if the standard error was too high, we're too uncertain of the real distribution in this case
    std.cutoff <- 10
    uncertain.class <- FALSE
    if (any(fit.beta$sd > std.cutoff)) {
        warning("Fit of distribution had a maximal standard error of: ",
                max(fit.beta$sd), ". Classification may be highly uncertain.")
        uncertain.class <- TRUE
    }
    nbr.fw.primers <- NA
    nbr.rev.primers <- NA
    if (primer.estimate) {
        # Identify the required nbr of primers
        primers.fw <- NULL
        primers.rev <- NULL
        if (mode.directionality == "fw") {
            primers.fw <- greedy.primers(cvg.df$fw, template.df, required.cvg)
        } else if (mode.directionality == "rev") {
            primers.rev <- greedy.primers(cvg.df$rev, template.df, required.cvg)
        } else {
            primers.fw <- greedy.primers(cvg.df$fw, template.df, required.cvg)
            # adjust rev requirement by forward cvg -> should capture the same templates!
            cvg.idx <- unique(unlist(covered.seqs.to.idx(primers.fw$Covered_Seqs, template.df)))
            required.nbr <- required.cvg * nrow(template.df)
            new.required.cvg <- min(length(cvg.idx) / required.nbr, 1)
            new.template.df <- template.df[cvg.idx[order(cvg.idx)],]
            primers.rev <- greedy.primers(cvg.df$rev, new.template.df, new.required.cvg)
        }
        if (length(primers.fw) != 0) {
            nbr.fw.primers <- nrow(primers.fw)
        } 
        if (length(primers.rev) != 0) {
            nbr.rev.primers <- nrow(primers.rev)
        } 
    }
    # Prepare output:
    result <- list("Classification" = classification,
                   "Class-Distances" = dists,
                   "Confidence" = confidence,
                   "Uncertain" = uncertain.class,
                   "Nbr_primers_fw" = nbr.fw.primers,
                   "Nbr_primers_rev" = nbr.rev.primers)

    return(result)
}
greedy.primers <- function(binding.df, template.df, required.cvg = 1) {
    # lower bound on the number of required primers
    # greedy takes quite some time ... improve code TODO
    o <- order(binding.df$Coverage_Ratio, decreasing = TRUE)
    binding.df <- binding.df[o,]
    selected.primers <- NULL
    cur.cvg <- 0
    while ((cur.cvg / nrow(template.df)) < required.cvg && nrow(binding.df) > 0) {
        best.covered <- as.numeric(unlist(strsplit(binding.df[1, "Covered_Seqs"], split = ",")))
        if (length(best.covered) == 0) {
            break
        }
        selected.primers <- my_rbind(selected.primers, binding.df[1,])
        marginal.cvg.gain <- length(best.covered)
        cur.cvg <- cur.cvg + marginal.cvg.gain
        binding.df <- binding.df[-1, ]  # remove the selected primer and selection candidates that were filtered this round
        binding.df <- evaluate.diff.primer.cvg(binding.df, best.covered, template.df)  
    }
    return(selected.primers)
}
#' Validates a Primers Object.
#'
#' Checks whether a Primers object is valid or not.
#'
#' @param object An input data frame to be checked for being a primer data frame.
#' @return \code{TRUE}, if the object is valid, FALSE otherwise.
#' @keywords internal
validate_primers <-  function(object) {
    # specify minimal set of columns that should be present in a primer data frame:
    required.fields <- list("Identifier" = c("factor"), "ID" = "factor", 
                   "Forward" = "character", "Reverse" = "character",
                   "primer_length_fw" = c("integer", "numeric"), 
                   "primer_length_rev" = c("integer", "numeric"), 
                   "Direction" = "character",
                   "Run" = "character")
    possible.fields <- NULL # don't check for additional fields here
    if (!is(object, "data.frame")) {
        return("Input was no data frame.")
    }
    check.fields <- check_setting(possible.fields, object, required.fields)
    if (check.fields) {
        # Check that "Run" is unique
        check.run <- length(unique(object$Run)) <= 1
        if (check.run) {
            return(TRUE)
        } else {
            msg <- "The 'Run' column may only contain a single unique value."
            return(msg)
        }
    } else {
        return(check.fields)
    }
}

#' Identification of Sequence Restriction Sites.
#'
#' Checks a set of primers for the presence of
#' restriction sites. To reduce the number of possible restriction sites,
#' only unambiguous restriction sites are taken into account and 
#' only common (typically used) restriction sites are checked if a common
#' restriction site can be found in a sequence.
#'
#' @param primer.df A \code{Primers} object containing the primer nucleotide
#' sequences to be checked for restriction sites.
#' @param template.df An object of class \code{Templates} containing
#' the templates corresponding to \code{primer.df}.
#' @param adapter.action The action to be performed when adapter sequences
#' are found. Either "warn" to issue a warning about adapter sequences or
#' "rm" to remove identified adapter sequences. Currently, only
#' the default setting ("warn") is supported.
#' @param selected Names of restriction sites that are to be checked.
#' By default \code{selected} is \code{NULL} in which case all REBASE 
#' restriction sites are taken into account.
#' @param only.confident.calls Whether only confident calls
#' of restriction sites are returned.
#' All restriction site call is considered \emph{confident} if the restriction site
#' is located in a region that does not match the template sequences.
#' Note that this classification requires that the provided primers
#' are somehow complementary to the provided templates.
#' In contrast, non-confident restriction site calls are 
#' based solely on the primer sequences and do not take the templates
#' into account, resulting in more false positive calls of restriction sites.
#' @param updateProgress A Shiny progress callback function. The default
#' is \code{NULL} meaning that no progress is tracked via the Shiny app.
#' @return A data frame with possible restriction sites found in every primer.
#' @references
#' Roberts, R.J., Vincze, T., Posfai, J., Macelis, D. (2010) REBASE–a database for DNA restriction
#' and modification: enzymes, genes and genomes. Nucl. Acids Res. 38: D234-D236. http://rebase.neb.com
#' @export
#' @family primer functions
#' @examples
#' data(Ippolito)
#' site.df <- check_restriction_sites(primer.df, template.df)
check_restriction_sites <- function(primer.df, template.df, 
                            adapter.action = c("warn", "rm"), 
                            selected = NULL, only.confident.calls = TRUE,
                            updateProgress = NULL) {
    if (length(template.df) == 0 || nrow(template.df) == 0 || 
        !is(template.df, "Templates")) {
        stop("Need correct templates to identify restriction sites.") 
    }
    if (length(primer.df) == 0 || nrow(primer.df) == 0 || 
        !is(primer.df, "Primers")) {
            stop("Need valid primers.")
    }
    adapter.action <- match.arg(adapter.action)
    fw.idx <- which(nchar(primer.df$Forward) != 0)
    rev.idx <- which(nchar(primer.df$Reverse) != 0)
    fw.sites <- NULL
    rev.sites <- NULL
    if (length(fw.idx) != 0) {
        seqs <- DNAStringSet(primer.df$Forward[fw.idx])
        names(seqs) <- primer.df$ID[fw.idx]
        primer.seqs <- DNAStringSet(seqs)
        template.seqs <- template.df$Sequence
        names(template.seqs) <- template.df$Identifier
        template.seqs <- DNAStringSet(template.seqs)
        fw.sites <- check_restriction_sites_single(primer.seqs, template.seqs, 
                    adapter.action, 
                    direction = "fw",
                    only.confident.calls = only.confident.calls,
                    updateProgress = updateProgress)
    }
    if (length(rev.idx) != 0) {
        seqs <- primer.df$Reverse[rev.idx]
        names(seqs) <- primer.df$ID[rev.idx]
        primer.seqs <- DNAStringSet(seqs)
        template.seqs <- Biostrings::reverseComplement(template.seqs)
        rev.sites <- check_restriction_sites_single(primer.seqs, template.seqs, 
                    adapter.action, 
                    direction = "rev", 
                    only.confident.calls = only.confident.calls,
                    updateProgress = updateProgress)
    }
    site.df <- rbind(fw.sites, rev.sites)
    return(site.df)
}
#' Update of Primer Binding Regions.
#'
#' Updates the relative primer binding sites in the templates
#' when the template binding regions have changed since the last
#' coverage computation.
#'
#' @param primer.df A \code{Primers} data frame.
#' @param template.df Templates with the new binding regions.
#' @param old.template.df Templates with the old binding regions.
#' @return A \code{Primers} object with updated relative binding positions.
#' @keywords internal
update_primer_binding_regions <- function(primer.df, template.df, old.template.df) {
    if (length(primer.df) == 0 || nrow(primer.df) == 0) {
        return(primer.df)
    }
    if (!is(primer.df, "Primers")) {
        stop("Please input a 'Primers' object.")
    }
    if (!"primer_coverage" %in% colnames(primer.df)) {
        warning("Cannot update primer binding regions as no coverage was available.")
        return(primer.df)
    }
    if (!is(template.df, "Templates") || !is(old.template.df, "Templates")) {
        stop("Need to have 'Templates' objects to update the primer regions.")
    }
    # check that new region annotations are for the same set of templates
    if (nrow(template.df) != nrow(old.template.df) || template.df$ID != old.template.df$ID) {
        # template sets do not match -> don't adjust anything
        return(primer.df)
    }
     # update the current primers object in rv_primers
    fw.region <- list(template.df$Allowed_Start_fw, template.df$Allowed_End_fw)
    rev.region <- list(template.df$Allowed_Start_rev, template.df$Allowed_End_rev)
    fw.region.old <- list(old.template.df$Allowed_Start_fw, old.template.df$Allowed_End_fw)
    rev.region.old <- list(old.template.df$Allowed_Start_rev, old.template.df$Allowed_End_rev)
    # update the relative binding regions of the primers according to the changes from new to old regions
    # for fw: 'end' defines the relatie position of interest
    fw.diff <- fw.region[[2]] - fw.region.old[[2]]
    # for reverse: start defines the position of relative interest
    rev.diff <- rev.region[[1]] - rev.region.old[[1]]
    dirs <- c("_fw", "_rev")
    # cols with relation to forward binding region:
    relative.cols.fw <- unlist(lapply(dirs, function(x) paste0("Relative_Forward_Binding_Position_", c("Start", "End"), x)))
    # this won't work if we have paired primers -> need to know covered_seqs_fw and covered_seqs_rev in this case
    template.ids <- lapply(strsplit(primer.df$Covered_Seqs, split = ","), as.numeric)
    # assume template order hasn't changed:
    template.idx <- lapply(template.ids, function(x) match(x, template.df$Identifier))
    for (i in seq_along(relative.cols.fw)) {
        col <- relative.cols.fw[i]
        vals <- lapply(strsplit(primer.df[, col], split = ","), as.numeric)
        len.ok <- unlist(lapply(seq_along(vals), function(x) length(vals[[x]]) == length(template.idx[[x]])))
        if (any(!len.ok) && any(primer.df[,col] != "")) {
            warning("Can't adjust binding range for paired primers at the moment.")
            return(primer.df)
        }
        # for adjustment, need to get the template-specific difference
        adj.vals <- unlist(lapply(seq_along(vals), function(x) paste(vals[[x]] - fw.diff[template.idx[[x]]], collapse = ",")))
        primer.df[, col] <- adj.vals
    }
    relative.cols.rev <- unlist(lapply(dirs, function(x) paste0("Relative_Reverse_Binding_Position_", c("Start", "End"), x)))
     for (i in seq_along(relative.cols.rev)) {
        col <- relative.cols.rev[i]
        vals <- lapply(strsplit(primer.df[, col], split = ","), as.numeric)
        len.ok <- unlist(lapply(seq_along(vals), function(x) length(vals[[x]]) == length(template.idx[[x]])))
        if (any(!len.ok) && any(primer.df[,col] != "")) {
            warning("Can't adjust binding range for paired primers at the moment.")
            return(primer.df)
        }
        # for adjustment, need to get the template-specific difference
        adj.vals <- unlist(lapply(seq_along(vals), function(x) paste(vals[[x]] + rev.diff[template.idx[[x]]], collapse = ",")))
        primer.df[, col] <- adj.vals
    }
    return(primer.df)
}
#' Identification of Restriction Sites.
#'
#' Identifies restriction sites in a list with putative
#' restriction sites provided by \code{bad.regions} using
#' a data frame of restriction sites given by \code{DB}.

#' @param bad.region IRanges with possible adapter sites.
#' @param DB A data frame with restriction enzyme sites.
#' @return A boolean data frame indicating the presence
#' of adapters for all considered restriction sites.
#' @keywords internal
restriction_hits <- function(bad.regions, DB) {
    bad.seqs <- unlist(lapply(bad.regions, function(x) as.character(unlist(x))))
    names(bad.seqs) <- unlist(lapply(seq_along(bad.regions), function(x) rep(names(bad.regions)[x], length(unlist(bad.regions[[x]])))))
    i <- NULL
    hit.db <- foreach (i = seq_len(nrow(DB)), .combine = cbind) %dopar% {
        enzyme <- as.character(DB$nam[i])
        temp <- Biostrings::vmatchPattern(DB$rep[i], bad.seqs, 
                max.mismatch = 0, with.indels = FALSE)
        hits <- sapply(temp, function(x) length(x) != 0)
        if (length(hits) == 0) {
            NULL
        } else {
            hit.idx <- which(hits)
            result <- data.frame("Hit" = hits, stringsAsFactors = FALSE) 
            colnames(result)[colnames(result) == "Hit"] <- enzyme
            result
        }
    }
    if (length(hit.db) != 0) {
        attr(hit.db, "ID") <- names(bad.seqs)
    }
    return(hit.db)
}

#' Identification of Badly Fitting Regions.
#'
#' Identify regions in the templates where the primers are
#' not very complementary. These regions indicate possible
#' restriction enzyme adapters.
#'
#' @param primer.seqs Primer sequences.
#' @param template.seqs Template sequences.
#' @param search.hits Template substrings that agree well
#' with the input primers.
#' @return A list with putative restriction sites for every primer.
#' @keywords internal
restriction_ali <- function(primer.seqs, template.seqs, search.hits) {
    bad.regions <- vector("list", length(search.hits))
    # define substitution matrix for alignment
    mat <- Biostrings::nucleotideSubstitutionMatrix(match = 1, mismatch = 0, 
                                                    baseOnly = FALSE, type = "DNA")
    for (i in seq_along(search.hits)) {
        # align primer with hit region in the templates
        if (is.na(search.hits[i])) {
            next
        }
        ali <- Biostrings::pairwiseAlignment(primer.seqs[i], search.hits[i],
                    substitutionMatrix = mat,
                     gapOpening = -3, gapExtension = -1, type = "global-local") # open up gaps directly
        # identify indels and non-matching regions
        # minimal number of bases of restriction 
        restrict.size <- 3
        if (length(ali@subject@indel[[1]]) != 0) {
            # there were indels -> look at the indel regions
            indel <- ali@subject@indel[[1]]
            # only select indels larger than cutoff
            sel <- indel@width >= restrict.size
            indel <- indel[sel]
            if (length(indel) != 0) {
                site <- Biostrings::extractAt(primer.seqs[i], indel)
                bad.regions[[i]] <- unname(site)
            }
        } else {
            # no indels -> look at the longest mismatch region
            mm.region <- unlist(Biostrings::mismatch(ali@pattern))
            # discard short mismatch regions:
            if (length(mm.region) > restrict.size) { 
                # more than 3 mismatches -> consider as adapter seq
                site <- Biostrings::extractAt(primer.seqs[i], 
                        IRanges(min(mm.region), max(mm.region)))
                bad.regions[[i]] <- unname(site)
            }
        }
    }
    names(bad.regions) <- names(primer.seqs)
    return(bad.regions)
}
#' Identification of Sequence Matches.
#'
#' Determines the most similar template sequence for every
#' input primer sequence. Used to identify regions for
#' alignment for the identification of restriction sites.
#'
#' @param primer.seqs Primer sequences.
#' @param template.seqs Template sequences.
#' @return A vector with the template regions matching the 
#' \code{primer.seqs} best.
#' @keywords internal
restriction_match <- function(primer.seqs, template.seqs) {
    search.hits <- rep(NA, length(primer.seqs))
    for (i in seq_along(primer.seqs)) {
        scores <- rep(Inf, length(template.seqs))
        best.ali <- NA
        for (j in seq_len(length(template.seqs))) { # n^2 * n^2 = n^4 performance ...
            # vmatchpattern doesn't support indels ...
            # max.mismatch: gives the detection limit for adapter sequences/overhang extent
            d <- Biostrings::matchPattern(as.character(primer.seqs[i]), 
                    as.character(template.seqs[j]), 
                    with.indels = TRUE, 
                    max.mismatch = 12)
            if (length(d@ranges) == 0) {
                # no hits
                next
            }
            # retrieve hits in the templates
            hits <- sapply(d, function(x) as.character(x))
            # determine best hit by computing distance between region of interest and pattern
            d <- stringdist::stringdist(hits, as.character(primer.seqs[i]))
            # select the hit with the smallest edit distance
            sel <- which.min(d)
            scores[j] <- d[sel]
            hit <- hits[sel]
            if (scores[j] <= min(scores)) {
                best.ali <- hit
            }
        }
        search.hits[i] <- best.ali
    }
    return(search.hits)
}
#' Identification of Sequence Restriction Sites.
#'
#' Checks the input sequences \code{seqs} for the presence of
#' restriction sites. By removing the restriction sites from a primer set, 
#' it is possible to identify the coverage of the primers
#' (e.g. using \code{\link{check_constraints}}) discounting for 
#' the impact of the mismatching bases caused by the insert.
#'
#' @param primer.seqs Nucleotide sequences of primers to be checked for restriction sites in terms of a \code{DNAStringSet} object.
#' @param template.seqs A \code{DNAStringSet} object with nucleotide sequences 
#' containing the templates corresponding to \code{seqs}.
#' @param adapter.action The action to be performed when adapter sequences
#' are found. Either "warn" to issue warning about adapter sequences or
#' "rm" to remove identified adapter sequences.
#' @param The primer direction that is checked.
#' @param selected Names of restriction sites that are to be checked.
#' By default \code{selected} is \code{NULL} in which case all REBASE 
#' restriction sites are checked.
#' @param only.confident.calls Only output confident calls of restriction sites.
#' @param updateProgress A Shiny progress callback function.
#' @return A data frame with restriction sites, if any could be found.
#' @references
#' Roberts, R.J., Vincze, T., Posfai, J., Macelis, D. (2010) REBASE–a database for DNA restriction
#' and modification: enzymes, genes and genomes. Nucl. Acids Res. 38: D234-D236. http://rebase.neb.com
#' @keywords internal
check_restriction_sites_single <- function(primer.seqs, template.seqs, adapter.action, 
                            direction = c("fw", "rev"),
                            selected = NULL, 
                            only.confident.calls = TRUE,
                            
                            updateProgress = NULL) {

    if (length(template.seqs) == 0 || length(primer.seqs) == 0) {
        return(NULL)
    }
    if (!is(primer.seqs, "DNAStringSet")) {
        stop("Primers should be a DNAStringSet.")
    }
    if (!is(template.seqs, "DNAStringSet")) {
        stop("Primers should be a DNAStringSet.")
    }
    if (length(direction) == 0) {
        stop("Please provide the 'direction' argument.")
    }
    direction <- match.arg(direction)
    #############
    # n.b.: this approach only works if the primer makes sense at all
    # if the primer isn't really complementary to the templates, 
    # we also won't be able to find restriction sites.
    # sysdata is available within the package automatically: enzdata
    #############
    # use only restriction sites without degeneracies for simplicity..
    counts <- sapply(regmatches(enzdata$site, gregexpr("[a-z]", enzdata$site, perl=TRUE)), length)
    sel <- which(counts == 0)
    DB <- enzdata[sel,]
    if (!is.null(selected)) {
        m <- match(selected, DB$nam)
        if (any(is.na(m))) {
            msg <- paste("Selected enzyme not available: ",
                paste(selected[is.na(m)], collapse = ", "), sep = "")
            stop(msg)
        }
        DB <- DB[m,]
    }
    if (is.function(updateProgress)) {
        detail <- "String matching"
        updateProgress(1/3, detail, "inc")
    }
    ##########
    # Confident calls: check for mismatches in the templates
    ##########
    # 1. Identify matches of primers in templates
    search.hits <- restriction_match(primer.seqs, template.seqs)
    if (is.function(updateProgress)) {
        detail <- "Aligning"
        updateProgress(1/3, detail, "inc")
    }
    # 2. Find the non-matching regions in the templates
    bad.regions <- restriction_ali(primer.seqs, template.seqs, search.hits)
    #####
    # Non-confident call: input primers
    #######
    input.regions <- vector("list", length(primer.seqs))
    names(input.regions) <- names(bad.regions)
    for (i in seq_along(primer.seqs)) {
        input.regions[[i]] <- Biostrings::DNAStringSetList(as.character(primer.seqs[i]))
    }
    all.regions <- list("Confident" = bad.regions, "Unconfident" = input.regions)
    if (is.function(updateProgress)) {
        detail <- "Identifying adapters"
        updateProgress(1/3, detail, "inc")
    }
    # 3. Identify whether the non-matching regions are restriction sites
    all.hits <- vector("list", length(all.regions))
    for (j in seq_along(all.regions)) {
        bad.regions <- all.regions[[j]]
        hit.db <- restriction_hits(bad.regions, DB) 
        if (is.null(hit.db)) {
            next
        }
        # 4. from all hits, only consider the common restriction sites
        common.DB <- DB[DB$common,]
        idx <- lapply(seq_len(nrow(hit.db)), function(x) which(as.logical(hit.db[x,])))
        hit.names <- lapply(idx, function(x) colnames(hit.db)[x])
        found.DB.idx <- lapply(seq_along(hit.names), function(x) idx[[x]][which(hit.names[[x]] %in% common.DB$nam)])
        no.hit.idx <- intersect(which(sapply(found.DB.idx, length) == 0), which(sapply(idx, length) != 0))
        if (length(no.hit.idx) != 0) {
            # no hit in the common restriction sites -> check for the most specific rare restriction site
            found.DB.idx[no.hit.idx] <- lapply(seq_along(no.hit.idx), function(x) idx[[no.hit.idx[x]]][which.max(nchar(DB$rep[idx[[no.hit.idx[x]]]]))])
        }
        if (length(found.DB.idx) == 0) {
            # still no hits? we're done
            next
        }
        names(found.DB.idx) <- attr(hit.db, "ID")
        hit.out <- vector("list", length(found.DB.idx))
        for (i in seq_along(found.DB.idx)) {
            if (length(found.DB.idx[[i]]) == 0) {
                # no hits found
                next
            }
            p.id <- names(found.DB.idx)[i]
            p.idx <- match(p.id, names(primer.seqs@ranges))
            out <- data.frame("Primer_Identifier" = p.id,
                              "Direction" = direction,
                              "Sequence" = unname(tolower(as.character(primer.seqs[p.idx]))),
                        "Enzyme" = as.character(DB$nam)[found.DB.idx[[i]]],
                        "Site" = tolower(DB$rep[found.DB.idx[[i]]]),
                        stringsAsFactors = FALSE)
            hit.out[[i]] <- out
        }
        hit.out <- do.call(rbind, hit.out)
        if (length(hit.out) != 0) {
            hit.out$Confidence <- names(all.regions)[j]
        }
        all.hits[[j]] <- hit.out
    }
    hit.out <- do.call(rbind, all.hits)
    # only return unique hits per primer
    hit.out <- plyr::ddply(hit.out, c("Primer_Identifier", "Sequence", "Enzyme", "Site"),
            function(x) plyr::arrange(x, substitute(Confidence))[1, ])
    if (only.confident.calls) {
        hit.out <- hit.out[hit.out$Confidence == "Confident",]
    }
    # warn about restriction sites:
    if (adapter.action == "warn") {
        # warn about restriction sites
        enzymes <- plyr::ddply(hit.out, c("Enzyme", "Site"), plyr::summarise, 
                    Identifiers = paste(substitute(Primer_Identifier), collapse = ","))
        if (nrow(enzymes) != 0) {
            out <- sapply(seq_len(nrow(enzymes)), function(x) 
                    paste("Found ", enzymes[x,"Enzyme"], " (", enzymes[x,"Site"], 
                    ") adapter in primers ", enzymes[x,"Identifiers"], 
                    ". Please check your sequences.", sep = ""))
            warning(paste(out, collapse = "\n"))
        }
    } else if (adapter.action == "rm") {
        # remove restriction sites from primers
        warning("Removal is not implemented yet.")
    }
    return(hit.out)
}
#' The Primers Class.
#'
#' The \code{Primers} class encapsulates a data frame
#' representing a set of primers. Objects of this class
#' store all properties associated with a set of primers,
#' for example the results from evaluating the properties
#' of a primer set or from determining its coverage.
#'
#' @section Basic columns:
#' In the following you can find a description of the most
#' important columns that can be found in objects of class \code{Primers}. 
#' Note that angular brackets indicate the existence of multiple possibilities.
#' The following columns are present when a set of primers
#' is loaded from a FASTA file using \code{\link{read_primers}}:
#' \describe{
#' \item{\code{ID}}{The identifiers of the primers.}
#' \item{\code{Identifier}}{The internal identifiers of the primers.}
#' \item{\code{Forward}}{The sequences of forward primers.}
#' \item{\code{Reverse}}{The sequences of reverse primers.}
#' \item{\code{primer_length<fw|rev>}}{The lengths of
#' forward and reverse primer sequences, respectively.}
#' \item{\code{Direction}}{Either 'fw' for forward primers,
#' 'rev' for reverse primers, or 'both' for a primer pair.}
#' \item{\code{Degeneracy_<fw|rev>}}{The degeneracy (ambiguity) of
#' forward and reverse primers, respectively.}
#' \item{\code{Run}}{An identifier describing the primer set.}
#' }
#'
#' @section Coverage-related columns:
#' The following columns are only available after primer coverage
#' has been computed, that is after \code{\link{check_constraints}}
#' has been called with the active \code{primer_coverage} constraint. Computed coverage
#' values relating solely to string matching are indicated by the prefix
#' \code{Basic_}, while columns without this prefix relate to the coverage after
#' applying the constraints formulated via \code{CoverageConstraints}.
#' Information on off-target coverage events are indicated by
#' the \code{Off_} prefix, while on-target coverage events do not carry
#' this prefix.
#' 
#' \describe{
#' \item{\code{primer_coverage}}{The number of templates that are
#' covered by the primers. Note that if a primer set contains
#' primers of both directions, a template is only considered covered
#' if it is covered by primers of both directions.}
#' \item{\code{Coverage_Ratio}}{The ratio of templates that are covered by the primers.}
#' \item{\code{Binding_Position_Start_<fw|rev>}}{The upstream position in the 
#' templates where forward and reverse primers respectively bind.}
#' \item{\code{Binding_Position_End_<fw|rev>}}{The downstream position in the templates where forward and reverse primers respectively bind.}
#' \item{\code{Relative_<Forward|Reverse>_Binding_Position_<Start|End>_<fw|rev>}}{
#' The binding upstream (\code{Start}) or downstream (\code{End}) positions 
#' of the primers relative to the forward (\code{Forward})
#' or reverse (\code{Reverse}) binding regions, either for 
#' forward (\code{fw}) or reverse primers (\code{rev}).}
#' \item{\code{Binding_Region_Allowed}}{Whether a coverage event
#' occurred in the target binding region or not. If the allowed
#' off-target ratio was set to 0 only coverage events within the 
#' the target region are reported.}
#' \item{\code{Nbr_of_mismatches_<fw|rev>}}{The number of mismatches
#' of forward and reverse primer coverage events, respectively.}
#' \item{\code{Mismatch_pos_<fw|rev>}}{The position of mismatches
#' for forward and reverse coverage events, respectively. Mismatch
#' positions are reported relative to the 3' end, that is, position
#' 1 indicates a mismatch in the last base of a primer.}
#' \item{\code{primer_specificity}}{The specificity of a primer
#' as determined by its ratio of off-target binding events.}
#' }
#' 
#' @section Constraint-related columns:
#' Each constraint that is considered when calling \code{\link{check_constraints}}
#' gives rise to at least one column in the provided \code{Primers} object.
#' Due to the large number of possible constraints, we will limit our description
#' to the \code{gc_clamp} constraint. Once the GC clamp property has been computed,
#' the \code{gc_clamp_fw} column contains the length of the GC clamp for forward 
#' primers and \code{gc_clamp_rev} the corresponding length for reverse primers.
#' Whether the desired extent of the GC clamp was obtained by a primer
#' is indicated by the \code{EVAL_gc_clamp} column. It contains \code{TRUE} when
#' the GC clamp constraint was fulfilled and \code{FALSE} when it was broken.
#' To identify whether all required constraints were fulfilled by a primer,
#' the \code{constraints_passed} column can be used. It contains \code{TRUE}
#' if all \code{active.constraints} used by \code{\link{check_constraints}} were fulfilled
#' and \code{FALSE} otherwise.
#'
#' @name Primers-class
#' @return A \code{Primers} object, an instance of a data frame.
#' @rdname Primers-class
#' @keywords Classes
#' @exportClass Primers
#' @export
#' @family primer functions
#' @seealso \code{\link{read_primers}} for loading a primer set,
#' \code{\link{score_degen}} for scoring the degeneracy of a primer,
#' \code{\link{primer_significance}} for determining the significance
#' of a primer set,
#' \code{\link{get_initial_primers}} for computing an initial set of primers,
#' \code{\link{design_primers}} for designing primer sets,
#' \code{\link{check_constraints}} for determining the properties of a primer set,
#' \code{\link{filter_primers}} for filtering a primer set,
#' \code{\link{check_restriction_sites}} to search for restriction sites,
#' \code{\link{get_cvg_ratio}} to determine the coverage ratio of a primer set,
#' \code{\link{create_report}} to create a PDF report for a primer set.
#'
#' @examples
#' primer.location <- system.file("extdata", "IMGT_data", "primers", "IGHV", 
#'                      "Ippolito2012.fasta", package = "openPrimeR")
#' primer.df <- read_primers(primer.location, "_fw", "_rev")
setClass("Primers",
   contains = c("data.frame"), validity = validate_primers)


setMethod("initialize", "Primers",
    function(.Object, ...) {
        # just call data frame constructor
        Object <- callNextMethod(.Object, ...)
    }
) 
#' @name Primers
#' @rdname Primers-class
#' @export
#' @param ... A data frame fulfilling the structural
#' requirements for initializing a \code{Primers} object.
Primers <- function(...) new("Primers", ...) 

my_format_df <- function (x) {
    if (is.atomic(x) && !is.null(x)) {
        stop("Internal structure doesn't seem to be a list. Possibly corrupt data.table.")
    }
    classes <- unlist(lapply(x, class))
    factor.idx <- which(classes == "factor")
    for (i in factor.idx) {
        x[, i] <- as.character(x[, i])
    }
    char.trunc <- function(x, trunc.char = 20) {
        trunc.char = max(0L, suppressWarnings(as.integer(trunc.char[1L])), na.rm=TRUE)
        if (!is.character(x) || trunc.char <= 0L) return(x)
        idx = which(nchar(x) > trunc.char)
        x[idx] = paste(substr(x[idx], 1L, as.integer(trunc.char)), "...", sep="")
        x
    }
    # truncate all strings
    classes <- unlist(lapply(x, class))
    string.idx <- which(classes == "character")
    for (i in string.idx) {
        x[, i] <- char.trunc(x[,i])
    }
    return(x)
}

my_show_df <- function(x, topn = 3,
                       nrows = 6) {
    # topn  - print the top topn and bottom topn rows with '---' inbetween 
    # nrows - if x is <= nrows, the whle table is printed
    if (!is.numeric(nrows)) nrows = 10
    if (!is.infinite(nrows)) nrows = as.integer(nrows)
    if (nrows <= 0L) return(invisible())   # ability to turn off printing
    if (!is.numeric(topn)) topn = 5
    topnmiss = missing(topn)
    topn = max(as.integer(topn),1L)
    if (nrow(x) == 0L) {
        if (length(x)==0L)
           cat("Empty object of class ", class(x), ".\n")
        else
           cat("Object of class ", class(x), " without any rows.")
        return(invisible())
    }
    if (topn*2<nrow(x) && (nrow(x)>nrows || !topnmiss)) {
        toprint = rbind(head(x, topn), tail(x, topn))
        rn = c(seq_len(topn), seq.int(to=nrow(x), length.out=topn))
        printdots = TRUE
    } else {
        toprint = x
        rn = seq_len(nrow(x))
        printdots = FALSE
    }
    toprint = my_format_df(x) # modify data frame to have only character classes 
    rownames(toprint) <- seq_len(nrow(toprint))
    if (is.null(names(x)) | all(names(x) == "")) colnames(toprint)=rep("", ncol(toprint))
    if (printdots) {
        toprint = rbind(head(toprint,topn),"---" = "",tail(toprint,topn))
        rownames(toprint) = format(rownames(toprint),justify="right")
        print(toprint,right=TRUE)
        return(invisible())
    }
    if (nrow(toprint)>20L)
        # repeat colnames at the bottom if more than 20 rows are printed
        toprint=rbind(toprint,matrix(colnames(toprint),nrow=1))
    print(toprint,right=TRUE)
    invisible()
}
setMethod("show", "Primers", function(object) {
    # overwrite the 'print' function using data.table's show
    my_show_df(asS3(object))
})
setMethod("summary", "Primers", function(object) {
    stats <- primer.set.parameter.stats(object, get.analysis.mode(object), NULL)
    return(stats)
})

#' cbind for Primers class.
#'
#' Ensures that the cbind result has the appropriate class.
#'
#' @param ... Parameters for cbind function.
#' @return Column binded Primers data frame.
#' @keywords internal
#' @export
#' @examples
#' data(Ippolito)
#' primer.df <- cbind(primer.df, primer.df)
cbind.Primers <- function(...) {
    df <- cbind.data.frame(...)
    df <- Primers(df)
    return(df)
}
#' S4 cbind for Primers.
#'
#' S4 cbind function for Primers.
#'
#' @export
#' @rdname Primers-method
#' @param x The Primers data frame.
#' @param y Another data frame.
#' @return Cbinded primer data frame.
#' @keywords internal
#' @examples
#' data(Ippolito)
#' primer.df <- cbind2(primer.df, seq_len(nrow(primer.df)))
setMethod("cbind2", c(x = "Primers", y = "ANY"),
    # need this in case of S3 dispatch fails.
    function(x, y, ...) {
        df <- cbind.data.frame(x, y, ...) # packages should not call .Internal()?
        df <- Primers(df)
        return(df)
    }
)
#' S4 rbind for Primers.
#'
#' S4 rbind function for Primers.
#'
#' @export
#' @return Rbinded primer data frame.
#' @rdname Primers-method
#' @keywords internal
setMethod("rbind2", c(x = "Primers", y = "ANY"),
    # need this in case of S3 dispatch fails.
    function(x, y, ...) {
        df <- rbind.data.frame(x, y, ...)
        df <- Primers(df)
        return(df)
    }
)
#' rbind for Primers class.
#'
#' Ensures that the rbind result has the appropriate class.
#'
#' @param ... Parameters for rbind function.
#' @return Row-binded Primers data frame.
#' @keywords internal
#' @export
#' @examples
#' data(Ippolito)
#' primer.df <- rbind(primer.df, primer.df)
rbind.Primers <- function(...) {
    df <- rbind.data.frame(...)
    df <- Primers(df)
    return(df)
}

#' Slicing Operator for Primers.
#'
#' Slices a Primers data frame.
#'
#' @param x The Primers data frame.
#' @param i The row index.
#' @param j The column index.
#' @param ... Other arguments to the slice operator.
#' @param drop Simplify data frame?
#' @exportMethod [
#' @rdname Primers-method
#' @return Subset of primer data frame.
#' @keywords internal
#' @aliases [,Primers-method
#' @examples
#' data(Ippolito)
#' primer.df <- primer.df[1:2,]
setMethod("[", c("Primers", "ANY", "ANY", "ANY"),
    function(x, i, j, ..., drop = TRUE) {
        if (missing(drop)) {
            df <- asS3(x)[i, j, ...]
        } else {
            df <- asS3(x)[i, j, ..., drop = drop]
        }
        if (class(df) == "data.frame") {
            # only set my class if we haven't simplified.
            #options("show.error.messages" = FALSE)
            p.df <- suppressWarnings(try(Primers(df), silent = TRUE))
            #options("show.error.messages" = TRUE)
            if (class(p.df) == "try-error") {
                # removed crucial columns -> turn into data frame
                p.df <- asS3(df)
            } 
            df <- p.df
        }
        return(df)
    }
)
#' Dollar Operator for Primers.
#'
#' Stores data in a column of a Primers data frame.
#'
#' @exportMethod $<-
#' @rdname Primers-method
#' @param name The name of the column.
#' @param value The values of the column.
#' @return Primer data frame with replaced column.
#' @keywords internal
#' @examples
#' data(Ippolito)
#' primer.df$Forward[1] <- "ctagcgggaccg"
setMethod("$<-", "Primers", 
    function(x, name, value) {
        df <- asS3(x)
        eval(parse(text = paste("df$", name, " <- value", sep = "")))
        df <- Primers(df)
        return(df)
    }
)
#' Updates the Primer Coverage.
#'
#' Updates the most important columns in a primer data frame according to the selected
#' coverage definition. Only coverage events with less or equal than the allowed number of
#' mismatches according to the selected coverage definition will be retained.
#'
#' @param primer.df A \code{Primers} object.
#' @param template.df A \code{Templates} object.
#' @param allowed.mismatches A numeric giving the maximal number of allowed.mismatches.
#' @param cvg.definition The definition of coverage to be used, either "constrained" or "basic".
#' @return A primer data frame with modified coverage information.
#' @keywords internal
update_primer_cvg <- function(primer.df, template.df, allowed.mismatches, cvg.definition = c("constrained", "basic")) {
    # note: only updates primer_coverage, covered_seqs, not binding positions etc.
    # Consider only coverage events with <= allowed.mismatches and according to the provided cvg.definition
    if (length(primer.df) == 0 || nrow(primer.df) == 0) {
        # nothing to update ..
        return(primer.df)
    }
    cvg.definition <- match.arg(cvg.definition)
    mm.info <- prepare_mm_plot(primer.df, template.df)
    # select events according to selected cvg.definition and number of mismatches only
    df <- mm.info[mm.info$Coverage_Type == cvg.definition & mm.info$Number_of_mismatches <= allowed.mismatches, ]
    # select one event per primer-template pair (worst-case terminal mismatch)
    ddf <- plyr::ddply(df, c("Primer", "Template"), plyr::summarize,
                        Position = unique(substitute(Position_3terminus)), 
                        Number_of_mismatches = unique(substitute(Number_of_mismatches)))
    #if (length(ddf) == 0 || nrow(ddf) == 0) {
        ## no coverage events to update available ..
        ## TODO: if no coverage events are found we have to set everything to 0...
        #return(primer.df)
    #}
    cvd <- plyr::ddply(ddf, "Primer", plyr::summarize, Covered_Seqs = paste(substitute(Template), collapse = ","))
	#print(cvd) -> Primer is 1,2,3,4 ... instead of its character rep in windows
    m <- match(cvd$Primer, primer.df$ID)
    new.df <- data.frame(ID = primer.df$ID, Covered_Seqs = "", primer_coverage = 0,
                         stringsAsFactors = FALSE)
    if (length(m) != 0) { # we found coverage events found for selection
        # set corresponding entries in new.df to updated values
        cvd.identifier <- unlist(lapply(strsplit(cvd$Covered_Seqs, split = ","), function(x) paste(template.df$Identifier[match(x, template.df$ID)], collapse = ",")))
        new.df[m, 'Covered_Seqs'] <- cvd.identifier
        new.df[m, 'primer_coverage'] <- sapply(strsplit(cvd.identifier, split = ","), length)
    }
    # update primer data frame with new values
    new.primer.df <- primer.df
    new.primer.df[, colnames(new.df)] <- new.df
    return(new.primer.df)
}

#' Plot of Primer Binding Regions.
#'
#' Visualizes the number of binding events of the input primers
#' with respect to the allowed binding regions in the templates.
#'
#' @param primers Either a single \code{Primers} object or 
#' a list with \code{Primers} objects.
#' @param templates If \code{primers} is a \code{primers} object,
#' please supply a \code{Templates} object.
#' If \code{primers} is a list, please supply a corresponding list of
#' \code{Templates} objects.
#' @param direction The directionality of primers to be plotted. This can either
#' be "both" to plot primers of both directions (the default), "fw" to plot
#' only forward primers, or "rev" to plot only reverse primers.
#' @param group Optional identifiers of template groups for which binding events should
#' be determined. By default, \code{group} is set to \code{NULL} such that
#' all templates are considered.
#' @param relation An optional character vector specifying whether binding region data shall
#' be plotted relative to the forward (\code{fw}) or reverse (\code{rev}) 
#' target binding regions.
#' @param region.names An optional, two-component character vector specifying
#' the identifiers for the primer binding region and the amplified region.
#' @param ... \code{highlight.set} (the identifiers of primer sets
#' to be highlighted, if \code{primers} is a list)
#' @return A plot for primer binding region comparison.
#' @export
#' @include templates.R
#' @family coverage visualizations
#' @examples
#' # Primer binding regions of a single primer set
#' data(Ippolito)
#' plot_primer_binding_regions(primer.df, template.df)
#' # Primer binding regions of multiple primer sets
#' data(Comparison)
#' plot_primer_binding_regions(primer.data[1:3], template.data[1:3])
setGeneric("plot_primer_binding_regions", 
    function(primers, templates, direction = c("both", "fw", "rev"), 
        group = NULL, relation = c("fw", "rev"), 
        region.names = c("Binding region", "Amplification region"), ...) {
        standardGeneric("plot_primer_binding_regions")
})
#' Plot of Primer Binding Regions for a Single Primer Set.
#'
#' Plots the primer binding regions in the templates.
#'
#' @param primers An object of class \code{Primers} with annotated
#' primer coverage.
#' @param templates An object of class \code{Templates} providing the
#' template sequences corresponding to \code{primers}. 
#' @param direction Primer direction
#' @param group The template groups for which binding events should
#' be determined. By default, \code{group} is set to \code{NULL} such that
#' all templates are considered.
#' @param relation A character vector specifying whether binding region data shall
#' be plotted relative to the forward (\code{fw}) or reverse (\code{rev}) 
#' target binding regions.
#' @param region.names Names for the primer binding region and the amplified region.
#' @return A histogram of primer binding regions.
#' @keywords internal
setMethod("plot_primer_binding_regions", 
    methods::signature(primers = "Primers", templates = "Templates"),
    function(primers, templates,
        direction = c("both", "fw", "rev"), group = NULL, 
        relation = c("fw", "rev"), 
        region.names = c("Binding region", "Amplification region")) {
        
        relation <- match.arg(relation)
        direction <- match.arg(direction)
        if (length(primers) == 0 && nrow(primers) == 0) {
            return(NULL)
        }
        if (!is(primers, "Primers")) {
            stop("Please input a valid primers object.")
        }
        primer.data <- list(primers)
        names(primer.data) <- unique(primers$Run)
        p <- plot_primer_binding_regions(primer.data, list(templates),
                direction, group, relation, region.names = region.names)
        return(p)
})
create_region_boxes <- function(primers, templates, relation, region.names, ymin, ymax, xmax) {
    regions <- vector("list", length(templates)) # binding regions
    if (relation == "fw") {
        for (i in seq_along(templates)) {
            template.df <- templates[[i]]
            region.extents <- template.df$Allowed_End_fw_initial - template.df$Allowed_Start_fw_initial
            idx <- which.max(region.extents)
            s <- template.df$Allowed_Start_fw_initial[idx]
            e <- template.df$Allowed_End_fw_initial[idx]
            e.now <- template.df$Allowed_End_fw[idx]
            s.now <- template.df$Allowed_Start_fw[idx]
            delta <- e.now - e
            rel.binding.start <- -(e.now - s.now) - 1
            rect.x.leader <- -((e-s) + delta) - 1
            rect.x.leader.end <- rect.x.leader + e - s
            region.info.l <- data.frame("Region" = region.names[1], "xmin" = rect.x.leader, "xmax" = rect.x.leader.end, "ymin" = ymin, "ymax" = ymax)
            region.info.t <- data.frame("Region" = region.names[2], "xmin" = rect.x.leader.end + 1, "xmax" = xmax, "ymin" = ymin, "ymax" = ymax)
            region.df <- rbind(region.info.l, region.info.t)
            region.df$Run <- unique(primers[i]$Run)
            region.df$RelStartPosition <- rel.binding.start
            regions[[i]] <- region.df 
            #binding.starts[[i]] <- data.frame("Run" = names(primer.data)[i], "Position" = rel.binding.start)
        }
    } else {
        for (i in seq_along(templates)) {
            template.df <- templates[[i]]
            region.extents <- template.df$Allowed_End_rev_initial - template.df$Allowed_Start_rev_initial
            idx <- which.max(region.extents)
            seq.len <- nchar(template.df$Sequence[idx])
            s <- seq.len - template.df$Allowed_Start_rev_initial[idx] + 1
            e <- seq.len - template.df$Allowed_End_rev_initial[idx] + 1
            e.now <- seq.len - template.df$Allowed_End_rev[idx] + 1
            s.now <- seq.len - template.df$Allowed_Start_rev[idx] + 1
            # positive delta: binding region has to be shifted back (negative), negative delta: binding region has to be shifted forward (positive)
            delta <- s - s.now
            rel.binding.start <- e.now - s.now - 1
            ########
             rect.x.leader <- -((e-s) + delta) - 1
            rect.x.leader.end <- rect.x.leader + e - s
            #########
            rect.x.leader <- -(s - e + delta) - 1
            rect.x.leader.end <- rect.x.leader - e + s
            region.info.l <- data.frame("Region" = region.names[1], "xmin" = rect.x.leader, "xmax" = rect.x.leader.end, "ymin" = ymin, "ymax" = ymax)
            region.info.t <- data.frame("Region" = region.names[2], "xmin" = rect.x.leader.end + 1, "xmax" = xmax, "ymin" = ymin, "ymax" = ymax)
            region.df <- rbind(region.info.l, region.info.t)
            region.df$Run <- unique(primers[i]$Run)
            region.df$RelStartPosition = rel.binding.start
            regions[[i]] <- region.df 
            #binding.starts[[i]] <- data.frame("Run" = names(primer.data)[i], "RelStartPosition" = rel.binding.start)
        }
    }
    region.df <- do.call(rbind, regions)
    return(region.df)
}
#' Plot of Primer Binding Regions for Multiple Sets.
#'
#' Plots the primer binding regions for every primer set.
#'
#' @param primers List with primer data frames.
#' @param templates List with template data frames.
#' @param direction Direction of primers.
#' @param group Template groups to plot.
#' This defaults to plotting all groups.
#' @param relation Plot binding region relative to forward binding region or reverse?
#' @param region.names Names for the primer binding region and the amplified region.
#' @param highlight.set Primer sets to highlight in the plot.
#' @return A plot for primer binding region comparison.
#' @keywords internal
setMethod("plot_primer_binding_regions", 
    methods::signature(primers = "list", templates = "list"),
    function(primers, templates,
        direction = c("both", "fw", "rev"), group = NULL, 
        relation = c("fw", "rev"), 
        region.names = c("Binding region", "Amplification region"), 
        highlight.set = NULL) {

    if (length(region.names) != 2) {
        stop("Need 2 region names.")
    }
    xlab <- "Location"
    direction <- match.arg(direction)
    relation <- match.arg(relation)
    xlab <- "Binding position"
    plot.data <- lapply(seq_along(primers), function(x) primer.binding.regions.data(primers[[x]], 
        templates[[x]], direction, group, relation))
    # annotate with run info
    runs <- get.run.names(primers)
    plot.data <- lapply(seq_along(plot.data), function(x) if (length(plot.data[[x]]) != 
        0 && nrow(plot.data[[x]]) != 0) {
        cbind(plot.data[[x]], Run = rep(primers[[x]]$Run[1], nrow(plot.data[[x]])))
    })
    plot.df <- do.call(rbind, plot.data)
    if (class(plot.df) != "data.frame") {
        return(NULL)
    }
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
    # new annotation function for rectangle plot
    plot.df <- plyr::ddply(plot.df, c("ID", "Start", "End", "Run"), plyr::summarize, Count = length(substitute(ID)))
    # cumulate counts for multiple primer IDs -> ymax value for rectangles
    plot.df$ymax <- stats::ave(plot.df$Count, plot.df$Start, plot.df$End, plot.df$Run, FUN = cumsum)
    # ymin: cumulative sum minus the current count
    plot.df$ymin <- plot.df$ymax - plot.df$Count
    colnames(plot.df)[colnames(plot.df) %in% c("Start", "End")] <- c("xmin", "xmax")
    #print(plot.df)
    levels <- unique(plot.df$Run)
    plot.df$Run <- factor(plot.df$Run, levels = levels[order(levels)])
    bwidth <- 10 # cover 10 positions with one bar
    ymin <- min(c(-1, -0.1 * max(plot.df$ymax)))
    ymax <- -0.2 # slightly negative to retain the border of the box
    # consider widest range of allowed regions from all template sets:
    ############################
    # set boundaries of the plot:
    xmin <- min(plot.df$xmin)
    xmax <- max(plot.df$xmax)
    ## require a minimal length for the region segments
    if (xmin > 50) {
        xmin <- -50
    }
    if (xmax < 50) {
        xmax <- 50
    }
    region.df <- create_region_boxes(primers, templates, relation, region.names, ymin, ymax, xmax)
    xmin <- min(xmin, min(region.df$xmin))
    x.ticks <- pretty(seq(xmin, xmax))
    x.labels <- x.ticks
    x.labels[x.labels > 0] <- paste0("+", x.labels[x.labels > 0])
    rect.colors <- c("#e5f4ff", "#ffefe5")
    names(rect.colors) <- region.names
    r.colors <- rep(rect.colors, length(levels(plot.df$Run)))
    plot.df$ID <- abbreviate(plot.df$ID, getOption("openPrimeR.plot_abbrev")) # shorten primer IDs
    pal <- getOption("openPrimeR.plot_colors")["Primer"] # the RColorBrewer palette to use
    primer.colors <- colorRampPalette(brewer.pal(8, pal))(length(unique(plot.df[, "ID"])))
    names(primer.colors) <- unique(plot.df[, "ID"])
    primer.colors <- c(rect.colors, primer.colors)
    p <- ggplot(plot.df) + ylab("Number of coverage events") + 
        xlab(xlab) + 
        ggtitle("Sites of primer binding in the templates")  +
        theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
         scale_x_continuous(limits = c(xmin, xmax),
                           breaks = x.ticks,
                           labels = x.labels) + 
        geom_vline(xintercept = -1, colour = "red") + # end of binding region
        geom_vline(data = region.df, aes_string(xintercept = "RelStartPosition"), colour = "red") + # start of binding region
        # rectangles for histogram of binding events:
        geom_rect(data = plot.df, 
                  mapping = aes_string(xmin = "xmin", xmax = "xmax",
                                      ymin = "ymin", ymax = "ymax",
                                      fill = "ID"),
                    linetype = "blank", alpha = 0.35) + # rectangles show some overlap due to the transparency!!
       scale_fill_manual(values = primer.colors, breaks = unique(plot.df$ID)) + 
        # x-axis rectangles to annotate binding/amplification region:
        geom_rect(data = region.df, 
            mapping = aes_string(xmin="xmin", xmax="xmax", 
                        ymin="ymin", ymax="ymax", fill = "Region"),
            alpha = 0.5,
            colour = "#3d3835", 
            size = 0.3, show.legend = FALSE)
    if (length(unique(plot.df$Run)) > 1) {
        # don't show facets and don't show individual primer legend
        p <- p + facet_wrap(~Run, ncol = 3) +
            guides(fill = FALSE)
    } else {
        # only show rectangle text for single plot
        p <- p + geom_text(data=region.df, 
            aes_string(x = "xmin+(xmax-xmin)/2", 
                       y = "ymin+(ymax-ymin)/2", 
                       label = "Region"), size = 4)

    }
    if (length(unique(plot.df$ID)) > 15) {
        # don't show legend for many primers
        p <- p + guides(fill = FALSE)
    }
    if (length(highlight.set) != 0) {
        # highlight selected sets
        # highlight of strip.text (facet part) individually doesn't work
        # -> this is the solution i came up with
        highlights <- data.frame(Run = highlight.set)
        p <- p + geom_rect(data=highlights,aes(xmin=-Inf, xmax=Inf, 
                    ymin=-Inf, ymax=Inf), fill='red', alpha=0.1)
    }
    return(p)
})

#' Input of Primers.
#'
#' Reads one or multiple input files with primer sequences. The input can either be in FASTA
#' or in CSV format.
#'
#' The input arguments \code{fw.id}, \code{rev.id}, \code{merge.ambig}, and \code{max.degen} are only used for loading primers from a FASTA file.
#' If you want to load a FASTA file, please ensure that \code{fw.id} and 
#' \code{rev.id} are set according to the keywords indicating
#' the primer directionalities in the FASTA file.
#' When loading a CSV file, the format of the file should adhere
#' to the structure defined by the \code{\link{Primers}} class. 
#' You can easily store a \code{Primers} objects as a CSV file using the
#' \code{\link{write_primers}} function.
#' 
#' @param primer.location Path to a single or multiple primer FASTA or CSV files.
#' @param fw.id For FASTA input, the identifier for forward primers in the FASTA headers.
#' @param rev.id For FASTA input, the identifier for reverse primers in the FASTA headers.
#' @param merge.ambig Indicates whether similar primers should be merged
#' ("merge") using IUPAC ambiguity codes or whether primers should 
#' be disambiguated ("unmerge").
#' By default \code{merge.ambig} is set to "none", leaving primers as they are.
#' @param max.degen A scalar numeric providing the maximum
#' allowed degeneracy for merging primers if
#' \code{merge.ambig} is set to "merge".
#' Degeneracy is defined by the number of disambiguated sequences 
#' that are represented by a degenerate primer.
#' @param template.df An object of class \code{Templates}.
#' If \code{template.df} is provided the primers are checked
#' for restriction sites upon input. By default
#' \code{template.df} is \code{NULL} such that the primers are not
#' checked for restriction sites.
#' @param adapter.action The action to be performed when \code{template.df} is
#' provided for identifying adapter sequences.
#' Either "warn" to issue warning about adapter sequences or
#' "rm" to remove identified adapter sequences. The default is "warn".
#' @param sample.name An identifier for the input primers. 
#' @param updateProgress A Shiny progress callback function. This is
#' \code{NULL} by default such that no progress is tracked.
#' @return An object of class \code{Primers} for a single 
#' \code{primer.location} or a list of such objects for multiple locations.
#' @export
#' @keywords Primers
#' @examples
#' primer.fasta <- system.file("extdata", "IMGT_data", "primers", "IGHV", 
#'                      "Ippolito2012.fasta", package = "openPrimeR")
#' primer.df <- read_primers(primer.fasta, "_fw", "_rev")
#' # Read multiple FASTA files
#' fasta.files <- list.files(system.file("extdata", "IMGT_data", "primers", 
#'                  "IGHV", package = "openPrimeR"), pattern = "*\\.fasta",
#'                  full.names = TRUE)[1:3]
#' primer.data <- read_primers(fasta.files)
#' # Read primers from a CSV file
#' primer.csv <- system.file("extdata", "IMGT_data", "comparison", 
#'              "primer_sets", "IGL", "IGL_openPrimeR2017.csv",  package = "openPrimeR")
#' primer.df <- read_primers(primer.csv)
#' # Read multiple primer CSV files
#' primer.files <- list.files(path = system.file("extdata", "IMGT_data", "comparison", 
#'                          "primer_sets", "IGH", package = "openPrimeR"),
#'                           pattern = "*\\.csv", full.names = TRUE)[1:3]
#' primer.data <- read_primers(primer.files)
#' # Read a mixture of FASTA/CSV files:
#' mixed.primers <- c(primer.fasta, primer.csv)
#' primer.data <- read_primers(mixed.primers)
read_primers <- function(primer.location, fw.id = "_fw", rev.id = "_rev", 
                    merge.ambig = c("none", "merge", "unmerge"), 
                    max.degen = 16, template.df = NULL, 
                    adapter.action = c("warn", "rm"), sample.name = NULL,
                    updateProgress = NULL) {

    adapter.action <- match.arg(adapter.action)
    if (is.function(updateProgress)) {
        detail <- "Reading primers"
        updateProgress(1/2, detail, "inc")
    }
    if (length(primer.location) > 1) {
        # load multiple primer FASTA/CSV files
        primers <- read_primers_multiple(primer.location, fw.id = fw.id, 
                    rev.id = rev.id, merge.ambig = merge.ambig, max.degen = max.degen, 
                    template.df = template.df,
                    adapter.action = adapter.action, sample.name = sample.name,
                    updateProgress = updateProgress) 
    } else {
        # load a single primer set from FASTA/CSV
        primers <- read_primers_single(primer.location, fw.id = fw.id, rev.id = rev.id,
                    merge.ambig = merge.ambig, max.degen = max.degen, 
                    template.df = template.df,
                    adapter.action = adapter.action, sample.name = sample.name,
                    updateProgress = updateProgress) 
    }
    return(primers)
}
read_primers_single <- function(primer.location, fw.id = "_fw", rev.id = "_rev", 
                    merge.ambig = c("none", "merge", "unmerge"), 
                    max.degen = 16, template.df = NULL, 
                    adapter.action = c("warn", "rm"), sample.name = NULL,
                    updateProgress = NULL) {

    if (!file.exists(primer.location)) {
        stop(paste("No file found at specified location: ",
                primer.location, sep = ""))
    }
    # load the primers and at the same time determine whether it's FASTA/CSV input
    primers <- try(my.read.fasta(primer.location, 
                tolower(names(Biostrings::IUPAC_CODE_MAP))), silent = TRUE)
    if (class(primers) == "try-error") {
        # let's try to read as CSV
        fasta.error <- attr(primers, "condition")
        primers <- try(read_primers_csv(primer.location), silent = TRUE)
        if (class(primers) == "try-error") {
            csv.error <- attr(primers, "condition")
            stop(paste("Could not read the file as FASTA or CSV.",
                "The FASTA error was: ",
                fasta.error, ".\n The CSV error was: ", csv.error))
       } else {
            # it's a CSV:
            ext <- "csv"
       }
    } else {
        # it's a FASTA:
        ext <- "fasta"
    }
    if (ext == "csv") {
        # nothing to change in the input
    } else if (ext == "fasta") {
        # sanitize the input
        if (length(primers) == 0) {
            return(NULL)
        }
        if (is.function(updateProgress)) {
            detail <- "Annotating primers"
            updateProgress(1/2, detail, "inc")
        }
        primer.seqs <- sapply(primers, function(x) paste(tolower(x), collapse = ""))
        headers <- sapply(primers, function(x) attr(x, "Annot"))
        if (length(sample.name) == 0) { 
            # use fasta name without file ending as sample name
            sample.name <- sub("^([^.]*).*", "\\1", basename(primer.location))
        }
        # create Primers object:
        primers <- read_primers.internal(primer.seqs, headers, fw.id, rev.id, 
                        merge.ambig, max.degen, sample.name)
    } else {
        stop("Unknown filetype.")
    }
    if (!is.null(template.df)) {
        # check for restriction sites if template.df was inputted
        sites <- check_restriction_sites(primers, template.df, 
                    adapter.action, updateProgress = updateProgress)
    }
    return(primers)
}
#' Input of Multiple Primer Sets.
#'
#' Reads multiple CSV files representing stored objects of class \code{Primers}.
#'
#' @param filenames The paths to multiple primer CSV/FASTA files.
#' @return A list containing objects of class \code{Primers}.
#' @keywords internal
read_primers_multiple <- function(filenames, fw.id, rev.id, merge.ambig,
                    max.degen, template.df, adapter.action, 
                    sample.name, updateProgress) {
    data <- lapply(filenames, function(x) read_primers_single(x, fw.id = fw.id, 
                    rev.id = rev.id, merge.ambig = merge.ambig,
                    max.degen = max.degen, template.df = template.df, 
                    adapter.action = adapter.action, sample.name = sample.name,
                    updateProgress = updateProgress))
    # annotate names of list/make 'Run' identifiers unique:
    data <- set.run.names(data, filenames)
    return(data)
}
set.run.names <- function(data, filenames = NULL) {
    # annotate names of list using set run identifiers
    if (length(data) == 0) {
        return(data)
    }
    if (length(filenames) != 0) {
        # filenames are available -> use stripped filename for empty files
        run.names <- sub("^([^.]*).*", "\\1", basename(filenames))
    } else {
        run.names <- rep("Unknown", length(data))
    }
    non.empty.idx <- which(unlist(lapply(data, nrow)) != 0)
    # use stored run identifier for non-empty files:
    run.names[non.empty.idx] <- sapply(data[non.empty.idx], function(x) x$Run[1])
    # ensure that run names are unique
    run.names <- uniqtag::uniqtag(run.names, k = Inf, make.unique)
    names(data) <- run.names
    # update run names in the loaded primer sets that aren't empty
    for (i in seq_along(non.empty.idx)) {
        data[[non.empty.idx[i]]]$Run <- run.names[non.empty.idx[i]]
    }
    return(data)
}
#' Read Primer CSV File.
#'
#' Reads a primer data frame stored in a CSV file.
#'
#' @param file The path to a csv file containing the primer data.
#' @return A \code{Primers} object, an instance of a data frame.
#' @keywords internal
read_primers_csv <- function(file) { 
    # turn some columns from integer to char if necessary 
    # (column with a single number will be interpreted as
    # numeric otherwise) ...
    if (!file.exists(file)) {
        stop(paste("Primer CSV file at '", file, "' could not be found."), sep = "")
    }
    # character columns:
    fix.cvg.cols.constrained <- c("Relative_Forward_Binding_Position_Start_fw", "Relative_Forward_Binding_Position_End_fw", 
        "Relative_Forward_Binding_Position_Start_rev", "Relative_Forward_Binding_Position_End_rev", 
        "Relative_Reverse_Binding_Position_Start_fw", "Relative_Reverse_Binding_Position_End_fw", 
        "Relative_Reverse_Binding_Position_Start_rev", "Relative_Reverse_Binding_Position_End_rev", 
        "Covered_Seqs", "Binding_Position_Start_rev", "Binding_Position_End_rev", 
        "Binding_Position_Start_fw", "Binding_Position_End_fw", "Binding_Region_Allowed", 
        "Nbr_of_mismatches_fw", "Nbr_of_mismatches_rev", "Binding_Region_Allowed_fw", 
        "Mismatch_pos_fw", "Mismatch_pos_rev", "Binding_Region_Allowed_rev")
    off.fix.cvg.cols.constrained <- paste0("Off_", fix.cvg.cols.constrained)
    fix.cvg.cols.constrained <- c(fix.cvg.cols.constrained, off.fix.cvg.cols.constrained)
    fix.cvg.cols.basic <- paste0("Basic_", fix.cvg.cols.constrained)
    fix.cvg.cols <- c(fix.cvg.cols.constrained, fix.cvg.cols.basic)
    fix.cols <- c(fix.cvg.cols,
                   # coverage constraints:
                   c("primer_efficiency", "T_EVAL_primer_efficiency", 
                     "terminal_mismatch_pos", "annealing_DeltaG",
                     "stop_codon", "substitution", "coverage_model",
                  "Forward", "Reverse", "constraints_passed_T", "Direction", "Run"))
    numeric.fix.cols <- c("primer_length_fw", "primer_length_rev",
                        "mean_primer_efficiency", "primer_specificity")
    factor.fix.cols <- c("ID", "Identifier")
    p <- try(read.csv(file, stringsAsFactors = FALSE, row.names = NULL), silent = TRUE)
    if (class(p) == "try-error") {
        msg <- paste0("The csv file: '", file, "' does not seem to represent a valid object of class 'Primers'",
                      " Please ensure the correct format of the file!")
        my.error("TemplateFormatIncorrect", msg)
    }
    for (i in 1:ncol(p)) {
        if (colnames(p)[i] %in% fix.cols) {
            p[, colnames(p)[i]] <- as.character(p[, colnames(p)[i]])
            na.idx <- which(is.na(p[, colnames(p)[i]]))
            if (length(na.idx) != 0) {
                # replace NA's in char cols with empty string:
                p[na.idx, colnames(p)[i]] <- ""
            }
        } else if (colnames(p)[i] %in% numeric.fix.cols) {
            p[, colnames(p)[i]] <- as.numeric(p[, colnames(p)[i]])
        }  else if (colnames(p)[i] %in% factor.fix.cols) {
            p[, colnames(p)[i]] <- factor(p[, colnames(p)[[i]]], 
                                    levels = p[, colnames(p)[[i]]])
        }
    }
    p <- try(Primers(p), silent = TRUE)
    if (class(p) == "try-error") {
        my.error("TemplateFormatIncorrect", 
            paste0("The loaded primer csv data from the file '",
                   file, "' did not represent a valid 'Primers' object.",
                   " Please check your input!"))
    }
    return(p)
}


#' Internal Function for Reading Primers
#'
#' Reads the given primer sequences into a data frame.
#'
#' @param primer.seqs Primer sequences.
#' @param headers Headers of the primer FASTA file.
#' @param fw.id Identifier of forward primers in the headers.
#' @param rev.id Identifier of reverse primers in the headers.
#' @param merge.ambig Should ambiguous primers be merged?
#' @param max.degen Maximum allowed degeneracy
#' @return A data frame with primer sequences.
#' @keywords internal
read_primers.internal <- function(primer.seqs, headers, fw.id, rev.id, merge.ambig = c("none", "merge", "unmerge"), 
    max.degen, sample.name) {
   
    merge.ambig <- match.arg(merge.ambig)
    mode <- "default"
    if (length(primer.seqs) == 0) {
        return(NULL)
    }
    if (fw.id == "" && rev.id == "") {
        # assume only fw primers are present
        fw.id <- "_fw"
        rev.id <- NA
        headers <- paste(headers, fw.id, sep = "")
        my.warning("NotifyPrimersMissingKeyword", "No keywords provided for primer directionalities. Primers were assumed to be sense primers.")
    } else if (fw.id == "" && rev.id != "") {
        # rev id is given -> assume complement of rev hits are forward primers
        my.warning("NotifyPrimersMissingKeyword", "A keyword was provided only for reverse primers. All primers without the reverse keyword were assumed to be forward primers.")
        mode <- "guess_forward"
        fw.id <- NA
    } else if (fw.id != "" && rev.id == "") {
        my.warning("NotifyPrimersMissingKeyword", "A keyword was provided only for forward primers. All primers without the forward keyword were assumed to be reverse primers.")
        mode <- "guess_reverse"
        rev.id <- NA
    }
    dirs <- rep(NA, length(headers))
    dir.fw.idx <- NULL
    dir.rev.idx <- NULL
    if (mode != "guess_forward") {
        dir.fw.idx <- grep(fw.id, headers, fixed = TRUE)
    } else {
        dir.fw.idx <- setdiff(seq_along(headers), grep(rev.id, headers, fixed = TRUE))
    }    
    if (mode != "guess_reverse") {
        dir.rev.idx <- grep(rev.id, headers, fixed = TRUE)
    } else {
        dir.rev.idx <- setdiff(seq_along(headers), grep(fw.id, headers, fixed = TRUE))
    }
    if (length(intersect(dir.fw.idx, dir.rev.idx) != 0)) {
        msg <- "Direction annotation is not distinct. Please use distinct keywords."
        my.error("NotifyPrimersNoDirection", msg)
    }
    dirs[dir.fw.idx] <- "Forward" # names of columns
    dirs[dir.rev.idx] <- "Reverse"
    use.idx <- seq_along(headers)
    if (any(is.na(dirs))) {
        warning.msg <- "Not all primer directions are known. Check your input primer fasta file for the correct direction keywords. Primers without known direction are ignored. The affected primers are:"
        p <- paste(headers[which(is.na(dirs))], collapse = ",")
        msg <- paste(warning.msg, p, sep = "")
        my.error("NotifyPrimersNoDirection", msg)
    }
    p.df <- data.frame(Header = headers, Direction = dirs, Sequence = primer.seqs, 
        Seq_Length = nchar(primer.seqs), stringsAsFactors = FALSE)[use.idx, ]
    rownames(p.df) <- NULL
    # try to match identifiers of primers
    h <- headers
    if (!is.na(fw.id)) {
        h <- sub(fw.id, "", headers[use.idx], fixed = TRUE)
    }
    if (!is.na(rev.id)) {
        h <- sub(rev.id, "", h, fixed = TRUE)
    }
    h <- sub(">", "", h)
    p.df$ID <- h
    # check that every sample has non-duplicate fw/rev primer
    unique.directions <- ddply(p.df, "ID", summarize, duplicate = any(duplicated(substitute(Direction))))
    if (any(unique.directions$duplicate)) {
        # duplicated directions found
        idx <- which(unique.directions$duplicate)
        warn <- paste("Found primers with duplicate directions! The following primers were affected:\n", 
            paste(unique.directions$ID[idx], collapse = "\n", sep = ""))
        my.error("PrimersDuplicateDirections", warn)
    }
    # n.b.: dcast changes variable order -> change back again
    ori.order <- p.df$ID
    p.df <- dcast(p.df, ID ~ Direction, value.var = "Sequence")  # wide-format :-)
    m <- match(p.df$ID, ori.order)
    p.df <- p.df[order(m),]
    # replace missing primers with ''
    p.df[is.na(p.df)] <- ""
    # if one direction is not present -> set it
    dir.mode <- "both"
    if (!"Forward" %in% colnames(p.df)) {
        dir.mode <- "rev"
        p.df$Forward <- ""
    }
    if (!"Reverse" %in% colnames(p.df)) {
        dir.mode <- "fw"
        p.df$Reverse <- ""
    }
    if (merge.ambig == "merge") {
        p.df <- merge.ambig.primers(p.df, dir.mode, max.degen)
        # message(p.df) disambiguate primers with ambiguities
    } else if (merge.ambig == "unmerge") {
        p.df <- disambiguate.primers(p.df)
    }
    # annotate direction of primers:
    directions <- ifelse(nchar(p.df$Forward) != 0 & nchar(p.df$Reverse) != 0, "both", 
                    ifelse(nchar(p.df$Forward) != 0, "fw", "rev"))
    p.df <- cbind(Identifier = paste(seq_len(nrow(p.df)), directions, sep = ""), 
        ID = p.df$ID, p.df[, c("Forward", "Reverse")], 
        primer_length_fw = nchar(p.df$Forward), 
		primer_length_rev = nchar(p.df$Reverse), 
        Direction = directions,
		Degeneracy_fw =  score_degen(strsplit(p.df$Forward, split = "")),
		Degeneracy_rev = score_degen(strsplit(p.df$Reverse, split = "")),
        Run = sample.name,
        stringsAsFactors = FALSE)
    sel <- which(p.df$Direction == "fw") 
    if (length(sel) != 0 && !is.na(fw.id)) {
        p.df$ID[sel] <- paste0(p.df$ID[sel], fw.id)
    }
    sel <- which(p.df$Direction == "rev") 
    if (length(sel) != 0 && !is.na(rev.id)) {
        p.df$ID[sel] <- paste0(p.df$ID[sel], rev.id)
    }
    p.df$ID <- factor(p.df$ID, levels = p.df$ID)
    p.df$Identifier <- factor(p.df$Identifier, levels = p.df$Identifier)
    p.df <- Primers(p.df) # create Primers object
    return(p.df)
}
#' Disambiguation of Primers.
#'
#' Disambiguates ambiguous primer sequences into all possible sequences.
#'
#' @param p.df Primer data frame.
#' @return Data frame with disambiguated primers.
#' @keywords internal
disambiguate.primers <- function(p.df) {
    f <- convert.from.iupac(p.df$Forward)
    r <- convert.from.iupac(p.df$Reverse)
    # combine
    combis <- lapply(seq_along(f), function(x) expand.grid(f[[x]], r[[x]], stringsAsFactors = FALSE))
    identifiers <- unlist(lapply(seq_along(combis), function(x) paste(p.df$ID[x], 
        "_", 1:(nrow(combis[[x]])), sep = "")))
    df <- do.call(rbind, combis)
    new.df <- data.frame(ID = identifiers, Forward = df[, 1], Reverse = df[, 2], 
        stringsAsFactors = FALSE)
    # message(new.df) message(class(new.df$Forward))
    return(new.df)
}
#' Pairing of Forward and Reverse Primers.
#'
#' Pairs forward and reverse primers such that coverage is maximized 
#' for every pair.
#'
#' @param primer.df An object of class \code{Primers}.
#' @return An object of class \code{Primers} containing the paired primers.
#' @keywords internal
pair_primers <- function(primer.df, template.df) {
    fw.idx <- which(primer.df$Forward != "")
    rev.idx <- which(primer.df$Reverse != "")
    if (length(fw.idx) == 0 || length(rev.idx) == 0) {
        return(primer.df) # nothing to pair
    }
    if (!"primer_coverage" %in% colnames(primer.df)) {
        stop("pair_primers requires primer coverage.")
    }
    # don't need to pair already paired primers ...
    both.df <- primer.df[primer.df$Direction == "both",]
    # try to pair all single forward primers
    fw.s.idx <- which(primer.df$Direction == "fw")
    rev.s.idx <- which(primer.df$Direction == "rev")
    # combinations of fw and rev primers
    combis <- expand.grid(fw.s.idx, rev.s.idx)
    # determine intersection between all combinations
    cvd <- covered.seqs.to.idx(primer.df[, "Covered_Seqs"], template.df)
    cvd.shared <- lapply(seq_len(nrow(combis)), 
        function(x) intersect(cvd[[combis[x,1]]], cvd[[combis[x,2]]])
    )
    cvg.shared <- sapply(cvd.shared, length)
    combis$Coverage <- cvg.shared
    sel.combis <- which(combis$Coverage != 0)
    combis <- combis[sel.combis,]
    # construct a new primer data frame using the pairs
    pair.df <- primer.df[combis[,1],]
    pair.df$Reverse <- primer.df$Reverse[combis[,2]]
    pair.df$primer_coverage <- combis$Coverage
    pair.df$Coverage_Ratio <- pair.df$primer_coverage / nrow(template.df)
    pair.df$Covered_Seqs <- unlist(lapply(cvd.shared[sel.combis], function(x) paste(template.df$Identifier[x], collapse = ",")))
    pair.df$Direction <- "both"
    pair.df$ID <- factor(paste0(abbreviate(primer.df$Identifier[combis[,1]], 5), "+", abbreviate(primer.df$Identifier[combis[,2]], 5)))
    cvg.matrix <- get.coverage.matrix(pair.df, template.df)
    sel.cols <- remove.redundant.cols(seq_len(nrow(pair.df)), cvg.matrix)
    pair.df <- pair.df[sel.cols,]
    # output the already paired primers as well as the new pairs
    out.df <- my_rbind(both.df, pair.df)
    return(out.df)
}

#' Storing Primers to Disk.
#'
#' Writes a set of primers to disk, either as a FASTA or CSV file.
#'
#' @param primer.df An object of class \code{Primers} to be stored to disk.
#' @param fname The path to the file where the primers should be stored.
#' @param ftype A character vector giving the type of the file.
#' This can either be "FASTA" or "CSV" (default: "FASTA").
#' @return Stores primers to \code{fname}.
#' @export
#' @keywords Primers
#' @family primer functions
#' @examples
#' data(Ippolito)
#' # Store primers as FASTA
#' fname.fasta <- tempfile("my_primers", fileext = ".fasta")
#' write_primers(primer.df, fname.fasta)
#' # Store primers as CSV
#' fname.csv <- tempfile("my_primers", fileext = ".csv")
#' write_primers(primer.df, fname.csv, "CSV")
write_primers <- function(primer.df, fname, ftype = c("FASTA", "CSV")) {
    ftype <- match.arg(ftype)
    if (ftype == "FASTA") {
        fw.idx <- which(primer.df$Forward != "")
        rev.idx <- which(primer.df$Reverse != "")
        seqs <- c(primer.df$Forward[fw.idx], primer.df$Reverse[rev.idx])
        labels <- c(as.character(primer.df$ID[fw.idx]), as.character(primer.df$ID[rev.idx]))
        seqinr::write.fasta(as.list(seqs), as.list(labels), fname)
    } else if (ftype == "CSV") {
        write.csv(primer.df, file = fname, row.names = FALSE)
    } else {
        warning("Unknown filetype: ", ftype)
    }
}
#' Direction of Primers.
#'
#' Identifies the directionality of the input primers.
#'
#' @param primers A primer data frame.
#' @return \code{both} if both, forward and reverse primers exist in \code{primers}.
#' Otherwise, if either only forward primers or reverse primers exist, returns \code{fw}
#' or \code{rev}, respectively.
#' @keywords internal
get.analysis.mode <- function(primers) {
    if (length(primers) == 0) {
        return(NULL)
    }
    mode <- "both"
    if (all(primers$Forward == "" | is.na(primers$Forward))) {
        mode <- "rev"
    } else if (all(primers$Reverse == "" | is.na(primers$Reverse))) {
        mode <- "fw"
    }
    return(mode)
}
#' Determination of the Ratio of Covered Templates. 
#'
#' Determines the ratio of template sequences 
#' that are covered by the evaluated input primers. The ratio
#' is in the interval [0,1] where 0 indicates 0\% coverage (no templates
#' covered) and 1 indicates 100\% coverage (all templates covered).
#'
#' The manner in which the coverage ratio is evaluated depends on 
#' the directionality of the input primers.
#' If either only forward or reverse primers are inputted, the individual
#' coverage of each primer is used to determine the overall coverage. 
#' If, however, forward and reverse primers are inputted at the same time, 
#' the coverage is defined by the intersection of binding events from both,
#' forward and reverse primers.
#'
#' @param primers A \code{Primers} object containing the primers
#' for which the coverage should be evaluated.
#' @param template.df A \code{Templates} object containing
#' the template sequences corresponding to \code{primers}. 
#' @param allowed.mismatches The number of allowed mismatches
#' for determining the coverage of the templates. By default,
#' \code{allowed.mismatches} is set to \code{NULL} such that
#' all annotated coverage events are considered.
#' @param cvg.definition Whether to output the coverage obtained
#' from string matching ("basic") or the expected coverage ("constrained"),
#' which is constructed by applying the coverage constraints. By default,
#' \code{cvg.definition} is set to "constrained".
#' @param mode.directionality If \code{mode.directionality} is provided,
#' the coverage of templates is computed for a specific direction of primers.
#' Either "fw" (forward coverage only), "rev" (reverse coverage only), or "both" for both directions. By default, \code{mode.directionality} is \code{NULL}
#' such that the directionality of the primers is determined automatically.
#' @param as.char Whether the coverage ratio should
#' be outputted as a percentage-formatted character vector. By default,
#' \code{as.char} is set to \code{FALSE} such that a numeric is returned.
#' @return By default, a numeric providing the expected primer coverage ratio. 
#' If \code{as.char} is \code{TRUE}, the output is provided as a
#' percentage-formatted character vector. 
#' The attributes \code{no_covered}, \code{no_templates}, 
#' and \code{covered_templates} provide the number of covered templates, 
#' the total number of templates, and the IDs of covered templates, respectively.
#' @export
#' @family Primers
#' @examples
#' data(Ippolito)
#' cvg.ratio <- get_cvg_ratio(primer.df, template.df)
#' # determine the identitity coverage ratio 
#' cvg.ratio.0 <- get_cvg_ratio(primer.df, template.df, allowed.mismatches = 0)
get_cvg_ratio <- function(primers, template.df, allowed.mismatches = NULL, 
                        cvg.definition = c("constrained", "basic"), 
                        mode.directionality = NULL, as.char = FALSE) {
    if (length(primers) == 0 || nrow(primers) == 0) {
        # no coverage if there's no input.
        return(0.0)
    }
    if (!is(primers, "Primers")) {
        stop("Please input a primer data frame.")
    }
    if (!is(template.df, "Templates")) {
        stop("Please provide a valid template data frame.")
    }
    if (!"Covered_Seqs" %in% colnames(primers) || length(template.df) == 0) {
        return(NA)
    }
    cvg.definition <- match.arg(cvg.definition)
    # check whether primers and templates correspond to one another
    if (!check_correspondence(primers, template.df)) {
        msg <- paste0("Could not match all coverage events to the input template data frame.\n",
            "Please verify whether the primer coverage was actually computed for the input template data frame.")
        stop(msg)
    }
    if (length(mode.directionality) == 0) {
        # determine directionality
        mode.directionality <- get.analysis.mode(primers)
    } else {
        # check input mode
        if (!mode.directionality %in% c("fw", "rev", "both")) {
            stop("Unknown mode directionality for cvg_ratio")
        }
    }
    if (cvg.definition == "constrained" && length(allowed.mismatches) == 0) {
        # use the vanilla coverage definition: faster and doesn't require mismatch information
        cvd <- get_covered.vanilla(primers, template.df, mode.directionality)
    } else {
        primer.df <- get_template_cvg_data(primers, template.df)
        # select only coverage events from the selected cvg definition
        primer.df <- primer.df[primer.df$Status == cvg.definition, ]
        # get only the unique events for every template (for both: fw & primer are found here)
        primer.df <- plyr::ddply(primer.df, c("Template", "Status"), function(x) plyr::arrange(x, substitute(Number_of_mismatches))[1,])
        # retrieve only the allowed coverage events
        if (length(allowed.mismatches) != 0) {
            # select only events with less or equal the allowed number of mismatches
            primer.df <- primer.df[primer.df$Number_of_mismatches <= allowed.mismatches,]
        } 
        cvd <- unique(primer.df$Template)
    }
    max.cvg <- length(cvd) / nrow(template.df)  # maximum theoretical coverage
    if (as.char) {
        max.cvg <- paste0(round(max.cvg *100, 2), "%")
    }
    # add number of templates covered and total number of templates as attributes
    attr(max.cvg, "no_covered") <- length(cvd)
    attr(max.cvg, "no_templates") <- nrow(template.df)
    attr(max.cvg, "covered_templates") <- as.character(cvd)
    return(max.cvg)
}

#' Determination of the Covered Sequences.
#'
#' Determines the covered template sequences given by \code{template.df} 
#' that are covered by the primers given by \code{primers}.
#'
#' The manner in which the coverage ratio is evaluated depends on 
#' the directionality of the input primers.
#' If either only forward or reverse primers are inputted, the individual
#' coverage of each primer is used to determine the overall coverage. If, however,
#' forward and reverse primers are inputted at the same time, 
#' the coverage is defined by the intersection of binding events from both,
#' forward and reverse primers.
#'
#' @param primers A \code{Primers} object containing the primers
#' for which the coverage should be evaluated.
#' @param template.df A \code{Templates} object containing
#' the template sequences corresponding to \code{primers}. 
#' @param mode.directionality If \code{mode.directionality} is provided,
#' the coverage of templates is computed for a specific direction of primers.
#' Either "fw" (forward coverage only), "rev" (reverse coverage only), or "both" for both directions. If \code{mode.directionality} is not provided
#' the direction is determined by the input primers.
#' @return The IDs of all covered templates.
#' @keywords internal
get_covered.vanilla <- function(primers, template.df, mode.directionality = NULL) {
    if (length(primers) == 0 || nrow(primers) == 0) {
        # no coverage if there's no input.
        return(0.0)
    }
    if (!is(primers, "Primers")) {
        stop("Please input a primer data frame.")
    }
    if (!is(template.df, "Templates")) {
        stop("Please provide a valid template data frame.")
    }
    if (!"Covered_Seqs" %in% colnames(primers) || length(template.df) == 0) {
        return(NA)
    }
    # check whether primers and templates correspond to one another
    if (!check_correspondence(primers, template.df)) {
        msg <- paste0("Could not match all coverage events to the input template data frame.\n",
            "Please verify whether the primer coverage was actually computed for the input template data frame.")
        stop(msg)
    }
    fw.idx <- which(primers$Forward != "")
    rev.idx <- which(primers$Reverse != "")
    if (length(mode.directionality) == 0) {
        # determine directionality
        mode.directionality <- get.analysis.mode(primers)
    } else {
        # check input mode
        if (!mode.directionality %in% c("fw", "rev", "both")) {
            stop("Unknown mode directionality for cvg_ratio")
        }
    }
    primer.df <- primers
    if (mode.directionality == "both") {
        # fw and rev primers -> intersection of fw and rev cvg events gives the total
        # coverage
        hits.fw <- unique(unlist(covered.seqs.to.idx(primer.df$Covered_Seqs[fw.idx], 
            template.df)))
        hits.rev <- unique(unlist(covered.seqs.to.idx(primer.df$Covered_Seqs[rev.idx], 
            template.df)))
        cvd <- intersect(hits.fw, hits.rev)  # idx of covered templates
    } else if (mode.directionality == "fw") {
        # only fw primers direction
        cvd <- unique(unlist(covered.seqs.to.idx(primer.df$Covered_Seqs[fw.idx], template.df)))
    } else {
        cvd <- unique(unlist(covered.seqs.to.idx(primer.df$Covered_Seqs[rev.idx], template.df)))
    }
    # convert to template IDs
    cvd <- template.df$ID[cvd]
    return(cvd)
}

