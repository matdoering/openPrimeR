#######
# Use alignments+trees to initalize primers
#######

#' Multiple Sequence Alignment
#'
#' Computes a multiple sequence alignment using MAFFT.
#'
#' @param seqs The sequences to be aligned.
#' @param names The identifiers of the sequences.
#'
#' @return An alignment object as created from seqinr's read.alignment method.
#' @keywords internal
#' @references 
#' Katoh, Misawa, Kuma, Miyata 2002 (Nucleic Acids Res. 30:3059-3066) 
#' MAFFT: a novel method for rapid multiple sequence alignment based
#' on fast Fourier transform. 
align.seqs <- function(seqs, names) {
    runID <- digest::digest(seqs, "md5")
    file.loc <- tempfile(paste0("in_ali_templates_", runID), fileext = ".fasta")
    seqinr::write.fasta(as.list(seqs), file.out = file.loc, names = as.list(names))
    out.loc <- tempfile(paste0("out_ali_templates_", runID), fileext = ".fasta")
    call <- paste("mafft", " --auto --quiet ", file.loc, " > ", out.loc, sep = "")
    if (Sys.which("mafft") == "") {
        stop("Cannot align sequences since MAFFT is not installed on your system.")
    }
    system(call)  #, ignore.stderr = TRUE) # cannot be called with suppressed output ..
    seq.ali <- try(seqinr::read.alignment(out.loc, format = "fasta"), silent = FALSE)
    if (class(seq.ali) == "try-error") {
        # alignment not possible
        warning("MAFFT wasn't able to align the templates. Returning.")
        return(NULL)
    }
	# clean up temporary files
	on.exit({
		file.remove(c(file.loc, out.loc))
	})
    return(seq.ali)
}
#' Identication of Short Primers.
#'
#' Identify initial primers that are too short.
#'
#' @param primer.candidates Primer alignment.
#' @param min.len Minimal primer length.
#' @return The index of proposed primers that are shorter than \code{min.len}.
#' @keywords internal
prefilter.primer.candidates <- function(primer.candidates, min.len) {
    n.gaps <- apply(as.character(primer.candidates), 1, function(x) length(which(x == 
        "-")))
    len <- ncol(primer.candidates) - n.gaps  # ungapped length
    rm.idx <- which(len < min.len)
    return(rm.idx)
}
#' Filtering of Primer Candidates
#'
#' Filters primer candidates according to length and duplications.
#"
#' @param primer.candidates Alignment of candidate primers.
#' @param min.len Minimal required length of primers.
#' @return Filtered alignment of candidate primers.
#' @keywords internal
filter.primer.candidates <- function(primer.candidates, min.len) {
    # remove entries in the primer candidate matrix whose ungapped sequences are
    # shorter than min.len
    
    # idea: count gaps
    rm.idx <- prefilter.primer.candidates(primer.candidates, min.len)
    if (length(rm.idx) != 0) {
        # message(paste('Removed ',length(rm.idx), ' candidates that were shorter than ',
        # min.len, ' bases.', sep = ''))
        primer.candidates <- primer.candidates[-rm.idx, ]
    }
    # remove duplicated candidates: we lose the implicit cvg info from the tree with
    # this, but we are faster.
    if (nrow(primer.candidates) != 0) {
        al <- ape::as.alignment(primer.candidates)
        rm.idx <- which(duplicated(al$seq))
        if (length(rm.idx) != 0) {
            al$nb <- al$nb - length(rm.idx)
            al$seq <- al$seq[-rm.idx]
            al$nam <- al$nam[-rm.idx]
        }
        primer.candidates <- ape::as.DNAbin(al)
    }
    return(primer.candidates)
}

#' Primer Degeneration Score.
#'
#' Determines the degeneration score of a sequence.
#'
#' The degeneration of an ambiguous sequence is defined as the number
#' of unambiguous sequences that the ambiguous sequence represents.
#' Let a sequence \code{S} of length \code{n} be represented by a collection of sets such that
#' \deqn{S = {s_1, s_2, \ldots, s_n}}
#' where \eqn{s_i} indicates the set of unambiguous bases found
#' at position \eqn{i} of the primer. Then the degeneracy \eqn{D} of a primer
#' can be defined as 
#' \deqn{D = \prod_i{|s_i|}}
#' where \eqn{|s_i|} provides the number of disambiguated bases at position \eqn{i}.
#'
#' @param seq A list of vectors with single characters.
#' @param gap.char The gap character in sequences.
#' @return The number of unambiguous sequences represented by \code{seq}.
#' @family primer functions
#' @export
score_degen <- function(seq, gap.char = "-") {
    if (length(seq) == 0) {
        return(NA)
    }
    if (class(seq) != "list") {
        stop("'seq' should be a list of split characters.")
    }
    gapped.codemap <- Biostrings::IUPAC_CODE_MAP
	if (length(gap.char) != 0 && nchar(gap.char) == 1) {
		gapped.codemap[gap.char] <- gap.char
	}
    degen <- unlist(lapply(seq, function(x) prod(nchar(gapped.codemap[toupper(x)]))))
    return(degen)
}
#' Tree Ancestry
#'
#' Checks whether \code{ancestor.node} is an ancestor to the nodes 
#' specified in \code{test.node}.
#'
#' @param tree The phylogenetic tree to be tested.
#' @param ancestor.node A node to be checked for being an ancestor to \code{test.node}.
#' @param test.node Possible descendants of \code{ancestor.node}.
#' @return TRUE, if \code{ancestor.node} is an ancestor to any node in \code{test.node}.
#' @keywords internal
ancestor_of <- function(tree, ancestor.node, test.node) {
    e <- tree$edge
    return(any(test.node %in% e[e[, 1] == ancestor.node, 2]))
}
#' Computation of Consensus.
#'
#' Computes the consensus of the sequences in the input alignment.
#'
#' @param ali An alignment object.
#' @return A consensus sequence without gap characters.
#' @keywords internal
get.consensus.seq <- function(ali) {
    c.seq <- seqinr::consensus(ali, method = "IUPAC")
    na.idx <- which(is.na(c.seq))
    if (length(na.idx) != 0) {
        c.seq.m <- seqinr::consensus(ali[, na.idx, drop = FALSE], method = "majority")
        # a) replace gaps in the consensus with gaps if they are the majority
        gap.idx <- which(c.seq.m == "-")
        c.seq[na.idx][gap.idx] <- c.seq.m[gap.idx]
        if (length(gap.idx) != 0) {
            na.idx <- na.idx[-gap.idx]
        }
        # b) form consensus for remaining positions where non-gap characters are the
        # majority, excluding gaps from consensus
        if (length(na.idx) != 0) {
            non.gap.idx <- apply(ali[, na.idx, drop = FALSE], 2, function(y) y != 
                "-")
            c.seq.m <- sapply(seq_len(ncol(non.gap.idx)), function(x) seqinr::consensus(ali[non.gap.idx[, 
                x], na.idx[x], drop = FALSE], method = "IUPAC"))
            c.seq[na.idx] <- c.seq.m
        }
    }
    # remove gaps from consensus
    c.seq <- c.seq[c.seq != "-"]
    return(c.seq)
}
#' Determine Tree Consensus Sequences
#'
#' Creates all possible consensus sequences from a phylogenetic tree.
#'
#' Ambiguous sequences are only generated with a degeneracy of at most \code{max.degen}. 
#' The tree is iterated from leaves to the top, i.e., starting from least degeneracy to most degeneracy. 
#' Merges only take place when the degeneracy of the resulting sequence would
#' be at most \code{max.degen}. Gaps are removed from the alignments.
#'
#' @param tree The phylogenetic tree.
#' @param max.degen The maximal degeneration of consensus primers.
#' @param primer.candidates Alignment of primers.
#' @return Data frame with consensus primers extracted from the tree.
#' @keywords internal
get.tree.seqs <- function(tree, max.degen, primer.candidates) {
    if (length(tree) == 0) {
        # too few seqs to build tree
        return(NULL)
    }
    N.tips <- length(tree$tip.label)
    selected.seqs <- NULL
    node.ignore.list <- NULL
    degeneracy.count <- NULL
    primer.IDs <- NULL
    merge.idx <- NULL
    for (i in seq_len(tree$Nnode)) {
        # message(i)
        node <- N.tips + tree$Nnode - i + 1
        stree <- ape::extract.clade(tree, node)  # get tree from current node
        ######### 
        # plot(stree) # visualize subtree 
        # nodelabels()
        # tiplabels() 
        ##############
        # ignore nodes that are
        # ancestors of trees that already exceeded the degeneracy count
        if (ancestor_of(tree, node, node.ignore.list)) {
            node.ignore.list <- c(node.ignore.list, node)
            next
        }
        seq.ids <- stree$tip.label  # seqs of this tree
        m <- match(seq.ids, rownames(primer.candidates))
        if (any(is.na(m))) {
            stop(paste("Integrity error: tree did not contain expected ID: ", seq.ids[is.na(m)], 
                " (hashing problem).", sep = ""))
        }
        ali <- primer.candidates[m, ]
        c.seq <- get.consensus.seq(as.character(ali))
        degen <- score_degen(list(c.seq))
        if (degen <= max.degen) {
            selected.seqs <- c(selected.seqs, paste(c.seq, collapse = ""))
            degeneracy.count <- c(degeneracy.count, degen)
            id <- paste(seq.ids, collapse = ",")
            primer.IDs <- c(primer.IDs, id)
            merge.idx <- c(merge.idx, list(m[order(m)]))
        } else {
            # message(paste('ignored: degeneracy was ', degen, sep = ''))
            node.ignore.list <- c(node.ignore.list, node)
        }
    }
    if (length(selected.seqs) != 0) {
        node.result <- data.frame(ID = paste(primer.IDs, sep = ""), Sequence = selected.seqs, 
            Degeneracy = degeneracy.count, Merge_Idx = sapply(merge.idx, function(x) paste(x, 
                collapse = ",", sep = "")), stringsAsFactors = FALSE)
        # remove duplicates
        node.result <- node.result[!duplicated(node.result$Sequence), ]
    } else {
        node.result <- NULL
    }
    return(node.result)
}
#' Ranges for Initial Primers.
#'
#' Creates a data frame indicating primer starts and ends.
#'
#' @param end.position End positions of primers.
#' @param p.lens Desired primer lengths.
#' @param start.posiion The start positions of primers.
#' @param step.size A numeric giving the steps with which start
#' positions are cycled. Should be 1 for primer design (evaluate all positions)
#' and higher values can be used for windowing.
#' @param groups Character vector with group annotation.
#' @return Data frame with ranges for initial primers.
#' @keywords internal
create.primer.ranges <- function(end.position, p.lens, start.position,
                                 step.size = 1, groups = NULL) {
    primer.ranges <- vector("list", length(p.lens))
    if (length(groups) == 0) {
        groups <- rep("default", length(end.position))
    }
    for (i in seq_along(p.lens)) {
        # ensure that we only select regions that aren't larger than the allowed binding region:
        p.len <- min(p.lens[i], end.position[i] - start.position[i])
        starts <- seq(start.position[i], end.position[i] - p.len + 1, step.size)
        ends <- starts + p.len - 1
        idx <- starts >= 1 & ends <= end.position[i]
        primer.ranges[[i]] <- data.frame(Start = starts[idx], End = ends[idx],
                                         Group = rep(groups[i], length(starts[idx])))
    }
    primer.range <- do.call(rbind, primer.ranges)
    primer.range <- plyr::ddply(primer.range, "Group", 
                        plyr::summarize,
                        Start = unique(substitute(Start)),
                        End = unique(substitute(End)))
    return(primer.range)
}
#' Conservation Scores
#'
#' Scores the conservation of alignment regions.
#'
#' @param primer.range A data frame with starts/ends of primers.
#' @param ali.entropy Entropies corresponding to the alignment
#' @return Entropies indicating conservation (similarity) of regions.
#' @keywords internal
score.conservation <- function(primer.range, ali.entropy) {
    # mean entropy per region
    scores <- unlist(lapply(seq_len(nrow(primer.range)), function(x) mean(ali.entropy[primer.range$Start[x]:primer.range$End[x]])))
    return(scores)
}
#' Selection by Conservation
#'
#' Selects primer regions for initialization of primers according
#' to their conservation scores.
#'
#' The conservation scores are computed using the entropies computed from
#' the alignment of the template sequence regions of interests.
#'
#' @param primer.range Data frame with primer starts/stops.
#' @param ali.entropy Entropy values for the alignment.
#' @param conservation Desired ratio of primer conservation.
#' Only regions with a conservation of at least \code{conservation}
#' are considered for the initialization of primers.
#' @param bin \code{DNABin} alignment of templates.
#' @return Updated primer regions according to the desired conservation.
#' @keywords internal
select.primer.region.by.conservation <- function(primer.range, ali.entropy, 
        conservation, bin) {
    # conservation: select only the top percent (ratio) of alignments 100% ->
    # consider all sub-alignments
    if (conservation == 1) {
        # conservation is 1 -> return top 100% of conserved ranges -> all candidate
        # ranges
        return(primer.range)
    }
    # otherwise: limit to conserved regions yielding mostly primers with required
    # length
    primer.candidates <- lapply(seq_len(nrow(primer.range)), function(x) bin[, primer.range[x, 
        "Start"]:primer.range[x, "End"]])  # sub-alignment
    # don't consider regions that are basically only gaps
    gap.char <- "-" # TODO as input
    # determine the gap positions in every possible binding region
    gap.idx <- detect.gap.columns(primer.candidates, gap.cutoff = 0.95, gap.char = gap.char) # TODO: improve the speed of gap detection
    # select only regions without major gap columns
    allowed.gaps <- 0
    sel.idx <- NULL
    selected <- FALSE
    while (!selected) { # i think this region selection doesn't work for 'rev'
        sel.idx <- which(unlist(lapply(gap.idx, function(x) length(x) <= allowed.gaps)))
        if (length(sel.idx) != 0) {
            # region found
            selected <- TRUE
        } else {
            # allow more gaps in selected regions
            allowed.gaps <- allowed.gaps + 1
        }
    }
    # turn selected index into idx for the whole alignment
    selected.range <- primer.range[sel.idx,]
    #warning("Need to check that the right regions are selected here!")
    c.score <- score.conservation(selected.range, ali.entropy)
    q <- quantile(c.score, seq(0, 1, 0.01))
    # entropy cutoff for the specified conservation level:
    sel.score <- q[paste(round(conservation, 2) * 100, "%", sep = "")]
    #print(sel.score)
    # index of regions below the entropy cutoff
    idx <- which(c.score <= sel.score)
    new.primer.range <- primer.range[sel.idx[idx], ]
    # add scores:
    new.primer.range$Entropy <- c.score[idx]
    return(new.primer.range)
}
#' Hierarchical Clustering.
#'
#' Performs hierarchical clustering on aligned primer sequences.
#'
#' The clustering is performed to identify similar groups of primer candidates
#' that can be merged to form degenerate primers.
#'
#' @param primer.candidates Alignment of primer candidates.
#' @return Phylogeny of the input \code{primer.candidates}.
#' @keywords internal
hclust.tree <- function(primer.candidates) {
    if (nrow(primer.candidates) < 2) {
        # cannot build a tree for a small number of seqs
        return(NULL)
    }
    # dist: use Hamming distance (nbr of required substitutions)
    s <- as.character(primer.candidates)
    seqs <- unlist(lapply(seq_len(nrow(s)), function(x) paste(s[x,], collapse = "")))
    names(seqs) <- rownames(primer.candidates)
    # don't parallelize here -> nthread =1, otherwise foreach loop hangs 
    dist <- stringdist::stringdistmatrix(seqs, method = "hamming",
                                        useNames = "names", nthread = 1) 
    #dist <- cluster::daisy(data.frame(as.character(primer.candidates)), metric = "gower")
    clusters <- stats::hclust(dist)
    tree <- ape::as.phylo(clusters)
    return(tree)
}
#' Tree-based Initialization of Primers.
#'
#' Creates a set of candidate primers using a tree-based algorithm.
#'
#' First, primers are aligned and their sequence similarity is determined to
#' compute a phylogenetic tree using hierarchical clustering.
#' Next, the tree is traversed from leaves to top in order to identify
#' groups of primers that can be merged (consensus) without exceeding the maximum
#' degeneracy of primers. 
#'
#' @param seqs Template sequences.
#' @param seq.IDs Identifiers of template sequences.
#' @param seq.groups Group identifiers of template sequences.
#' @param start For each template the start of the interval where primers 
#' should be created.
#' @param end For each template the end of the interval where primers
#' should be created.
#' @param primer.lengths Vector of desired primer lengths.
#' @param allowed.region.definition Definition of allowed regions.
#' @param max.degen Maximal degeneracy of primers.
#' @param conservation Required conservation of template regions considered
#' for the generation of primers. Conservation identifies the top conserved 
#' percentile of possible primers. 
#' @param sample Sample name for the analysis.
#' @param identifier Identifier (e.g. for directionality).
#' @param updateProgress Shiny progress object.
#' @return A vector with initialized primers.
#' @keywords internal
create.primers.tree <- function(seqs, seq.IDs, seq.groups, start, end, primer.lengths, allowed.region.definition, 
    max.degen, conservation, sample = "", identifier = "", updateProgress = NULL) {
    # use template group as primer identifier or individual primers?
    use.group.identifier <- TRUE
    if (length(unique(seq.groups)) <= 1) {
        use.group.identifier <- FALSE
    } 
    ##### select region of interest
    primer.lengths <- unique(primer.lengths) # ensure that we don't create duplicate primers
    input.start <- start
    input.end <- end
    if (allowed.region.definition == "any") {
        # we're selecting the minimal extension (min primer length) here to ensure that
        # no shorter primers are generated that do not overlap with the target region.
        start <- sapply(start, function(x) max(x - min(primer.lengths) + 1, 1))
        end <- sapply(seq_along(end), function(x) min(end[x] + min(primer.lengths) - 
            1, nchar(seqs)[x]))
    }
    input.seqs <- seqs
    seqs <- substring(seqs, start, end)
    ######### align templates
	message("i) Aligning sequences")
    lex.ali <- align.seqs(seqs, seq.IDs)
    if (is.null(lex.ali)) {
        # alignment wasn't possible
        return(NULL)
    }
	# sanitize the sequences: (carriage returns cause problems in windows)
	lex.ali$seq <- gsub("\r", "", lex.ali$seq)
    bin <- ape::as.DNAbin(lex.ali)
    ali.entropy <- shannon.entropy(as.character(bin))  # determine conservation according to Shannon entropy
    ####### define primer extents
    primer.range <- create.primer.ranges(rep(ncol(bin), length(primer.lengths)), primer.lengths,
                                         rep(1, length(primer.lengths)))
    #print("primer range:")
    #print(primer.range)
    primer.range <- select.primer.region.by.conservation(primer.range, ali.entropy, 
                        conservation, bin)
    ####### create primers from subalignments for every possible primer placement
    i <- NULL
	message("ii) Hierarchical clustering and tree construction")
    min.len <- min(primer.lengths)
    tree.data <- vector("list", nrow(primer.range))
    for (i in seq_len(nrow(primer.range))) { 
        # n.b. removed foreach loop due to memory issues -> no big impact on runtime here
    #tree.seqs <- foreach(i = seq_len(nrow(primer.range)), .combine = "rbind") %dopar% {
        if (is.function(updateProgress)) {
            detail <- identifier
            updateProgress(1/nrow(primer.range), detail, "inc")
        }
        range <- primer.range[i, "Start"]:primer.range[i, "End"]
        primer.candidates <- bin[, range]  # sub-alignment
        t.seqs <- NULL
        # n.b.: disadvantage -> for gappy regions, no primers may be generated here since primers would be too short!
        primer.candidates <- filter.primer.candidates(primer.candidates, min.len)
        ####### construct primers from tree
        tree <- hclust.tree(primer.candidates)
        t.seqs <- get.tree.seqs(tree, max.degen, primer.candidates)
        if (length(t.seqs) != 0 && nrow(t.seqs) != 0) {
            t.seqs <- t.seqs[nchar(t.seqs$Sequence) >= min.len, ]  # select only seqs of target length (can be shorter after merge & gap removal)
            node.ids <- strsplit(t.seqs$ID, split = ",")
            if (use.group.identifier) {
                # use group identifiers 
                m <- lapply(node.ids, function(x) match(x, seq.IDs))
                seq.ids <- sapply(m, function(x) paste0(unique(seq.groups[x]), collapse = ","))
            } else {
                # use the first sequence identifier
                seq.ids <- unlist(lapply(node.ids, function(x) x[1]))
            }
            primer.identifiers <- get.primer.identifier.string(sample, 
                    paste("D_", seq.ids, "_", t.seqs$Degeneracy, sep = ""), 
                    i, primer.range[i, "Start"], primer.range[i, "End"], 
                    identifier, t.seqs$Sequence)
            t.seqs$ID <- primer.identifiers
        }
        tree.data[[i]] <- t.seqs
    }
    tree.seqs <- do.call(rbind, tree.data)
    if (length(tree.seqs) == 0) {
        warning("No degenerate primers were created (gappy alignments, degeneration cutoff?)!")
    }
	message("iii) Augmentation with naive primers")
    # A) create naive primers in the target regions using the alignment to respect conservation:
    ali.matrix <- as.character(bin)
    rownames(ali.matrix) <- NULL
    #print("primer range:")
    #print(primer.range)
    naive.primers <- unlist(lapply(seq_len(nrow(primer.range)), function(x) {
                    range <- primer.range[x, "Start"]:primer.range[x, "End"]
                    if (use.group.identifier) {
                        # use group identifiers 
                        use.ids <- seq.groups
                    } else {
                        # use the first sequence identifier
                        use.ids <- seq.IDs
                    }
                    id <- get.primer.identifier.string(sample, paste0("N_", abbreviate(use.ids, 5), "_", seq_along(use.ids)), 
                                 x, min(range), max(range), identifier, 
                                 paste(rep("X", max(range) - min(range) + 1), collapse = ""))
                    out <- apply(ali.matrix, 1, function(s)  {
                        s <- s[range] # select substring
                        non.gap.idx <- which(s != "-") # remove gaps
                        paste(s[non.gap.idx], collapse = "")
                    })
                    # filter for length
                    len.ok <- which(nchar(out) >= min(primer.lengths))
                    names(out) <- id
                    out <- out[len.ok]
                }))
    # B) create primers for templates for which no primers have been generated yet (due to gappy alignments)
    if (length(naive.primers) != 0) { 
        naive.IDs <- sapply(strsplit(names(naive.primers), split = "\\|"), function(x) x[[1]])
    } else {
        naive.IDs <- NULL
        warning("Naive primer generation from alignment failed. Creating primers from sequences instead.")
    }
    missing.idx <- match(setdiff(seq.IDs, naive.IDs), seq.IDs)
    if (length(missing.idx) != 0) {
        # create primers that are still missing
        missing.primers <- create.primers.naive(input.seqs[missing.idx], seq.groups[missing.idx], 
            seq.IDs[missing.idx], input.start[missing.idx], input.end[missing.idx], primer.lengths, 
            allowed.region.definition, max.degen, sample, identifier, updateProgress)
        naive.primers <- c(naive.primers, missing.primers)
    }
    tree.primers <- tree.seqs$Sequence
    names(tree.primers) <- tree.seqs$ID
    all.seqs <- c(naive.primers, tree.primers)
    # remove duplicates
    all.seqs <- all.seqs[!duplicated(all.seqs)]
    return(all.seqs)
}
#' Shannon Entropy
#' 
#' Computation of Shannon entropy for an alignment.
#'
#' @param ali An alignment of primer sequences.
#' @return The Shannon entropy for the alignment.
#' @keywords internal
shannon.entropy <- function(ali) {
    #gapped.codemap <- Biostrings::IUPAC_CODE_MAP
    #gapped.codemap["-"] <- "-"
    # don't consider ambiguous positions to keep the score tight.
    symbols <- c("a", "c", "g", "t")
    all.counts <- apply(ali, 2, function(x) table(factor(x, levels = symbols)))
    freqs <- apply(all.counts, 2, function(x) x/nrow(ali))
    entropy <- apply(freqs, 2, function(x) -sum(x[x!=0] * log(x[x != 0], base = length(symbols))))
    return(entropy)
}
