#######
# Secondary structures (viennaRNA)
########

#' Read a Secondary Structure
#'
#' Reads the secondary structure output of ViennaRNA.
#'
#' @param fw.out Path to a ViennaRNA secondary structure output file.
#' @return Data frame with information on secondary structures.
#' @keywords internal
read.secondary.structure.raw <- function(fw.out) {
	if (!file.exists(fw.out)) {
		stop(paste("The input file:", fw.out,
					"does not exist."))
	}
    con <- file(fw.out)
	lines <- readLines(con)
    #lines <- try(suppressWarnings(readLines(con)), silent = TRUE)
    close(con)
    if (class(lines) == "try-error" || length(lines) == 0) {
        return(data.frame(Structure = character(0), deltaG = numeric(0), 
            stringsAsFactors = FALSE))
    }
    seq.idx <- seq(1, length(lines), 2)
    seqs <- lines[seq.idx]
    str.idx <- seq.idx + 1
    str <- lines[str.idx]
    str.s <- strsplit(str, split = "  \\(  ")
    expr <- "([\\.\\(\\)|x><]+)[[:space:]]\\([[[:space:]]]*(.[^\\)^[[:space:]]]*)"
    m <- stringr::str_match(str, expr)
    delta.G <- as.numeric(m[, 3])
    structure <- m[, 2]
    result <- data.frame(Structure = structure, deltaG = delta.G, stringsAsFactors = FALSE)
    # message(result)
    return(result)
}


#' Computation of Secondary Structures with ViennaRNA.
#'
#' Computes secondary structures using ViennaRNA.
#'
#' @param seqs The input sequences for which structures shall be computed.
#' @param annealing.temperature The temperature in degree Celsius 
#' at which to compute secondary structures.
#' @param constraints If provide
#' @param id An identifier for storing the files
#' @param folding.constraints Character vector specifying the folding
#' conditions for every input sequence.
#' For example the constraint \code{xxxxxx...} would forbid 
#' folding in the first 6 bases and allow folding in the last 3 bases.
#' @return A data frame with secondary structures.
#' @keywords internal
compute.structure.vienna <- function(seqs, annealing.temperature, 
                                    folding.constraints = NULL, id = "") {
    
    if (length(annealing.temperature) != 1) {
        stop("Can only call viennaRNA with a single annealing temperature!")
    }
    if (length(seqs) == 0) {
        return(NULL)
    }
	runID <- digest::digest(seqs, "md5")
    id <- paste0("viennaRNA_", id, "_", runID)
    #message("Computing secondary structures @ annealing temp: ", annealing.temperature)
    rna.fold <- Sys.which("RNAfold")
    if (rna.fold == "") {
        stop("ViennaRNA is not available on your system. Cannot compute secondary structures.")
    }
    # set parameters for the call to RNAfold:
    # we removed --noLP here, too slow now?
    param.string <- paste("--temp", annealing.temperature, "--noconv", "--noPS")
	old.dir <- NULL
    # out.file can be equal to inputfile as RNAfold appends '.fold' to the specified file
	if (grepl("windows", .Platform$OS.type)) {
		# viennaRNA can't handle windows filename ('sanitizes' them), 'filename-delim' arg doesnt help. 
		# -> use local files instead. need to be relative to the executable!!
		old.dir <- setwd(dirname(rna.fold))  # i have to do this such that i can avoid using absolute path specifiers
		#dir.create("temp", showWarnings=FALSE) # also doesnt work
		input.file <- paste0(id, ".txt")
		out.file <- input.file
	} else {
		input.file <- tempfile(pattern = id, fileext = ".txt")
		out.file <- input.file
	}
	# clean up result files to prevent problem with viennaRNA appending to existing file.
	on.exit({
		if (file.exists(out.file)) {
			file.remove(out.file)
		}
		# revert back to old working directory
		if (length(old.dir) != 0) {
			setwd(old.dir)
		}
	})
	file.string <- paste0("--infile=", input.file, " --outfile=", out.file)
    seqs <- toupper(seqs)
    if (length(folding.constraints) != 0) {
        # add constraints to input file for viennaRNA 
        seq.records <- sapply(seq_along(seqs), function(x) paste(seqs[x], 
            folding.constraints[x], sep = "\n"))
        writeLines(seq.records, input.file)
    } else {
        # temporary file containing the input sequences for viennaRNA
        writeLines(seqs, input.file)
    }
    if (!file.exists(input.file)) {
        stop(paste("IO error for secondary structure computation. ",
                "Could not write fw input file for viennaRNA at:",
                input.file))
    }
    if (length(folding.constraints) != 0) {
        # add constrained folding params:
        param.string <- paste(param.string, "-C", "--batch")
    }
    vienna.call <- paste(rna.fold, param.string, file.string)
	system2(rna.fold, args = c(param.string, file.string))
    # read the computed structures and energies from out.file:
    out.file <- paste0(out.file, ".fold")
	if (file.exists(out.file)) {
		result <- read.secondary.structure.raw(out.file)
	} else {
		result <- NULL
	}
    return(result)
}
#' Secondary Structure Computations.
#'
#' Computes the secondary structures of the input primers using ViennaRNA.
#' 
#' @param primer.df Primer data frame.
#' @param mode.directionality Direction of primers.
#' @param annealing.temperature Temperatures at which to compute secondary structures for every primer
#' @return Data frame with secondary structure information.
#' @keywords internal
#'
#' @references Lorenz, Ronny and Bernhart, Stephan H. and Höner zu Siederdissen, 
#' Christian and Tafer, Hakim and Flamm, Christoph and Stadler, Peter F. and Hofacker, Ivo L.
#' ViennaRNA Package 2.0
#' Algorithms for Molecular Biology, 6:1 26, 2011, doi:10.1186/1748-7188-6-26
compute.secondary.structures <- function(primer.df, mode.directionality = c("fw", "rev", "both"), annealing.temperature) {
	if (!is(primer.df, "Primers")) {
		stop("Please supply a valid set of primers.")
	}
	if (length(annealing.temperature) == 0) {
		stop("Please supply an annealing temperature.")
	}
    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' argument.")
    }
    mode.directionality <- match.arg(mode.directionality)
    fw.idx <- which(primer.df$Forward != "")
    rev.idx <- which(primer.df$Reverse != "")
    # create groups of annealing temperatures:
    # parallelize -> split into batches
    if (length(annealing.temperature) == 1)  {
        # split up seqs independent of annealing temp
        batches <- batchify(seq_len(nrow(primer.df)))
        annealing.temperature <- rep(annealing.temperature, length(batches))
    } else {
        # split up seqs dependent on annealing temp
        batches <- batchify(seq_len(nrow(primer.df)), annealing.temperature)
	    annealing.temperature <- as.numeric(names(batches))
    }
    # compute structures with ViennaRNA:
	i <- NULL
    str.fw <- foreach(i = seq_along(batches), .combine = "rbind") %dopar% {
        batch.idx <- batches[[i]]
        cur.temp <- annealing.temperature[i]
        str.fw <- compute.structure.vienna(primer.df$Forward[intersect(fw.idx, batch.idx)], cur.temp, id = "fw")
        #str.rev <- compute.structure.vienna(primer.df$Reverse[rev.idx], annealing.temperature, id = "rev")
    }
    str.rev <- foreach(i = seq_along(batches), .combine = "rbind") %dopar% {
        batch.idx <- batches[[i]]
        cur.temp <- annealing.temperature[i]
        str.rev <- compute.structure.vienna(primer.df$Reverse[intersect(rev.idx, batch.idx)], cur.temp, id = "rev")
    }

    # create order of indices for retrieval of results
    fw.idx.out <- intersect(unlist(batches), fw.idx)
    rev.idx.out <- intersect(unlist(batches), rev.idx)
    str.fw <- str.fw[match(fw.idx, fw.idx.out), ]
    str.rev <- str.rev[match(rev.idx, rev.idx.out), ]
    out <- data.frame(Structure_fw = rep(NA, nrow(primer.df)), Structure_rev = rep(NA, 
                nrow(primer.df)), Structure_deltaG_fw = rep(0, nrow(primer.df)), Structure_deltaG_rev = rep(0, 
                nrow(primer.df)), stringsAsFactors = FALSE)
    if (mode.directionality == "both") {
        # add info for missing primers
        res.fw <- str.fw[FALSE, ]
        res.rev <- str.rev[FALSE, ]
        res.fw[1:nrow(primer.df), ] <- NA
        res.rev[1:nrow(primer.df), ] <- NA
        res.fw[fw.idx, ] <- str.fw
        res.rev[rev.idx, ] <- str.rev
        out$Structure_fw <- res.fw$Structure
        out$Structure_deltaG_fw <- res.fw$deltaG
        out$Structure_rev <- res.rev$Structure
        out$Structure_deltaG_rev <- res.rev$deltaG
    } else if (mode.directionality == "fw") {
        out$Structure_fw[fw.idx] <- str.fw$Structure
        out$Structure_deltaG_fw[fw.idx] <- str.fw$deltaG
    } else {
        out$Structure_rev[rev.idx] <- str.rev$Structure
        out$Structure_deltaG_rev[rev.idx] <- str.rev$deltaG
    }
    if (length(out) != 0) {
        out$Structure_deltaG <- sapply(seq_len(nrow(out)), 
            function(x) min(out$Structure_deltaG_fw[x], 
                        out$Structure_deltaG_rev[x], na.rm = TRUE))
    }
    return(out)
}
