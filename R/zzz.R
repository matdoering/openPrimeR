####
# Actions to perform on loading/attaching the package
######
#' @import ggplot2 lpSolveAPI
#' @importFrom reshape2 melt dcast
#' @importFrom plyr ddply summarize arrange .
#' @importFrom foreach foreach %dopar%
#' @importFrom IRanges IRanges as.matrix
#' @importFrom Biostrings DNAStringSet
#' @importFrom RColorBrewer brewer.pal 
#' @importFrom grDevices colorRampPalette
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom BiocGenerics unlist start end
#' @importFrom magrittr %>%
#' @importFrom stats na.omit qnorm quantile sd 
#' @importFrom utils head read.csv read.delim setTxtProgressBar tail txtProgressBar write.csv write.table
#' @importFrom methods as callNextMethod new validObject setClass setMethod setGeneric setReplaceMethod is rbind2 cbind2 slot prototype
NULL # need to have some evaluated code here

#' Determination if Selenium is installed.
#'
#' Checks whether selenium module for python is installed on the system.
#'
#' @return \code{TRUE} is selenium for python is available,
#' \code{FALSE} otherwise.
#'
#' @keywords internal
selenium.installed <- function() {
	if (Sys.which("python") == "") {
		# python not available
		return(FALSE)
	}
    cmd <- "python -c 'import selenium'"
    ret <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (ret == 0) {
        # installed
        return(TRUE)
    } else {
        return(FALSE)
    }
}
#' Check Tool Installation
#'
#' Checks whether all required tools are installed.
#'
#' @param frontend Whether tool installation shall be checked for the frontend.
#' If \code{TRUE}, dependencies that are required only by the frontend are considered additionally.
#' @return \code{TRUE} for each installed tool, \code{FALSE} otherwise.
#' @keywords internal
check.tool.installation <- function(frontend = FALSE) {
    available.tools <- NULL
    # for melting temperatures
    available.tools["MELTING"] <- Sys.which("melting-batch") != ""
    # for secondary structures
    available.tools["ViennaRNA"] <- Sys.which("RNAfold") != ""
    # for DECIPHER
    available.tools["OligoArrayAux"] <- Sys.which("hybrid-min") != ""
     # for multiple sequence alignments
    available.tools["MAFFT"] <-  Sys.which("mafft") != "" 
    available.tools["Pandoc"] <- Sys.which("pandoc") != ""
    ## for IMGT data retrieval in frontend
    if (frontend) {
        available.tools["Selenium"] <- selenium.installed()
        available.tools["PhantomJS"] <- Sys.which("phantomjs") != ""
    }
    return(available.tools)
}
#' Check Functionality of Third-Party Tools.
#'
#' Checks whether all required tools should work.
#'
#' @param frontend Whether tool functionality shall be checked for the frontend.
#' @return \code{TRUE} for each functioning tool, \code{FALSE} for non-functioning tools.
#' @keywords internal
check.tool.function <- function(frontend = FALSE) {
    available.tools <- check.tool.installation(frontend)
    if (available.tools["OligoArrayAux"] && Sys.getenv("UNAFOLDDAT") == "") {
        available.tools["OligoArrayAux"] <- FALSE
        warning("Please define the UNAFOLDDAT variable in your path before using OligoArrayAux.")
    } else if (available.tools["OligoArrayAux"] && Sys.getenv("UNAFOLDDAT") != "") {
        out <- system("hybrid-min -n DNA -t 50 -T 50 -N 0.05 -E -q ACAGGTGCCCACTCCCAGGTGCAG CTGCACCTGGGAGTGGGCACCTGT", 
                    intern = FALSE, ignore.stdout = TRUE)
       if (out != 0) {
            # there was an error
            warning("oligoArrayAux failed checks: disabled. ")
            available.tools["OligoArrayAux"] <- FALSE
        }
    }
    # check whether custom MELTING parameter file exists
    melt.bin <- Sys.which("melting-batch")
    # custom tandem file removed: ...
    ##################
    #if (melt.bin != "") {
        #tandem.mm.file <- system.file("extdata", 
                        #"AllawiSantaluciaPeyret1997_1998_1999tanmm_mod.xml",
                        #package = "openPrimeR")
       #melt.config.file <- file.path(dirname(melt.bin), "..", "Data", 
                        #basename(tandem.mm.file))
       #if (!file.exists(melt.config.file)) {
            #warning("Could not find MELTING parameter file '", melt.config.file,  
                     #"', which should have been copied by openPrimeR from ",
                     #tandem.mm.file, ". Check your write permissions!")
            #available.tools["MELTING"] <- FALSE
        #}
    #}
    #################
    if (available.tools["Pandoc"] && Sys.which("pdflatex") == "") {
        # don't warn here, otherwise too many warnings are generated
        #warning("Cannot create reports with pandoc since LateX is missing.")
        available.tools["Pandoc"] <- FALSE
    }
    return(available.tools)
}
#' Copy MELTING Config File
#'
#' Copies modified MELTING tandem mismatch file to the MELTING data folder.
#'
#' @return TRUE if the file is available in the MELTING folder, FALSE otherwise.
#' @keywords internal
copy.melt.config <- function(melt.bin = NULL) {
    print("DEPRECATED")
	if (length(melt.bin) == 0) {
		melt.bin <- Sys.which("melting-batch")[1]
	}
    tandem.mm.file <- system.file("extdata", 
                        "AllawiSantaluciaPeyret1997_1998_1999tanmm_mod.xml",
                        package = "openPrimeR")
    if (tandem.mm.file == "") {
        warning("The MELTING config file is not present in the openPrimeR package.")
        return(FALSE)
    }
    if (melt.bin != "" ) {
        melt.config.file <- file.path(dirname(melt.bin),
                                "..", "Data", basename(tandem.mm.file))
        if (!file.exists(melt.config.file)) {
			message("Copying MELTING config to: ", melt.config.file)
            s <- file.copy(tandem.mm.file, melt.config.file)
            if (any(!s)) {
                warning("Could not copy MELTING config file to destination.")
            }
            return(all(s))
        } else {
            # file is available
            return(TRUE)
        }
    } else {
        return(FALSE)
    }
}
# actions to be performed when loading the package namespace
.onLoad <- function(libname, pkgname) {
    # copy melting config file if required
    #copy.melt.config() # not necessary anymore
    ################
    # Define package options
    #######
    # order in which constraints are computed (least runtime to highest)
    con.order <- c("primer_length", "gc_clamp", "gc_ratio", "no_runs", "no_repeats", 
               "melting_temp_range", "self_dimerization", "secondary_structure", 
               "primer_coverage", "primer_specificity", 
               "melting_temp_diff", "cross_dimerization")
    # order in which constraints are relaxed (start with least important constraint)
    relax.order <- c("primer_length", "primer_coverage", "no_repeats", "no_runs", "gc_clamp", "primer_specificity", "secondary_structure", "self_dimerization", "cross_dimerization", "gc_ratio", "melting_temp_range", "melting_temp_diff")
    available.constraints <- select.constraints(con.order) # select constraints that can be computed by installed software
    con.order <- con.order[con.order %in% available.constraints]
    relax.order <- relax.order[relax.order %in% available.constraints]
    plot.colors <- c("Constraint" = "Set1", "Group" = "Set2", 
                     "Run" = "Set3", "Primer" = "Accent")
    op <- options()
    op.openPrimeR <- list(
        openPrimeR.constraint_order = con.order,
        openPrimeR.relax_order = relax.order,
        openPrimeR.plot_colors = plot.colors,
        openPrimeR.plot_abbrev = 15 # limit label extent for plots
    )
    # only set options once:
    toset <- !(names(op.openPrimeR) %in% names(op))
    if (any(toset)) {
        options(op.openPrimeR[toset])
    }
    # do not provide any output, even if no var is assigned to this function: 
    invisible()
}
# actions to be performed when attaching the package
.onAttach <- function(libname, pkgname) {
    # add start up message
    available.tools <- check.tool.function()
    if (any(!available.tools)) {
        tool.df <- build.tool.overview(available.tools)
        out <- paste0("There are missing/non-functioning external tools.\n",
                "To use the full potential of openPrimeR, please make sure\n",
                "that all of the listed tools are installed:\n")
        idx <- which(!available.tools)
        tools <- paste0("o ", names(available.tools)[idx], 
              " (", tool.df$URL[match(names(available.tools[idx]), tool.df$Tool)], ")", 
              collapse = "\n")
        out <- paste(out, tools, sep = "")
        packageStartupMessage(out)
		# special warning for Pandoc (latex dependency)
		if (Sys.which("pdflatex") == "") {
			warning("'Pandoc' is non-functional, since 'pdflatex' is not installed on your system.")
		}
    }
    doParallel.available <- requireNamespace("doParallel", quietly = TRUE)
	is.win <- grepl("windows", .Platform$OS.type)
    if (is.win) {
        warning("Disabling parallel computations for Windows. Sorry.")
		return(NULL)
    }
    if (doParallel.available) { # no parallel support for windows at the moment (unserialize errors if we don't do clusterexport/clustercall)
        if (!foreach::getDoParRegistered()) {
			avail.cores <- parallel::detectCores()
            # always keep 1 core free to ensure that system doesn't freeze
            cores <- max(1, avail.cores - 1)
            doParallel::registerDoParallel(cores = cores)
            # also set mc.cores for 'mclapply'
            if (!"mc.cores" %in% names(options())) {
                options(list("mc.cores" = cores))
            }
			#message("REGISTERED PSOCK CLUSTER")
			# as soon as we register here, we get unserialize errors ...
        }
    } else {
        warning("Please install 'doParallel' to use multiple cores.")
    }
}
