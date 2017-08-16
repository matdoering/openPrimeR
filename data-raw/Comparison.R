####################
# Create data for comparison of primers
#############
remove.all.duplicates.from.fasta <- function(primers, primers.fixed) {
    for (i in seq_along(primers)) {
        sample <- names(primers)[i]
        message(sample)
        primer.loc <- primers[i]
        print(primer.loc)
        primer.files <- list.files(path = primer.loc, pattern = ".fasta", full.names = TRUE)
        out.loc <- primers.fixed[sample]
        dir.create(out.loc, showWarnings=FALSE)
        for (j in seq_along(primer.files)) {
            primer.location <- primer.files[j]
            cur.primers <- read.fasta(primer.location)
            seqs <- sapply(cur.primers, function(x) paste(as.character(x), collapse = ""))
            # check for duplicates
            dup.idx <- which(duplicated(seqs))
            # remove > character at the start from annotation
            ids <- sapply(cur.primers, function(x) substring(attr(x, "Annot"), 2))
            if (length(dup.idx) != 0) {
                warning(paste("Duplicates found in: ", primer.location, sep = ""))
                seqs <- seqs[-dup.idx]
                ids <- ids[-dup.idx]
            }
            # remove reverse primers
            rev.idx <- which(grepl("_rev", ids))
            if (length(rev.idx) != 0) {
                message(paste("Removing reverse primers: ", paste(ids[rev.idx], collapse = ","), sep = ""))
                seqs <- seqs[-rev.idx]
                ids <- ids[-rev.idx]
            }
            out.file <- file.path(out.loc, basename(primer.location))
            message(out.file)
            write.fasta(as.list(seqs), names = as.list(ids), file.out = out.file)
        }
    }
}
create.primer.set.evaluation.data <- function(data, settings, comp.results.loc, 
                                              primer.locations, rm.partials, max.degen, fw.id, rev.id) {

    primer.folders <- list.files(path = primer.locations, full.names = TRUE)
    constraint.settings <- constraints(settings)
    na_salt_conc <- PCR(settings)$Na_concentration
    mg_salt_conc <- PCR(settings)$Mg_concentration
    k_salt_conc <- PCR(settings)$K_concentration
    tris_salt_conc <- PCR(settings)$Tris_concentration
    primer_conc <- PCR(settings)$primer_concentration
    for (i in seq_along(data)) {
        sample <- names(data)[i]
        message(sample)
        lex.df <- data[[i]]
        cur.results.loc <- file.path(comp.results.loc, sample)
        dir.create(cur.results.loc, showWarnings=FALSE)
        write.csv(lex.df, file = file.path(comp.results.loc, paste(sample, "_templates.csv", sep = "")))
        # list all primer sets in folder
        primer.files <- list.files(path = primer.folders[i], full.names = TRUE)
        for (j in seq_along(primer.files)) {
            primer.location <- primer.files[j]
            bname <- basename(primer.location)
            run.info <- unlist(strsplit(bname, split = "_"))
            if (grepl("IPS", bname)) {
                run <- run.info[2]
            } else {
                run <- run.info[1]
            }
            if (any(grepl("1st", run.info))) {
                run <- paste(run, "_1st", sep = "")
            } 
            if (any(grepl("2nd", run.info))) {
                run <- paste(run, "_2nd", sep = "")
            }
            run <- sub(".fasta", "", run)
            message(paste("Loading: ", primer.location, sep = ""))
            message(paste("Primer set: ", run, sep = ""))
            set.name <- sub(".fasta", "", basename(primer.location))
            primer.df  <- read_primers(primer.location, fw.id, rev.id, merge.ambig = "none", max.degen)
            mode.directionality <- get.analysis.mode(primer.df)
            # need primer cvg to compute annealing  temperature estimate
            active.constraints <- names(constraint.settings)
            constraint.df <- check_constraints(primer.df, lex.df, settings, 
                                active.constraints)
            # add Run info to constraint.df
            fname <- file.path(cur.results.loc, paste(sample, "_", set.name, ".csv", sep = ""))
            message(paste("Writing: ", fname, sep = ""))
            write.csv(constraint.df, file = fname, row.names =FALSE)
        }
    }
}

get_templates <- function() {
    # initial retrieval templates from IMGT
    script.loc <- system.file("inst", "shiny", "shiny_server", 
                              "extra_IMGT_template_set_extractor.py",
                    package = "openPrimeR")
    # define what we want to extract:
    species <- "Homo sapiens"
    loci <- c("IGH", "IGK", "IGL")
    funcs <- c("any", "functional")
    for (i in seq_along(loci)) {
        locus <- loci[i]
        for (j in seq_along(funcs)) {
            func <- funcs[i]
            # retrieve & store leader+variable in extdata folder
            retrieve.IMGT.templates(species, locus, func, TRUE)
        }
    }
    message("templates reside in extdata target dir now :-)")
}
get_primers <- function() {
    # initial retrieval of primers from IMGT
    # n.b.: primers are not written to extdata here because they were processed by hand
    script.loc <- system.file("data-raw", "IMGT_primer_set_extractor.py",
                    package = "openPrimeR")
    cmd <- paste("Python", script.loc)
    res <- system(cmd) # stored in data/IMGT_data/IMGT_primers
    # combine literature and IMGT primer sets
    data.dir <- file.path("data", "IMGT_data")
    IMGT.loc <- file.path(data.dir, "IMGT_primers")
    # n.b.: for literature primers: apply mac2unix first to change formatting from mac
    lit.loc <- file.path(data.dir, "literature_primers")
    primer.loc <- file.path(data.dir, "primers_raw")
    # copy primers from both locations
    IMGT.dirs <- list.files(IMGT.loc, full.names = TRUE)
    out.dirs <- file.path(primer.loc, basename(IMGT.dirs))
    names(out.dirs) <- basename(IMGT.dirs)
    lit.dirs <- list.files(lit.loc, full.names = TRUE)
    for (i in seq_along(out.dirs)) {
        out.dir <- out.dirs[i]
        dir.create(out.dir, recursive = TRUE)
        res <- dir.copy(IMGT.dirs[i], out.dir, overwrite = TRUE)
        res <- dir.copy(lit.dirs[i], out.dir, overwrite = TRUE)
    }
    warning(paste("extdata was not updated with primers, need to manually modify adaptors.",
                "If you don't want this, use the existing processed primers."))
    # adaptor removal was done by hand:
    ## Modify primers (remove duplicates, remove adapters: restriction sites, overhangs, etc)
    #fixed.primers.loc <- file.path(data.dir, "primers")
    #dir.create(fixed.primers.loc, showWarnings=FALSE)
    #primers.fixed <- file.path(fixed.primers.loc, names(PRIMERS))
    #names(primers.fixed) <- names(PRIMERS)
    #raw.primers <- file.path(out.dirs)
    #remove.all.duplicates.from.fasta(out.dirs, primers.fixed) # only needs to be called once to write the primers without duplicates
    ## TODO: fix adapters when we know the stats of each set (low cvg primers -> probably have adaptors .. )
    ## now, copy primers to extdata/ folder of the openPrimeR package
    #target.dir <- system.file("extdata", "IMGT_data", "primers")
    #res <- dir.copy(primer.loc, target.dir)
    #message("primers reside in extdata target dir now :-)")
}
generate_evaluation_primers <- function(settings) {
    # evaluates the comparison primer sets according  to the provided settings
    # load templates:
    template.loc <- system.file("extdata", "IMGT_data", "templates",
                                package = "openPrimeR")
    hdr.structure <- c("ACCESSION", "GROUP", "SPECIES", "FUNCTION")

    template.locs <- file.path(template.loc, c("Homo_sapiens_IGH_functional_exon_exp.fasta",
                                               "Homo_sapiens_IGK_functional_exon.fasta",
                                               "Homo_sapiens_IGL_functional_exon.fasta"))
    leader.locs <- file.path(template.loc, c("Homo_sapiens_IGH_functional_leader_exp.fasta",
                                               "Homo_sapiens_IGK_functional_leader.fasta",
                                               "Homo_sapiens_IGL_functional_leader.fasta"))
    functional.data <- lapply(seq_along(template.locs), function(x) {
                            template.file <- template.locs[x]
                            leader.file <- leader.locs[x]
                            template.df <- read_templates(template.file, hdr.structure, "|", "GROUP", 
                                                          rm.keywords = "partial")
                            template.df <- assign_binding_regions(template.df, leader.file, NULL)
                            template.df <- adjust_binding_regions(template.df, 
                                           c(-max(template.df$Allowed_End_fw), 0), NULL) # allow binding in the first exon base as well
    })
    names(functional.data) <- c("IGH", "IGK", "IGL") # change names for shiny app comparison data
    # store template data in comparison folder:
    for (i in seq_along(functional.data)) {
        fname <- file.path(system.file("extdata", "IMGT_data", "comparison", "templates", package = "openPrimeR"),
                    paste0(names(functional.data)[i], "_templates.csv"))
        write_templates(functional.data[[i]], fname, "CSV")
    }
    # primers to load for comparison:
    primer.locations <- system.file("extdata", "IMGT_data", "primers",
                                    package = "openPrimeR")
    # output of set evaluations:
    comp.primer.out <- system.file("extdata", "IMGT_data", "comparison",
                                    "primer_sets", package = "openPrimeR")
    create.primer.set.evaluation.data(functional.data, settings,
                                      comp.primer.out, primer.locations, 
                                      rm.partials = TRUE, max.degen = max.degen,
                                      fw.id = fw.id, rev.id = rev.id)
}
generate_Comparison <- function() {
    # store comparison data for IGH: primers and templates  
    comp.files.primer <- list.files(path = system.file("extdata", "IMGT_data", 
                            "comparison", "primer_sets", "IGH", package = "openPrimeR"),
                            pattern = "*\\.csv", full.names = TRUE)
    comp.files.seqs <- rep(system.file("extdata", "IMGT_data", "comparison", 
                        "templates", "IGH_templates.csv", package = "openPrimeR"),
                        length(comp.files.primer))
    primer.data <- read_primers(comp.files.primer)
    template.data <- read_templates(comp.files.seqs)
    # update template cvg:
    template.data <- lapply(seq_along(template.data), function(x) {
                            update_template_cvg(template.data[[x]], primer.data[[x]])
    })
    out.loc <- system.file("data", "Comparison.rda", package = "openPrimeR")
    save(primer.data, template.data, file = out.loc, compress = "xz")
}
###########
# SOURCES
#############
devtools::load_all("src/openPrimeR")
source("src/analyses/load_data.R")
source("src/analyses/set_default_parameters.R") # load settings
source(system.file("inst", "shiny", "shiny_server", "extra_IO_shiny.R", package = "openPrimeR"))

refresh_src <- FALSE # should source primers/templates be retrieved once again?
if (refresh_src) {
    # WARNING: should never be run as primers were selected by hand afterwards ..
    get_templates()
    get_primers()
}
generate_evaluation_primers(settings) # evaluate primers
# n.b. need to copy templates from primer to template folder manually atm
generate_Comparison() # create rda file
