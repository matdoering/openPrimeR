################
# Generate procedures for creating internal Rdata here.
# o Computes the reference distribution for Fisher's exact test on the number of
# fulfilled constraints
#################


#' Creation of Reference Data for Significance Tests
#'
#' Creates reference data for significance test on the ratio of
#' fulfilled constraints.
#'
#' @param set.name Primer set to exclude from the data.
#' @param top.k Select only the top k coverage sets for every locus. 
#' If \code{top.k} is \code{NULL}, all primer sets 
#' are selected to form the reference distribution.
#' @param settings Settings for evaluating constraints that are not part of the original data.
#' @keywords internal
sig_eval_ref_data <- function(set.name = "openPrimeR2017", top.k = NULL, settings) {
    # load comparison data as reference frequencies
    # other options: kolmogorov smirnov (compare 2 distributions)
    # F-test to compare two variances (mine should be more narrow -> close to the target)
    comp.loc <- system.file("extdata", "IMGT_data", "comparison", "primer_sets",
                    package = "openPrimeR")
    comp.folders <- list.dirs(comp.loc, full.names = TRUE, recursive = FALSE)
    comp.loc.t <- list.files(
                    system.file("extdata", "IMGT_data", "comparison", 
                                "templates",package = "openPrimeR"),
                    full.names = TRUE)
    count.data <- NULL
    template.data <- read_templates(comp.loc.t)
    for (i in seq_along(comp.folders)) { # for each Ig group (IGHV etc)
        primer.set.location <- comp.folders[i]
        primer.files <- list.files(path = primer.set.location, pattern = ".csv", full.names = TRUE)
        primer.data <- read_primers(primer.files)
        # add melting_temp_diff column here
        # TODO: use settings here?
        primer.data <- lapply(primer.data, function(x) cbind(x, "EVAL_melting_temp_diff" = x$melting_temp_diff < constraints(settings)$melting_temp_diff))
        template.df <- template.data[[i]]
        locus <- template.df$Run[1]
        cvg <- sapply(primer.data, function(x) get_cvg_ratio(x, template.df))
        # select reference set for this locus:
        if (!is.null(top.k)) {
            # selet top-k coverage sets
            o <- order(cvg, decreasing = TRUE)
            ref.sets <- primer.data[o[1:top.k]]
            if (set.name %in% names(ref.sets)) {
                # need another reference to have k sets (will remove `set.name` later
                ref.sets <- primer.data[o[1:(top.k+1)]]
            } 
        } else {
            # select all cvg sets
            ref.sets <- primer.data
        }

        # count number of fulfilled constraints of all other sets
        # construct count data
        ref.df <- asS3(do.call(my_rbind, ref.sets))
        count.df <- create_fulfilled_counts(ref.df)
        count.df <- cbind(Locus = locus, count.df)
        count.data <- my_rbind(count.data, count.df)
    }       
    # exclude my set from ref.sets
    ref.data <- count.data[count.data$Run != set.name,]
    # summarize across all primer sets
    ref.data <- plyr::numcolwise(function(x) sum(x, na.rm = TRUE))(ref.data)
    return(ref.data)
}
generate.country.dists <- function() {
    # compute distances between countries for determining the fastest mirror:
    # computed once and stored library(cshapes) library(countrycode)
    library(cshapes)
    library(countrycode)
    dmat <- distmatrix(as.Date("2002-1-1"), type = "capdist", useGW = FALSE)  # correlates of war
    my.names <- tolower(countrycode(rownames(dmat), "cown", "iso2c"))  # turn countries to iso with 2 character encoding as used by R
    country.dists <- dmat
    colnames(country.dists) <- my.names
    rownames(country.dists) <- my.names
    return(country.dists)
}
#' Explore the properties of different coverage distribtuions.
#' @param template.df List with templates.
#' @param primer.data List with primers.
#' @param k Length of primers.
#' @return Does not return anything but performs exploration.
#' @keywords internal
explore.dists <- function(template.data, k = 18) {
    # library(logspline)
    for (i in seq_along(template.data)) {
        template.df <- template.data[[i]]
        sample <- names(template.data)[i]
        print(sample)
        cvg.df <- estimate.cvg(template.df, k = 18, "both", sample = sample)
        x <- c(cvg.df$fw$Coverage_Ratio, cvg.df$rev$Coverage_Ratio)
        png(paste0(sample, "_fit_discovery.png"))
        descdist(x, discrete = FALSE)
        dev.off()
        fit.beta <- fitdist(x, "beta") # beta has high error
        print(fit.beta)
        fit.gamma <- fitdist(x, "gamma")
        print(fit.gamma)
        png(paste0(sample, "_gamma_fit.png"))
        plot(fit.gamma)
        dev.off()
        png(paste0(sample, "_beta_fit.png"))
        plot(fit.beta)
        dev.off()
    }
}
#' Creation of Reference Coverage Ratio Distributions.
#'
#' Creates reference coverage ratio distributions in order to classify
#' primer design problems into three categories: easy, medium, and hard.
#'
#' @return A list with three reference beta distributions of coverage ratios.
#' @keywords internal
create.ref.dists <- function() {
    #library(fitdistrplus)
    library(distrEx)
    # reference dists were determined empirically
    # using the commented code from explore.dists() to determine good fits
    ###################
    #data(Ippolito)
    #template.data <- list("HCV" = read_templates("data/HCV/target_seqs.fasta"), 
                          #"IGH" = template.df)
    #explore.dists(template.data, k = 18)
    #x.easy <- hist(rbeta(10000, shape1 = 0.9, shape2 = 15))
    #x.medium <- hist(rbeta(10000, shape1 = 0.60, shape2 = 40))
    # we made the hard distribution a bit like a normal distribution
    # since in this case the parameter estimation may not give the 
    # actual distribution of the data. this is only the case for hard problems
    # where the space of possible coverage ratios is tight.
    #x.hard.before <- hist(rbeta(10000, shape1 = 0.3, shape2 = 200))
    #x.hard <- hist(rbeta(10000, shape1 = 2, shape2 = 1000))
    #x.HCV <- rbeta(10000, shape1 = 9.96, shape2 = 1606)
    #x.IGH <- rbeta(10000, shape1 = 0.798, shape2 = 23.825)
    #hist(x.easy)
    #hist(x.medium)
    #hist(x.hard)
    #hist(x.HCV)
    #hist(x.IGH)
    #KL.matrix <- cbind(Easy = ref.dist$easy, Medium = ref.dist$medium, Hard = ref.dist$hard, TestVeryEasy = Beta(0.8, 10), TestEasy = Beta(0.7, 30), TestMedium = Beta(1, 100), TestHard = Beta(0.3, 1000))
    #combis <- expand.grid(seq_len(ncol(KL.matrix)), seq_len(ncol(KL.matrix)))
    #dist <- matrix(rep(NA, nrow(combis)), nrow = ncol(KL.matrix))
    #colnames(dist) <- colnames(KL.matrix)
    #rownames(dist) <- colnames(KL.matrix)
    #for (i in 1:nrow(combis)) {
        #x <- combis[i,1]
        #y <- combis[i,2]
        #dist[x,y] <-  TotalVarDist(KL.matrix[1,x][[1]], unlist(KL.matrix[1,y][[1]])) 
    #}
    # store three reference distributions
    ref.dist <- list("very_easy" = Beta(1, 10),
                     "easy" = Beta(0.8, 20), 
                    "medium" = Beta(0.60, 40), 
                    "hard" = Beta(0.3, 200),
                    "very_hard" = Beta(10, 1000))
    return(ref.dist)
}
get_model_formula <- function() {
    # formula used for creating coverage model
    # use only 3' hexamer mutation position
    return(Experimental_Coverage ~ Position_3terminusLocal + annealing_DeltaG)
}
get_model_formula_log <- function() {
    # formula used for creating coverage model
    # log2 transformation: ensure that long primers are not overrated!
    return(Experimental_Coverage ~ log2(Position_3terminus) + annealing_DeltaG)
}
get.train.idx <- function(feature.matrix) {
    # sample from training set to keep the data more balanced: terminal mismatch pos & low DeltaG are overrepresented! -> sample by deltaG for stratification
    deltaG.cuts <- -c(0, 5, 9, 13, Inf) 
    strat <- cut(feature.matrix$annealing_DeltaG, deltaG.cuts)
    strat.types <- unique(strat)
    strat.counts <- sapply(strat.types, function(x) length(which(x == strat)))
    set.seed(1)
    n.samples <- min(strat.counts)
    # controlling for equal distribution of "Run": from 300 (equal distribution of DeltaG) to 168 samples
    train.idx <- vector("list", length(strat.types))
    set.seed(1)
    for (i in seq_along(strat.types)) {
        s <- strat.types[i]
        idx <- which(strat == s)
        cur.runs <- feature.matrix$Run[idx]
        # select from each run randomly:
        run.dist <- table(cur.runs)
        run.sample.size <- min(run.dist, n.samples/2)
        run.ids <- names(run.dist)
        sel <- unlist(lapply(run.ids, function(x) sample(idx[which(feature.matrix$Run[idx] == x)], run.sample.size)))
        #sel <- sapply(sample(idx, n.samples)
        train.idx[[i]] <- sel
    }
    train.idx <- unlist(train.idx)
    return(train.idx)
}
create.cvg.model <- function(feature.matrix) {
    # create model for predicting coverage
    train.idx <- get.train.idx(feature.matrix) # remove bias from training data
    cvg.model <- glm(get_model_formula(), family = "binomial", data = feature.matrix[train.idx,])
    return(cvg.model)
}
create.FPR.table <- function(feature.matrix) {
    # perform 10 repeats of 5-fold cross-validation using logistic regression models
    library(caret)
    library(e1071)
    train.idx <- get.train.idx(feature.matrix) # remove bias from training data
    ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 10, savePredictions = TRUE, classProbs = TRUE)
    model <- train(get_model_formula(), data = feature.matrix, method = "glm", 
                family="binomial", trControl = ctrl, tuneLength = 5, subset = train.idx)
    # retrieve all fold-predictions (coverage probabilities) across all repetitions
    cutoffs <- seq(0, 1, 0.0001)
    n.pos <- length(which(model$pred$obs == "Covered"))
    n.neg <- length(which(model$pred$obs == "Uncovered"))
    FPR.data <- lapply(cutoffs, function(x) {
        print(x)
        cut.pred <- factor(ifelse(model$pred$Covered > x, "Covered", "Uncovered"), levels = c("Uncovered", "Covered"))
        tpr <- length(which(model$pred$obs == "Covered" & cut.pred == "Covered")) / n.pos
        fpr <- length(which(model$pred$obs == "Uncovered" & cut.pred == "Covered")) / n.neg
        data.frame(Cutoff = x, FPR = fpr, TPR = tpr)
    })
    FPR.df <- do.call(rbind, FPR.data)
    # verify:
    # FPR.df[500,]
    return(FPR.df)
}
#########
# generate p-value ref data
#######
# use default settings for evaluating melting temp diff:
filename <- system.file("extdata", "settings", 
                  "C_Taq_PCR_high_stringency.xml", package = "openPrimeR")
settings <- read_settings(filename)
ref.data <- sig_eval_ref_data(settings = settings)
REF.pass.counts <- ref.data[, !grepl("failure", colnames(ref.data))]
REF.fail.counts <- ref.data[, grepl("failure", colnames(ref.data))]
colnames(REF.fail.counts) <- gsub("_failure", "", colnames(REF.fail.counts))
###
# generate common restriction site Rdata
####
# seqRFLP data from REBASE:
# pkg: seqRFLP (>= 1.0.1)
utils::data("enzdata", package = "seqRFLP")
# common restriction enzymes retrieved from:
# https://www.addgene.org/mol-bio-reference/restriction-enzymes/
# annotate common sites in enzdata
common.file <- file.path("src", "openPrimeR", "data-raw", "common_restriction_sites.txt")
common.sites <- as.character(read.table(common.file)[,1])
m <- match(common.sites, enzdata$nam)
enzdata$common <- FALSE
enzdata$common[m[!is.na(m)]] <- TRUE
# re-format search strings
# apostrophe -> indicates the sense cutting position, 
    # e.g. G_ACGT'C indicates a cut into ...GACGT and C...
# underscore -> indicates the antisense cutting position, 
    # e.g. G_ACGT'C indicates a cut into ...G and ACGTC....
enzdata$rep <- toupper(
            gsub("_", "", 
                gsub("'", "", enzdata$site, fixed=TRUE), 
                fixed = TRUE))
####
# design problem classification -> coverage distribution references
####
cvg.ref.dists <- create.ref.dists()
#######
# build supervised model for predicting primer coverage:
#######
load(RefCoverage) # load ref.data and feature.matrix
CVG_MODEL <- create.cvg.model(feature.matrix)
FPR_TABLE <- create.FPR.table(feature.matrix)
########
# create country dists for install helper function for the shiny app
########
country.dists <- generate.country.dists() # selection of CRAN mirrors using closest distance to countries
# store country dists to extdata as shiny app isn't exactly in the package (cannot access sysdata)
save(country.dists, file = file.path(system.file("extdata", package = "openPrimeR"), "country_dists.Rdata"))
############
# store sysdata
##########
devtools::use_data(REF.pass.counts,  # p-value: pass counts for constraint fulfillment
                   REF.fail.counts,  # p-value: fail counts for constraint fulfillment
                   cvg.ref.dists, # reference coverage distributions (problem difficulty estimation)
                   enzdata, # enzyme restriction sites
                   CVG_MODEL, # logistic model for predicting coveragee
                   FPR_TABLE, # conversion table from coverage-probabililties to FPR
                   internal = TRUE, 
                   pkg = "src/openPrimeR", overwrite = TRUE)
