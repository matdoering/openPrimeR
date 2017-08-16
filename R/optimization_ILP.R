#' Coverage Matrix
#'
#' Constructs a coverage matrix where rows indicate templates and columns indicate primers.
#'
#' Entry (i,j) in the matrix is equal to 1 if primer j covers template i and otherwise 0.
#'
#' @param primer.df Primer data frame.
#' @param template.df Template data frame.
#' @param constraints A character vector of coverage constraints
#' to be used as entries for the coverage matrix instead of the 0/1 encoding.
#' At its default setting (\code{NULL}), the 0/1 encoding is used.
#' @return The binary coverage matrix.
#' @keywords internal
get.coverage.matrix <- function(primer.df, template.df, constraints = NULL) {
    # create an mxn matrix where templates are rows and primers are columns a_ij is 1
    # if template i is covered by primer j, 0 else
    if (nrow(primer.df) == 0) {
        return(NULL)
    }
    template.coverage <- get.covered.templates(primer.df, template.df)
    if (any(!constraints  %in% c("primer_efficiency", "annealing_DeltaG"))) {
        stop("Incorrect coverage constraint: ", constraints)
    }
    if (length(constraints) > 1) {
        warning("Only one constraint possible at a time. Selecting the first one.")
        constraints <- constraints[1]
        print(constraints)
    }
    # transform to matrix
    if (length(constraints) == 0) {
        A <- matrix(rep(0, nrow(template.df) * nrow(primer.df)), nrow = nrow(template.df), ncol = nrow(primer.df))
    } else {
        # insert NA to identify templates that weren't covered
        A <- matrix(rep(NA, nrow(template.df) * nrow(primer.df)), nrow = nrow(template.df), ncol = nrow(primer.df))
        cvd <- covered.seqs.to.idx(primer.df$Covered_Seqs, template.df)
        con.vals <- lapply(strsplit(primer.df[, constraints], split = ","), function(vals) as.numeric(vals))
    }
    for (i in seq_along(template.coverage)) {
        # for every template
        primer.idx <- template.coverage[[i]]
        if (length(constraints) == 0) {
            A[i, primer.idx] <- 1
        } else {
            cur.cvd <- cvd[primer.idx]
            cur.idx <- sapply(cur.cvd, function(x) which(i == x))
            if (length(cur.idx) > 0) {
                # retrieve value @ cur.idx
                cur.values <- con.vals[primer.idx]
                sel.values <- sapply(seq_along(cur.values), function(x) round(as.numeric(cur.values[[x]][cur.idx[x]]), 3))
                A[i, primer.idx] <- sel.values
            }
        }
    }
    rownames(A) <- template.df$ID
    colnames(A) <- primer.df$ID
    return(A)
}
#' Addition of Dimerization Constraints
#'
#' Updates ILP formulation with dimerization constraints.
#'
#' @param lprec An ILP instance.
#' @param D.idx Data frame giving the indices of dimerizing primer pairs.
#' @param indices Row indices for setting the dimerization constraints in \code{lprec}.
#' @return \code{lprec} with added dimerization constraints.
#' @keywords internal
add.dimerization.constraints <- function(lprec, D.idx, indices) {
    # D.idx: pairs of primers that dimerize indices: rows in which to insert the
    # dimerization constraints in the lprec
    if (nrow(D.idx) > 0) {
        #pb <- txtProgressBar(min = 0, max = nrow(D.idx), style = 3)
        for (j in seq_len(nrow(D.idx))) {
            # pairs of dimerizing primers
            row <- D.idx[j, 1]
            col <- D.idx[j, 2]
            set.row(lprec, indices[j], rep(1, 2), indices = c(row, col))  # don't include pairs of dimerizing primers
            #setTxtProgressBar(pb, j)
        }
    }
    if (length(indices) != 0) {
        set.constr.type(lprec, rep("<=", length(indices)), indices)
        set.rhs(lprec, rep(1, length(indices)), indices)
        # don't set rownames -> takes quite long rownames <- paste('Dimer', indices, sep
        # = '') rownames(lprec)[indices] <- rownames
    }
    return(lprec)
}
#' Addition of Coverage Constraints.
#'
#' Adds coverage constraints to ILP instance.
#'
#' @param lprec An ILP instance.
#' @param covered.templates Indices of covered template sequences.
#' @param template.coverage List containing the indices of covering primers for each template.
#' @return \code{lprec} with coverage constraints.
#' @keywords internal
add.coverage.constraints <- function(lprec, covered.templates, template.coverage) {
    # add coverage constraints to ILP 'lprec'
    nbr.coverage.constraints <- length(covered.templates)
    indices <- 1:nbr.coverage.constraints
    cur.nbr.constraints <- 0
    if (nbr.coverage.constraints == 0) {
        indices <- NULL
    }
    #cat(paste("Adding ", length(covered.templates), " coverage constraints. \n", 
        #sep = ""))
    #if (length(covered.templates) > 0) {
        #pb <- txtProgressBar(min = 0, max = length(covered.templates), style = 3)
    #}
    for (j in seq_along(covered.templates)) {
        # select the idx of those primers that are covering the current template
        idx <- template.coverage[[covered.templates[j]]]
        set.row(lprec, j, rep(1, length(idx)), indices = idx)
        #setTxtProgressBar(pb, j)
    }
    if (length(indices) != 0) {
        set.constr.type(lprec, rep(">=", length(indices)), indices)
        set.rhs(lprec, rep(1, length(indices)), indices)
        # rownames <- paste('Cvg', indices, sep = '') rownames(lprec)[indices] <-
        # rownames
    }
    return(lprec)
}

#' Template Coverage.
#'
#' Determines the indices of covering primers for every template.
#'
#' @param cvg.matrix Binary matrix of covering events.
#' @return A list with covering primers for every template.
#' @keywords internal
J.cvg <- function(cvg.matrix) {
    if (length(cvg.matrix) == 0) {
        return(NULL)
    }
    idx <- lapply(seq_len(nrow(cvg.matrix)), function(x) which(cvg.matrix[x, ] > 
        0))
    return(idx)
}
#' Primer Coverage.
#'
#' Determines the indices of covered templates for every primer.
#'
#' @param cvg.matrix Binary matrix of covering events.
#' @return A list with covered templates for every primer.
#' @keywords internal
I.cvg <- function(cvg.matrix) {
    if (length(cvg.matrix) == 0) {
        return(NULL)
    }
    idx <- lapply(seq_len(ncol(cvg.matrix)), function(x) which(cvg.matrix[, x] > 
        0))
    return(idx)
}
#' Construct Coverage ILP.
#'
#' Constructs an ILP modeling the primer set cover problem.
#'
#' @param D Binary dimerization matrix.
#' @param cvg.matrix Binary coverage matrix.
#' @param time.limit Time limit for ILP optimization in seconds.
#' @param presolve.active Whether the ILP presolver should be used.
#' This is set to \code{FALSE} by default, since presolving may lead
#' to inferior solutions. However, for large problems presolving
#' might be useful.
#' @return An instance of the set cover ILP.
#' @keywords internal
ILPConstrained <- function(D, cvg.matrix, time.limit = NULL, presolve.active = FALSE) {
    # melting temperature and binding region conditions are relaxed. here, binding
    # region is modeled with another aux var.
    
    # D: dimerization with entry d_ij = 1 if primers i and j dimerize, 0 else
    # time.limit: stop model computation when timeout is surpassed cvg.matrix:
    # templates x primers. 1 if template is covered
    X <- ncol(cvg.matrix)  # nbr vars
    template.coverage <- J.cvg(cvg.matrix)
    covered.templates <- which(apply(cvg.matrix, 1, sum) != 0)  # number of templates we can cover
    nbr.coverage.constraints <- length(covered.templates)
    # remove symmetric entries in D -> halve the number of dimerization constraints
    D[lower.tri(D)] <- 0  # lower triangle doesn't contain additional info
    D.idx <- which(D == 1, arr.ind = TRUE)  # all pairs of dimerizing primers
    nbr.dimer.constraints <- nrow(D.idx)
    total.nbr.constraints <- nbr.coverage.constraints + nbr.dimer.constraints
    ####### create MILP with binary variables for each primer and real auxiliary variables
    ####### for the relaxed conditions
    lprec <- make.lp(total.nbr.constraints, X)  # 0 constraints and X vars
    name.lp(lprec, "openPrimeR_set_cover_ILP")
    for (col in 1:X) {
        # define primer vars as integer (bounds of 0,1)
        set.type(lprec, col, "binary")
    }
    coefs <- c(rep(1, X))  # costs: uniform
    set.objfn(lprec, coefs)
    # set ILP options
    if (presolve.active) {
        # timeout does not seem to work for loading ILP data, but only for simplex
        # optimization
        lp.info <- lp.control(lprec, sense = "min", timeout = ifelse(is.finite(time.limit), 
            time.limit, 0), verbose = "neutral", epslevel = "tight", presolve = c("lindep", 
            "rows", "probefix", "probereduce", "rowdominate", "impliedfree", "reducegcd", 
            "bounds", "duals", "sensduals", "mergerows", "cols", "colfixdual"))
    } else {
        # timeout does not seem to work for loading ILP data, but only for simplex
        # optimization
        lp.info <- lp.control(lprec, sense = "min", timeout = ifelse(is.finite(time.limit), 
            time.limit, 0), verbose = "neutral", epslevel = "tight")
    }
    rownames <- NULL
    ####### coverage constraint: every template that can be covered should have at least 1
    ####### covering primer
    cur.nbr.constraints <- 0
    lprec <- add.coverage.constraints(lprec, covered.templates, template.coverage)
    cur.nbr.constraints <- cur.nbr.constraints + length(covered.templates)
    # set new indices for lprec access of dimer constraint:
    if (nbr.dimer.constraints != 0) {
        indices <- (cur.nbr.constraints + 1):(cur.nbr.constraints + nbr.dimer.constraints)
    } else {
        indices <- NULL
    }
    ############### dimerization constraint
    lprec <- add.dimerization.constraints(lprec, D.idx, indices)
    cur.nbr.constraints <- cur.nbr.constraints + length(nbr.dimer.constraints)
    # write.lp(lprec, 'mynewLP.lp', 'lp') #write it to a file in LP format
    return(lprec)
}

#' Retrieval of ILP Decisions
#'
#' Retrieves ILP decision variables. 
#'
#' The original dimension of the ILP is required to determine the correct decisions
#' when presolve has been active and dimensions of the ILP might have changed.
#'
#' @param ILP A solved ILP instance.
#' @param original.dim Dimension of \code{ILP} before using presolve.
#' @return The ILP decision variables.
#' @keywords internal
get.ILP.vars <- function(ILP, original.dim = NULL) {
    if (length(original.dim) == 0) {
        return(get.variables(ILP))  # does not contain variables that were removed by presolve
    } else {
        # if presolve removed some variables, we need to know the original ILP
        # dimensionality, to find the 'original' decision variable solutions
        v <- get.primal.solution(ILP, orig = TRUE)  # contains also the variables that were removed by presolve
        vars <- v[(original.dim[1] + 1):length(v)]
        return(vars)
    }
}
#' Selection of Best ILP
#'
#' Selects the best solution from multiple solved ILP instances.
#'
#' @param ILP.df Data frame with ILP result properties.
#' @return The index of the best solution.
#' @keywords internal
select.best.ILP <- function(ILP.df) {
    # select best ILP according to the smallest objective function value after
    # selecting only those ILPs maximizing the coverage ILP.df: overview of solved LP
    # data
    if (length(ILP.df) == 0 || nrow(ILP.df) == 0) {
        warning("Cannot select best ILP as ILP.df was empty.")
        message(ILP.df)
        return(NULL)
    }
    if (all(is.na(ILP.df$Coverage))) {
        # arbitrarily select a set
        sel.idx <- 1
    } else {
        max.cvg <- max(ILP.df$Coverage, na.rm = TRUE)
        sel.idx <- which(ILP.df$Coverage == max.cvg)
        if (any(!is.na(ILP.df$MaxTempDiffToTargetTemp))) {
            sel.idx <- sel.idx[which.min(ILP.df$MaxTempDiffToTargetTemp[sel.idx])]
        } else {
            sel.idx <- sel.idx[1]
        }
    }
    return(sel.idx)
}
#' Construction of ILP Results.
#'
#' Constructs a data frame summarizing the properties of an ILP solution.
#'
#' @param ILP A solved ILP instance.
#' @param vars The ILP decision variables.
#' @param primer.df The primer data frame correspdong to the \code{ILP}.
#' @param template.df The template data frame.
#' @param i Index for the ILP.
#' @param target.temp Target melting temperature in Celsius.
#' @param time Runtime of the ILP.
#' @param deltaG_Cutoff Free energy cutoff used for the dimerization constraint.
#' @param deltaG_Limit The free energy boundary for dimerization.
#' @return Data frame summarizing the ILP solution.
#' @keywords internal
build.ILP.df <- function(ILP, vars, primer.df, template.df, i, target.temp, time = NA, 
    deltaG_Cutoff = NA, deltaG_Limit = NA) {
    # Builds a data frame summarizing a solved ILP
    
    # i: iteration time: time measurement of run need non-null entries for df...
    if (is.null(deltaG_Cutoff)) {
        deltaG_Cutoff <- NA
    }
    if (is.null(deltaG_Limit)) {
        deltaG_Limit <- NA 
    }
    data <- data.frame(Identifier = NA, Iteration = i, Set_Size = NA, Penalty = NA, 
        Coverage = NA, Iterations = NA, Nodes = NA, Time = time, MaxTempDiffToTargetTemp = NA, 
        Objective = NA, deltaG_Cutoff = deltaG_Cutoff, Variables = NA, 
        stringsAsFactors = FALSE)
    if (length(ILP) == 0 || nrow(primer.df) == 0) {
        # no solution found :-(
        return(data)
    }
    # vars <- get.ILP.vars(ILP, original.dim) # this should be used for external
    # representation
    primer.vars <- vars[seq_len(nrow(primer.df))]
    solution.idx <- which(primer.vars > 0.5)
    solution.size <- length(solution.idx)
    cur.primers <- primer.df[solution.idx, ]  # solution primers
    cur.cvg <- get_cvg_ratio(cur.primers, template.df)
    iter <- get.total.iter(ILP)
    val <- get.objective(ILP)
    nodes <- get.total.nodes(ILP)
    # data-associated properties
    max.temp.diff <- NA
    if ("melting_temp" %in% colnames(cur.primers)) {
        max.temp.diff <- max(abs(cur.primers$melting_temp - target.temp))
    }
    id <- paste(template.df$Group[1], "_", target.temp, sep = "")
    data <- data.frame(Identifier = id, Target_Temperature = target.temp, Iteration = i, 
        Set_Size = solution.size, Coverage = cur.cvg, Iterations = iter, Nodes = nodes, 
        Time = time, MaxTempDiffToTargetTemp = max.temp.diff,
        Objective = val, DeltaG_Cutoff = deltaG_Cutoff, DeltaG_Limit = deltaG_Limit, Variables = paste(vars, collapse = ","), 
        stringsAsFactors = FALSE)
    return(data)
}
#' Solve an ILP
#'
#' Constructs and solves an ILP and outputs a list with the reuslts.
#'
#' @param cur.D Binary dimerization matrix.
#' @param cur.G Free energy matrix for cross-dimerization.
#' @param cur.settings Current \code{DesignSettings} object.
#' @param deltaG.cutoff Cutoff for dimerization free energy.
#' @param deltaG.limit Relaxation limit for free energy cutoff.
#' @param cur.cvg.matrix Binary coverage matrix.
#' @param time.limit Time limit for solving the ILP in seconds.
#' @param required.cvg The target coverage of the designed primer set.
#' @param primer.df A \code{Primers} object.
#' @param template.df A \code{Templates} object.
#' @return List with ILP solution data.
#' @keywords internal
solve.ILP <- function(cur.D, cur.G, cur.settings, cur.cvg.matrix, 
                      time.limit, required.cvg, primer.df, template.df) {
    if (length(cur.cvg.matrix) == 0) {
        # no data to optimize
        return(NULL)
    }
    target.cvg <- min(get_cvg_ratio(primer.df, template.df), required.cvg)
    message("Solving ILP; Target cvg: ", target.cvg)
    deltaG.cutoff <- opti(cur.settings)$cross_dimerization
    deltaG.limit <- optiLimits(cur.settings)$cross_dimerization
    ILP <- ILPConstrained(cur.D, cur.cvg.matrix, time.limit)
    original.dim <- dim(ILP)
    return.val <- solve(ILP)  # presolve can also (even when only 'rows' is active, remove variables from the model)
    cvg.ratio <- 0
    if (return.val == 1) {
        # suboptimal solution (due to timeout)
        info.string <- paste("lpsolve: Stopped some of the iterations due to reaching the timeout.", 
            sep = "")
        warning(info.string)
    } else if (return.val != 0) {
        info.string <- paste("lpsolve could not find a solution. The return value was: ", 
            return.val, ". Please check the lpsolve documentation for more details.", 
            sep = "")
        if (return.val != 2) {
            stop(info.string)
        }
    } else if (return.val == 0) {
        vars <- get.ILP.vars(ILP, original.dim)
        primer.idx <- which(vars[seq_len(nrow(primer.df))] > 0.5)
        # get cvg ratio
        cur.cvg <- get_cvg_ratio(primer.df[primer.idx,], template.df)
    }
    initial.deltaG.cutoff <- deltaG.cutoff
    initial.deltaG.limit <- deltaG.limit
    relax.counter <- 0
    # keep on iterating until we're solveable OR we've reached the target cvg
    # OR cvg cannot be improved by changing dimerization cutoff
    while (return.val == 2 || cur.cvg < target.cvg && !all(cur.D == 0)) { 
        relax.counter <- relax.counter + 1
        # relax dimerization constraint
        deltaG.cutoff <- relax.opti.constraints(list("cross_dimerization" = deltaG.cutoff), list("cross_dimerization" = initial.deltaG.limit), 
                                                list("cross_dimerization" = initial.deltaG.cutoff))$cross_dimerization
        deltaG.limit <- relax.opti.constraints(list("cross_dimerization" = deltaG.limit), list("cross_dimerization" = initial.deltaG.limit), 
                                                list("cross_dimerization" = initial.deltaG.cutoff))$cross_dimerization
        message(paste("Setting deltaG.cutoff to ", deltaG.cutoff, sep = ""))
        cur.D <- compute.dimer.matrix(cur.G, deltaG.cutoff["min"])  # 1 for dimers, 0 else
        # update ILP with new dimerization values
        ILP <- ILPConstrained(cur.D, cur.cvg.matrix, time.limit)
        return.val <- solve(ILP)  # presolve can also (even when only 'rows' is active, remove variables from the model)
        original.dim <- dim(ILP)  # keep original.dim updated 
        if (return.val == 0) {
            # ILP never selects all primers, why? it should to reach cvg of all templates ..!  TODO
            vars <- get.ILP.vars(ILP, original.dim) 
            primer.idx <- which(vars[seq_len(nrow(primer.df))] > 0.5)
            cur.cvg <- get_cvg_ratio(primer.df[primer.idx,], template.df)
            message("Current cvg is: ", cur.cvg, ". Required cvg is: ", target.cvg)
        }
        #write.lp(ILP, file.path(cur.results.loc, paste0("LP_", relax.counter, ".lp")), 'lp') #write it to a file in LP format
    }
    result <- list(ILP = ILP, DeltaG_Cutoff = deltaG.cutoff["min"], DeltaG_Limit = deltaG.limit["min"], original_dim = original.dim)
    return(result)
}
#' Identification of Redudant Primers.
#'
#' Identifies primers that are redundant.
#'
#' Redundant primers do not reduce the coverage when removed.
#'
#' @param cvg.matrix Binary matrix of coverage events.
#' @return TRUE for redundant primers, FALSE otherwise.
#' @keywords internal
get.redundant.cols <- function(cvg.matrix) {
    template.coverage <- I.cvg(cvg.matrix)
    total.cvg <- length(unique(unlist(template.coverage)))
    is.redundant.col <- unlist(lapply(seq_along(template.coverage), function(x) length(unique(unlist(template.coverage[-x]))) == 
        total.cvg))
    return(is.redundant.col)
}
#' Removal of Redundant Primers.
#'
#' Removes redundant primers from an optimal solution. 
#'
#' An optimal solution can contain primers with redundant coverage when using presolve or greedy optimization.
#'
#' @param S Indices of primers that are selected in an optimal solution.
#' @param cvg.matrix Binary matrix of coverage events.
#' @return Updated indices of selected primers \code{S} such that indices
#' representing primers with redundant coverage are removed.
#' @keywords internal
remove.redundant.cols <- function(S, cvg.matrix) {
    redundant.cols <- get.redundant.cols(cvg.matrix[, S, drop = FALSE])  # boolean: redundant primers in solution
    
    R <- S[redundant.cols]  # index of redundant columns
    non.R <- setdiff(S, R)  # index of non-redundant columns
    if (length(R) == 0) {
        # no redundant columns -> no need to solve another SCP instance
        return(S)
    }
    M <- setdiff(unique(unlist(I.cvg(cvg.matrix[, R, drop = FALSE]))), unique(unlist(I.cvg(cvg.matrix[, 
        non.R, drop = FALSE]))))  # template indices covered by redundant cols only
    if (length(M) == 0) {
        # redundant primers don't afford any additional coverage (this shouldn't happen)
        S <- S[-R]  # remove all redundant primers        
    } else {
        cur.D <- matrix(rep(0, length(R) * length(R)), nrow = length(R))  # matrix without dimerization events
        # solve set cover for redundant primers and corresponding templates
        ILP <- ILPConstrained(cur.D, cvg.matrix[M, R, drop = FALSE], time.limit = Inf, 
                             presolve.active = FALSE) 
        return.val <- solve(ILP)
        vars <- lpSolveAPI::get.variables(ILP)
        S <- S[-R[vars == 0]]
    }
    return(S)
}
#' Solver for ILP Set Cover
#'
#' Solves the primer set cover problem using an ILP formulation.
#'
#' @param primer.df Primer data frame to be optimized.
#' @param template.df Template data frame with sequences.
#' @param settings A \code{DesignSettings} object.
#' @param primer_conc Primer concentration.
#' @param template_conc Template concentration.
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @param allowed.mismatches The number of mismatches primers are allowed to have with the templates.
#' @param allowed.other.binding.ratio Ratio of primers allowed to bind to non-target regions.
#' @param allowed.stop.codons  Consider mismatch binding events that induce stop codons.
#' @param allowed.region.definition Definition of the target binding sites used for evaluating the coverage.
#' If \code{allowed.region.definition} is \code{within}, primers have to lie within the allowed binding region.
#' If \code{allowed.region.definition} is \code{any}, primers have to overlap with the allowed binding region.
#' The default is that primers have to bind within the target binding region.
#' @param disallowed.mismatch.pos The number of positions from the primer 3' end where mismatches should not be allowed.
#' All primers binding templates with mismatches within \code{disallowed.mismatch.pos} from the 3' end are disregarded.
#' @param target.temps Target melting temperatures for primer sets in Celsius.
#' @param required.cvg Target coverage ratio of the templates by the primers.
#' @param fw.primers List with optimized primer data frames corresponding to \code{target.temps}. 
#' Only required for optimizing both strand directions and only 
#' in the second optimization run in order to check for cross dimerization.
#' @param diagnostic.location Directory for storing results.
#' @param timeout Timeout in seconds for the optimization with ILPs.
#' @param updateProgress Shiny progress callback function.
#' @return List with optimization results.
#' @keywords internal
optimize.ILP <- function(primer.df, template.df, settings, primer_conc, 
    template_conc, na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc, 
    allowed.mismatches, allowed.other.binding.ratio, allowed.stop.codons, 
    allowed.region.definition, disallowed.mismatch.pos, target.temps, 
    required.cvg, fw.primers = NULL, diagnostic.location = NULL, 
    timeout = Inf, updateProgress = NULL) {
    # todo: add updateprogress arg in server.R
    
    # sanity check nothing to optimize
    if (length(primer.df) == 0) {
        return(NULL)
    }
    #if (!"melting_temp" %in% colnames(primer.df)) {
        #stop("Melting temp has to be computed before optimizing")
    #}
    opti.constraints <- opti(settings)
    opti.limits <- optiLimits(settings)
    cvg.matrix <- get.coverage.matrix(primer.df, template.df)
    mode.directionality <- get.analysis.mode(primer.df)
    Tm.brackets <- create.Tm.brackets(primer.df, template.df, settings, target.temps)
    primer.df <- Tm.brackets$primers
    target.temps <- Tm.brackets$df$target_Tm
    #### COMPUTE CONSTRAINTS cross-dimerization
    if ("cross_dimerization" %in% names(opti.constraints)) {
        deltaG.cutoff <- opti.constraints$cross_dimerization["min"]
        deltaG.limit <- opti.limits$cross_dimerization["min"]
        if (length(deltaG.cutoff) == 0 || length(deltaG.limit) == 0) {
            stop("DeltaG constraint/limit was active but values were not specified.")
        }
        message("Computing cross dimers. This may take a while ...")
        G.df <- compute.all.cross.dimers(primer.df, primer_conc, na_salt_conc, 
                            mg_salt_conc, k_salt_conc, tris_salt_conc, 
                            min(Tm.brackets$df$annealing_temp),
                            no.structures = TRUE) # TODO: compute G matrix for every Tm bracket at some point, think about how we can relax the matrix then?
        G.matrix <- create.G.matrix(primer.df, G.df)  # min deltaG values of cross-dimerization conformations for every pair of primers
        # load G matrix for debugging: G.matrix <-
        # read.csv(file.path(diagnostic.location, 'G_matrix.csv'))
        if (length(diagnostic.location) != 0) {
            write.csv(G.matrix, file.path(diagnostic.location, "G_matrix.csv"), row.names = FALSE)
        }
        D <- compute.dimer.matrix(G.matrix, deltaG.cutoff)  # 1 for dimers, 0 else
    } else {
        # set dimerization D-matrix to zeros -> all primers are compatible
        D <- matrix(rep(0, nrow(primer.df) * nrow(primer.df)), nrow = nrow(primer.df), 
            ncol = nrow(primer.df))
        G.matrix <- matrix(rep(0, nrow(primer.df) * nrow(primer.df)), nrow = nrow(primer.df), 
            ncol = nrow(primer.df))
        deltaG.cutoff <- NULL
        deltaG.limit <- NULL
    }
    # determine annealing-temperature dependent constraints and filter according to them
    Tm.data <- compute.Tm.sets(primer.df, template.df, Tm.brackets, settings, 
                mode.directionality, primer_conc, template_conc, na_salt_conc, mg_salt_conc, 
                k_salt_conc, tris_salt_conc, allowed.mismatches, allowed.other.binding.ratio, 
                allowed.stop.codons, allowed.region.definition, disallowed.mismatch.pos, 
                TRUE, required.cvg, fw.primers, diagnostic.location, updateProgress)  
    Tm.sets <- Tm.data$sets
    Tm.settings <- Tm.data$settings # relaxed settings
    ######### 
    i <- NULL
    #for (i in seq_along(target.temps)) { # TODO: for debug
    ILP.df <- foreach(i = seq_along(target.temps), .combine = my_rbind) %dopar% {
            target.temp <- target.temps[i]
            Tm.set <- Tm.sets[[i]]
            cur.settings <- Tm.settings[[i]]
            if (nrow(Tm.set) == 0) {
                # create empty entry
                data <- build.ILP.df(NULL, NULL, Tm.set, template.df, i, target.temp)
            }
            # Tm.set <- read_primers(file.path(opti.results.loc,
            # paste('filtered_opti_data_target_temp_', target.temp, '.csv', sep = '')))
            #cat(paste("Optimization for target temp ", target.temp, " commences ...\n", 
                #sep = ""))
            sel.idx <- match(Tm.set$Identifier, primer.df$Identifier)  # index for accessing global properties of primer.df
            cur.G <- G.matrix[sel.idx, sel.idx, drop = FALSE]
            cur.D <- D[sel.idx, sel.idx, drop = FALSE]
            cur.cvg.matrix <- cvg.matrix[, sel.idx, drop = FALSE]
            time <- Sys.time()  # measure optimization time
            #### test sample.size <- 200 cur.G <- cur.G[1:sample.size,1:sample.size] Tm.set <-
            #### Tm.set[1:sample.size,]
            solution.data <- solve.ILP(cur.D, cur.G, cur.settings, 
                                       cur.cvg.matrix, timeout, required.cvg, Tm.set, template.df)
            if (length(solution.data) == 0 || class(solution.data) == "try-error") {
                # optimization was not possible for some reason (most likely due to dimerization
                # constraints)
                data <- build.ILP.df(NULL, NULL, Tm.set, template.df, i, target.temp)
            } else {
                vars <- get.ILP.vars(solution.data$ILP, solution.data$original_dim)
                primer.idx <- which(vars[seq_len(nrow(Tm.set))] > 0.5)
                primer.idx <- remove.redundant.cols(primer.idx, cur.cvg.matrix)
                vars[!seq_len(nrow(Tm.set)) %in% primer.idx] <- 0  # remove redundant cols from vars
                solution <- solution.data$ILP
                relaxed.deltaG.cutoff <- solution.data$DeltaG_Cutoff
                relaxed.deltaG.limit <- solution.data$DeltaG_Limit
                # for parallelization, need to extract the relevant infos from the ILP solutions,
                # because objects are destroyed when threads end (solutions are only ptrs to
                # lprec objects)
                time <- as.numeric(difftime(Sys.time(), time, units = "secs"))
                data <- build.ILP.df(solution, vars, Tm.set, template.df, i, target.temp, 
                  time = time, deltaG_Cutoff = relaxed.deltaG.cutoff, deltaG_Limit = relaxed.deltaG.limit)
            }
            data
    }
    solution.idx <- select.best.ILP(ILP.df)
    if (length(solution.idx) == 0) {
        warning("No optimal primers found. No primer set could be optimized with ILPs. Check your constraint settings (dimerization!).")
        return(NULL)
    }
    all.optimal.sets <- get.sets.from.decisions(ILP.df, Tm.sets)
    if (length(diagnostic.location) != 0) {
        # output info
        optimal <- rep(FALSE, nrow(ILP.df))
        optimal[solution.idx] <- TRUE
        lambda.idx <- which(colnames(ILP.df) == "Objective")
        ILP.df <- cbind(ILP.df[, 1:lambda.idx], Optimal = optimal, ILP.df[, ((lambda.idx + 
            1):ncol(ILP.df))], stringsAsFactors = FALSE)
        write.csv(ILP.df, file = file.path(diagnostic.location, "ILP_summary.csv"), 
            row.names = FALSE)
    }
    optimal.primers <- all.optimal.sets[[solution.idx]]
    relaxed.deltaG_Cutoff <- ILP.df$deltaG_Cutoff[solution.idx]
    # output the used opti constraints for each Tm set
    if ("cross_dimerization" %in% names(opti.constraints)) {
        # deltaG is the only relaxed opti constraint for ILPs
        for (i in seq_along(Tm.settings)) {
            if (!is.na(ILP.df$Identifier[i])) {
                constraintLimits(Tm.settings[[i]])$cross_dimerization["min"] <- ILP.df$DeltaG_Limit[i]
                constraints(Tm.settings[[i]])$cross_dimerization["min"] <- ILP.df$DeltaG_Cutoff[i]
            }
        }
    }
    result <- list(opti = optimal.primers, all_results = all.optimal.sets, 
        all_used_constraints = Tm.settings, 
        used_constraints = Tm.settings[[solution.idx]]) 
    return(result)
}
#' Optimal Sets from Decision Variables
#'
#' Determines primer sets from decision variables from ILP.
#'
#' @param ILP.df Data frame with ILP optimization results.
#' @param Tm.sets List with primer data frames for every target melting temperature.
#' @return A list with optimal primer data sets for every target temperature.
#' @keywords internal
get.sets.from.decisions <- function(ILP.df, Tm.sets) {
    results <- vector("list", length(Tm.sets))
    for (i in seq_along(Tm.sets)) {
        Tm.set <- Tm.sets[[i]]  # solution filtered data (with efficiencies!)
        primer.idx <- which(as.numeric(strsplit(as.character(ILP.df$Variables[i]), 
            split = ",")[[1]])[seq_len(nrow(Tm.set))] > 0.5)
        optimal.primers <- Tm.set[primer.idx, ]
        results[[i]] <- optimal.primers
    }
    names(results) <- ILP.df$Target_Temperature
    return(results)
}
#' Subset ILP Constructor
#'
#' Constructs an ILP for selecting optimal primer subsets.
#'
#' Here, "optimal" refers to a subset of a certain size that maximizes the coverage.
#'
#' @param primer.df Primer data frame to be subsetted.
#' @param template.df Template data frame.
#' @param k Required number of primers to be selected.
#' @return An ILP for choosing the primer subset of size \code{k} with the largest coverage.
#' @keywords internal
subset.ILP <- function(primer.df, template.df, k) {
    if (length(primer.df) == 0 || nrow(primer.df) == 0)  {
        return(NULL)
    }
    if (length(template.df) == 0 || nrow(template.df) == 0) {
        return(NULL)
    }
    if (k <= 0) {
        stop("Subset size 'k' should be positive.")
    }
    # template.coverage: for each template, the indices of covering primers
    template.coverage <- get.covered.templates(primer.df, template.df)
    covered.templates <- which(sapply(template.coverage, length) != 0)  # nbr of coverage constraints
    nbr.coverage.constraints <- length(covered.templates)
    xp <- nrow(primer.df)  # primer decision variables
    xt <- nbr.coverage.constraints  # template decision variables: only those templates that are covered are considered
    # we have two variables (x_i: primer selected, y_i: template covered)
    X <- xp + xt  # nbr of decision variables
    # no coverage constraint, instead constraint on the set size
    total.nbr.constraints <- 1 + nbr.coverage.constraints  # constraint for set size + template coverage constraint
    lprec <- make.lp(total.nbr.constraints, X)
    name.lp(lprec, "openPrimeR_SubsetSelection_ILP")
    for (col in 1:X) {
        # define primer vars as integer (bounds of 0,1)
        set.type(lprec, col, "binary")
    }
    coefs <- c(rep(0, xp), rep(1, xt))  # maximize the sum of templates that are covered
    set.objfn(lprec, coefs)
    # set ILP options
    lp.info <- lp.control(lprec, sense = "max", verbose = "neutral", 
        epslevel = "tight")
    rownames <- NULL
    # add constraints subset size constraint
    set.row(lprec, 1, c(rep(1, xp), rep(0, xt)))
    set.constr.type(lprec, "=", 1) # set should be of specified size
    set.rhs(lprec, k, 1)
    rownames(lprec)[1] <- "SetSize_Constraint"
    ####### coverage constraint: 
    #cat(paste("Adding ", length(covered.templates), " coverage constraints. \n", 
        #sep = ""))
    #if (length(covered.templates) > 0) {
        #pb <- txtProgressBar(min = 0, max = length(covered.templates), style = 3)
    #}
    for (j in seq_along(covered.templates)) {
        # select the idx of those primers that are covering the current template
        template.lp.idx <- j + xp
        idx <- template.coverage[[covered.templates[j]]]  # index of primers covering the j-th coverable template
        set.row(lprec, j + 1, c(-1, rep(1, length(idx))), indices = c(template.lp.idx, 
            idx))  # j+1 because we already have the coverage constraint in lprec row 1
        #setTxtProgressBar(pb, j)
    }
    #cat("\n")
    if (length(covered.templates) != 0) {
        set.constr.type(lprec, rep(">=", length(covered.templates)), 2:(nbr.coverage.constraints + 
            1))
        set.rhs(lprec, rep(0, length(covered.templates)), 2:(nbr.coverage.constraints + 
            1))
        rownames(lprec)[2:(nbr.coverage.constraints + 1)] <- "Coverage_Constraint"
    }
    cnames.p <- paste("Primer_", seq_len(nrow(primer.df)), sep = "")
    if (length(covered.templates) != 0) {
        cnames.t <- paste("Template_", covered.templates, sep = "")
    } else {
        cnames.t <- NULL
    }
    colnames(lprec) <- c(cnames.p, cnames.t)
    return(lprec)
}
#' Determination of Primer Coverage for Groups.
#'
#' Modifies a primer data frame to retain only coverage events
#' relating to the selected groups of templates.
#'
#' @param primer.df Primer data frame.
#' @param template.df Template data frame.
#' @param groups Template groups for which coverage should be determined.
#' @return \code{primer.df} with coverage considered only for \code{groups}.
#' @keywords internal
primer.coverage.for.groups <- function(primer.df, template.df, groups) {
    if (!is.null(groups)) {
        # filter template.df for groups
        template.df <- template.df[template.df$Group %in% groups, ]
    }
    idx <- covered.seqs.to.idx(primer.df$Covered_Seqs, template.df)
    cvg.counts <- sapply(idx, function(x) ifelse(length(x) == 0, 0, length(x[!is.na(x)]))) # na check
    names(cvg.counts) <- seq_len(nrow(primer.df))
    primer.df$primer_coverage <- cvg.counts
    # modify covered_seqs of primers
    template.IDs <- unlist(lapply(idx, function(x) {
        if (is.null(x)) {
            ""
        } else {
            paste(template.df$Identifier[x[!is.na(x)]], 
                    collapse = ",")
        }
   }))
    primer.df$Covered_Seqs <- template.IDs
    primer.df$Coverage_Ratio <- cvg.counts/nrow(template.df)
    ## remove all fields that weren't updated
    # if other fields are required -> re-evaluate the primers
    keep.cols <- c("Forward", "Reverse", "ID", "Identifier", 
                    "primer_length_fw", "primer_length_rev",
                    "Run", "primer_coverage", "Coverage_Ratio", 
                    "Covered_Seqs", "Direction")
    primer.df <- primer.df[, colnames(primer.df) %in% keep.cols]
    return(primer.df)
}
#' Optimal Primer Subsets. 
#' 
#' Determines subsets of the input primer set that 
#' are optimal with regard to the number of covered template 
#' sequences.
#'
#' The optimal subsets are identified by solving an integer-linear program.
#' Since the quality of the primers (in terms of properties) is not taken into
#' account when creating the subsets, this method should only be used
#' for primer sets that are already of high quality. 
#' 
#' @param primer.df An objectc of class \code{Primers} providing the
#' primers for which optimal subsets should be constructed.
#' @param template.df An object of class \code{Templates} providing the
#' template sequences that are targeted by \code{primer.df}.
#' @param k The spacing between generated primer subset sizes. By default,
#' \code{k} is set to 1 such that all primer subsets are constructed.
#' @param groups The identifiers of template groups according to which 
#' coverage should be determined. By default, \code{groups} is set to 
#' \code{NULL} such that all all covered templates are considered.
#' @param identifier An identifier for storing the primer set. By default,
#' \code{identifier} is set to  \code{NULL}.
#' @param cur.results.loc Directory for storing the results. By default,
#' \code{cur.results.loc} is set to \code{NULL}, which means that
#' no results are stored.
#' @return A list with optimal primer subsets, each of class \code{Primers}.
#' @export
#' @examples
#' data(Ippolito)
#' primer.subsets <- subset_primer_set(primer.df, template.df)
subset_primer_set <- function(primer.df, template.df, k = 1, groups = NULL, identifier = NULL, cur.results.loc = NULL) {
    if (k <= 0) {
        stop("k has to be positive ...")
    }
    if (length(primer.df) == 0 || nrow(primer.df) == 0)  {
        return(NULL)
    }
    if (length(template.df) == 0 || nrow(template.df) == 0) {
        return(NULL)
    }
    if (k > nrow(primer.df)) {
        stop("k shouldn't exceed size of primer set")
    }
    if (!is(primer.df, "Primers")) {
        stop("Please supply a valid primer data frame.")
    }
    if (!"primer_coverage" %in% colnames(primer.df)) {
        return(NULL)
    }
    if (!is(template.df, "Templates")) {
        stop("Please provide a valid template data frame.")    
    }
    # need either individual primers or only primer pairs for subsets to make sense
    if ((!all(primer.df$Forward == "") && any(primer.df$Reverse != "")) || (!all(primer.df$Reverse == 
        "") && any(primer.df$Forward != "")) || (any(primer.df$Forward != "") && 
        any(primer.df$Reverse != "") && !all(primer.df$Forward != "" & primer.df$Reverse != 
        ""))) {
            primer.df <- pair_primers(primer.df, template.df)
    }
    
    p.df <- primer.coverage.for.groups(primer.df, template.df, groups)
    out.loc <- file.path(cur.results.loc, identifier)
    if (length(out.loc) != 0) {
        dir.create(out.loc, showWarnings = FALSE)
    }
    K <- 1:100 * k
    K <- K[K <= nrow(primer.df)]
    # for (i in seq_along(K)) {
    k <- NULL
    top.primers <- foreach(k = seq_along(K), .combine = c) %dopar% {
        ILP <- subset.ILP(p.df, template.df, k)
        return.val <- solve(ILP)
        if (return.val != 0) {
            warning(return.val)
        }
        if (length(out.loc) != 0) {
            write.lp(ILP, file.path(out.loc, paste("ILP_", k, sep = "")), "lp")  #write it to a file in LP format
        }
        vars <- get.variables(ILP)
        primer.vars <- vars[seq_len(nrow(p.df))]
        sel.idx <- which(primer.vars == 1)
        cur.top.primers <- p.df[sel.idx, ]
        fname <- file.path(out.loc, paste("subset_k=", k, ".csv", sep = ""))
        if (length(fname) != 0) {
            write.csv(cur.top.primers, fname, row.names = FALSE)
        }
        list(cur.top.primers)
    }
    return(top.primers)
}
