context("Primer Design") 

test_that("primer_filtering_constraints", {
    # check that the filtering procedure excludes the correct primers
    data(Ippolito)
    constraints(settings)$primer_length <- c("min" = 18, "max" = 22)
    constraints(settings)$gc_clamp <- c("min" = 1, "max" = 2)
    constraints(settings)$gc_ratio <- c("min" = 0.4, "max" = 0.6)
    constraints(settings)$no_runs <- c("min" = 0, "max" = 3)
    constraints(settings)$no_repeats <- c("min" = 0, "max" = 3)
    constraints(settings)$melting_temp_range <- c("min" = 55, "max" = 65)
    constraints(settings)$melting_temp_diff <- c("min" = 0, "max" = 5)
    constraints(settings)$self_dimerization <- c("min" = -5)
    constraints(settings)$secondary_structure <- c("min" = -2)
    constraints(settings)$primer_coverage <- c("min" = 1)
    constraints(settings)$primer_specificity <- c("min" = 0.9, "max" = 1)
    constraints(settings)$cross_dimerization <- c("min" = -5)
    conOptions(settings)$allowed_mismatches <- 5
    cvg_constraints(settings) <- list("annealing_DeltaG" = c("max" = -5))
    # evaluate constraints first to have reference values
    active.constraints <- c("primer_length", "gc_clamp", "gc_ratio", 
                            "no_runs", "no_repeats", 
                            "primer_coverage", "primer_specificity")
    #constraint.df <- check_constraints(primer.df, template.df, settings)
    filter.result <- suppressWarnings(cascaded.filter.quick(primer.df, template.df, settings,
                        to.compute.constraints = active.constraints))
    excluded.df <- filter.result$excluded
    # filtered by primer length:
    filter.reasons <- rep(NA, nrow(primer.df))
    filter.reasons[1:2] <- "primer_length"
    filter.reasons[3:6] <- "gc_clamp"
    filter.reasons[7] <- "gc_ratio"
    filter.reasons[8] <- "primer_specificity"
    expect_equal(filter.reasons, excluded.df$Exclusion_Reason)
})
test_that("primer_degeneration", {
    primers <- c("ctccaaggt", "ktccaaggt",
                 "ntccaaggt", "nktccaggt",
                 "nkwsccaggt", "nkw-svhdb")
    degens <- c(1, 2, 4, 4 * 2, 4*2*2*2, 4*2*2*2*3*3*3*3)
    score <- score_degen(strsplit(primers, split = ""), gap.char = "-")
    expect_equal(score, degens)
})
test_that("primer_initialization_naive", {
    data(Ippolito)
    allowed.fw <- c(20,50)
    allowed.rev <- c(10,40)
    template.df <- assign_binding_regions(template.df, fw = allowed.fw, rev = allowed.rev)
    # nb: position check may fail by accident if primer binds earlier/later in some sequence by chance than intended ..
    # -> don't change selected templates!
    template.df <- template.df[c(1, 3),]
    # naive primer initialization fw
    primer.lengths <- 18:20
    primers <- create.initial.primer.set(template.df, primer.lengths, 
                mode.directionality = "fw", "test", allowed.region.definition = "within",
                init.algo = "naive", 16, 1)
    # ensure that all primers are unique
    expect_equal(length(unique(primers)), length(primers))
    # ensure that primer identifiers are unique
    expect_equal(length(unique(names(primers))), length(names(primers)))
    # check that primers match templates WITHIN the target region
    for (i in seq_along(primers)) {
        x <- primers[i] 
        hits <- regexpr(x, template.df$Sequence)
        match.len <- attr(hits, "match.length")
        hit.len.idx <- which(match.len != -1)
        # require at least one hit per primer:
        expect_gte(length(hit.len.idx), 1) 
        # check length of complementary region
        expect_lte(max(match.len[hit.len.idx]), max(primer.lengths))
        expect_gte(min(match.len[hit.len.idx]), min(primer.lengths))
        # check whether primer binds in allowed region
        binding.pos.s <- hits[hit.len.idx] # start of binding
        binding.pos.e <- hits[hit.len.idx] + (match.len[hit.len.idx] - 1) # end of binding
        expect_gte(min(binding.pos.s), allowed.fw[1])
        expect_lte(max(binding.pos.e), allowed.fw[2])
    }
    # check for "any" positions
    primers <- create.initial.primer.set(template.df, primer.lengths, 
                mode.directionality = "fw", "test", allowed.region.definition = "any",
                init.algo = "naive", 16, 1)
    min.pos.s <- 1000
    max.pos.e <- -1
    for (i in seq_along(primers)) {
        x <- primers[i] 
        hits <- regexpr(x, template.df$Sequence)
        match.len <- attr(hits, "match.length")
        hit.len.idx <- which(match.len != -1)
        # require at least one hit per primer:
        expect_gte(length(hit.len.idx), 1) 
        # check length of complementary region
        expect_lte(max(match.len[hit.len.idx]), max(primer.lengths))
        expect_gte(min(match.len[hit.len.idx]), min(primer.lengths))
        # check whether primer intersects with allowed region
        binding.pos.s <- hits[hit.len.idx] # start of binding
        binding.pos.e <- hits[hit.len.idx] + (match.len[hit.len.idx] - 1) # end of binding
        expect_lte(max(binding.pos.s), allowed.fw[2]) # don't overextend
        if (min(binding.pos.s) < min.pos.s) {
            min.pos.s <- min(binding.pos.s)
        }
        if (max(binding.pos.e) > max.pos.e) {
            max.pos.e <- max(binding.pos.e)
        }
    }
    expect_lt(min.pos.s, allowed.fw[1]) # there should be at least one primer binding earlier
    expect_gt(max.pos.e, allowed.fw[2]) # there should be at least one primer binding later
    # naive rev
    primers <- create.initial.primer.set(template.df, primer.lengths, 
                mode.directionality = "rev", "test", allowed.region.definition = "within",
                init.algo = "naive", 16, 1)
    for (i in seq_along(primers)) {
        x <- primers[i] 
        hits <- regexpr(x, rev.comp.sequence(template.df$Sequence))
        match.len <- attr(hits, "match.length")
        hit.len.idx <- which(match.len != -1)
        # require at least one hit per primer:
        expect_gte(length(hit.len.idx), 1) 
        # check length of complementary region
        expect_lte(max(match.len[hit.len.idx]), max(primer.lengths))
        expect_gte(min(match.len[hit.len.idx]), min(primer.lengths))
        # check whether primer binds in allowed region
        binding.pos.s <- hits[hit.len.idx] # start of binding
        binding.pos.e <- hits[hit.len.idx] + (match.len[hit.len.idx] - 1) # end of binding
        expect_gte(min(binding.pos.s), allowed.rev[1])
        expect_lte(max(binding.pos.e), allowed.rev[2])
    }
    # check for MAFFT:
    if (!check.tool.function()["MAFFT"]) {
        # cannot test without OligoArrayAux
        skip("MAFFT not available.")
    }
    # test tree primer creation
    primer.lengths <- 30
    data(Ippolito)
    template.df <- assign_binding_regions(template.df, fw = allowed.fw, rev = allowed.rev)
    max.degen <- 1
    idx <- 50:52
})


test_that("full_design_function", {
    data(Ippolito)
    template.df <- template.df[1:5,]
    constraints(settings)$primer_length[2] <- 18
    # test design without melting temp range but with melting temp diff:
    constraints(settings) <- constraints(settings)[names(constraints(settings)) != "melting_temp_range"]
    optimal.primers.greedy <- design_primers(template.df, "both", settings)
    # check that melting temp diff and coverage of the optimal set are ok
    opti.set <- optimal.primers.greedy$opti
    cvg.ratio <- as.numeric(get_cvg_ratio(opti.set, template.df))
    expect_equal(cvg.ratio, 1.0) # 100% cvg?
    # check Tm diff, should be small ...
    expect_lt(max(opti.set$melting_temp_diff), 5)
    # test without any melting temp
    constraints(settings) <- constraints(settings)[names(constraints(settings)) != "melting_temp_diff"]
    optimal.primers.greedy <- design_primers(template.df, "both", settings)
    # coverage ok?
    opti.set <- optimal.primers.greedy$opti
    cvg.ratio <- as.numeric(get_cvg_ratio(opti.set, template.df))
    expect_equal(cvg.ratio, 1.0) # 100% cvg?
    # only one result?
    expect_equal(length(optimal.primers.greedy$all_results), 1)
    # design with no constraints
    constraints(settings) <- list()
    optimal.primers.greedy <- try(design_primers(template.df, "both", settings), silent = TRUE)
    # shouldn't be possible to design without any constraints
    expect_that(class(optimal.primers.greedy), equals("try-error"))
    # but it should be possible with just primer length & coverage
    constraints(settings) <- list("primer_coverage" = c("min" = 1),
                                  "primer_length" = c("min" = 18, "max" = 18))
    optimal.primers.greedy <- design_primers(template.df, "both", settings)
    expect_equal(length(optimal.primers.greedy$all_results), 1.0)
    opti.set <- optimal.primers.greedy$opti
    cvg.ratio <- as.numeric(get_cvg_ratio(opti.set, template.df))
    expect_equal(cvg.ratio, 1.0) # 100% cvg?
    # test what happens if we forbid any relaxation:
    data(Ippolito)
    template.df <- template.df[1:5,]
    constraintLimits(settings) <- constraints(settings)
    # no relaxation should happen in this case ... coverage should be lower
    optimal.primers.greedy <- design_primers(template.df, "both", settings)
    opti.set <- optimal.primers.greedy$opti
    cvg.ratio <- as.numeric(get_cvg_ratio(opti.set, template.df))
    expect_lte(cvg.ratio, 1.0) # less/equal to 100% cvg?
    # test whehther we are also happy with a lower target cvg if specified so
    data(Ippolito)
    template.df <- template.df[1:5,]
    # ramp up the constraints a bit to ensure we don't get a too high cvg
    constraints(settings)$no_runs <- c("min" = 2, "max" = 2)
    constraints(settings)$no_repeats <- c("min" = 2, "max" = 2)
    constraints(settings)$gc_clamp <- c("min" = 1, "max" = 1)
    optimal.primers.greedy <- design_primers(template.df, "both", settings,
                              required.cvg = 0.2)
    opti.set <- optimal.primers.greedy$opti
    expect_gte(as.numeric(get_cvg_ratio(opti.set, template.df)), 0.2)
    # test this with ILP as well 
    optimal.primers.ILP <- design_primers(template.df, "both", settings,
                              required.cvg = 0.2, opti.algo = "ILP")
    opti.set <- optimal.primers.ILP$opti
    expect_gte(as.numeric(get_cvg_ratio(opti.set, template.df)), 0.2)
})

