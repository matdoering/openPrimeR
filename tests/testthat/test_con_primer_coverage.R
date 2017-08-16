context("Primer coverage")

test_that("Stop codon filter", {
    # test coverage stop codon filter
    data(Ippolito)
    template.df <- template.df[1:3,]
    # reading frame: aga|tac|...
    template.df$InputSequence <- c("agatacacccccccccccccc",
                                   "agatacagattacagattaca",
                                   "agataaagattacagattaca")
    template.df <- assign_binding_regions(template.df, fw = c(1,14), rev = c(1,5))
    stop.codons <- c("TAG", "TAA", "TGA")
    # first primer has gaTAGa: should bind to 1st, 2nd, 3rd with stop codon
    # second primer has gaTAAa: should bind to 1st, 2nd, 3rd with stop codon
    # third primer has gaTGA: should bind to the third template with stop codon
    seqs <- c("gataga", "gataaa", "gatgaa")
    fname <- tempfile("test.fasta")
    seqinr::write.fasta(as.list(seqs), as.list(c("testPrimer1_fw", "testPrimer2_fw", "testPrimer3_fw")), fname)
    primer.df <- read_primers(fname)
    # allow a single mismatch to introduce stop codons
    conOptions(settings)$allowed_mismatches <- 1
    # do not use any additional cvg constraints:
    cvg_constraints(settings) <- list()
    constraint.df <- check_constraints(primer.df, template.df, settings,
                        active.constraints = "primer_coverage")
#######
    #mutation.types <- c("stop_codon")
    #mutation.check <- mismatch.mutation.check(constraint.df, template.df, mutation.types)
###
    expect_that(constraint.df$primer_coverage,
                equals(c(3, 3, 1), tolerance = 0.0))
    expect_that(constraint.df$primer_specificity,
                equals(c(1, 1, 1), tolerance = 0.01))
    # now, add stop codon constraint to filter out stop codon events
    cvg_constraints(settings) <- list("stop_codon" = c("min" = 0, "max" = 0))
    constraint.df <- suppressWarnings(check_constraints(primer.df, template.df, settings,
                        active.constraints = "primer_coverage"))
    expect_that(constraint.df$primer_coverage,
                equals(c(0, 0, 0), tolerance = 0.0))
    # ensure that specificity is also updated after filtering for stop codons
    expect_that(constraint.df$primer_specificity,
                equals(c(0, 0, 0), tolerance = 0.01))
}) 
test_that("Free energy filter", {
    # check whether binding events are filtered correctly according to their free energies
    if (!check.tool.function()["OligoArrayAux"]) {
        # cannot test without OligoArrayAux
        skip("OligoArrayAux not available.")
    }
    data(Ippolito)
    primer.df <- primer.df[1:2,]
    # get baseline results if we don't filter
    cvg_constraints(settings) <- list("annealing_DeltaG" = c("max" = 0)) # no filtering
    constraint.df <- check_constraints(primer.df, template.df, settings,
                        active.constraints = "primer_coverage")
    # decide for a cutoff:
    deltaG.cutoff <- -15 # everything greater than -15 should be filtered
    # retrieve observed values
    deltaG.values <- lapply(strsplit(constraint.df$annealing_DeltaG, split = ","), as.numeric)
    sel.idx <- lapply(deltaG.values, function(x) which(x <= deltaG.cutoff))
    filtered.cvd.seqs <- sapply(seq_along(sel.idx), function(x) paste0(strsplit(constraint.df$Covered_Seqs[x], split = ",")[[1]][sel.idx[[x]]], collapse = ","))
    filtered.cvg <- sapply(sel.idx, length)
    # determine results from tool:
    cvg_constraints(settings) <- list("annealing_DeltaG" = c("max" = deltaG.cutoff)) # filter!
    constraint.df <- check_constraints(primer.df, template.df, settings,
                        active.constraints = "primer_coverage")
    expect_equal(filtered.cvg, constraint.df$primer_coverage) 
    expect_equal(filtered.cvd.seqs, constraint.df$Covered_Seqs)
})

test_that("Efficiency Filter", {
    # check whether binding events are filtered correctly according to their free energies
    if (!check.tool.function()["OligoArrayAux"]) {
        # cannot test without OligoArrayAux
        skip("OligoArrayAux not available.")
    }
    data(Ippolito)
    primer.df <- primer.df[1:2,]
    # get baseline results if we don't filter
    cvg_constraints(settings) <- list("primer_efficiency" = c("min" = 0)) # no filtering
    constraint.df <- check_constraints(primer.df, template.df, settings,
                        active.constraints = "primer_coverage")
    # decide for a cutoff:
    eff.cutoff <- 0.1 # everything smalelr than 0.1 efficiency should be filtered
    # retrieve observed values
    eff.values <- lapply(strsplit(constraint.df$primer_efficiency, split = ","), as.numeric)
    sel.idx <- lapply(eff.values, function(x) which(x >= eff.cutoff))
    filtered.cvd.seqs <- sapply(seq_along(sel.idx), function(x) paste0(strsplit(constraint.df$Covered_Seqs[x], split = ",")[[1]][sel.idx[[x]]], collapse = ","))
    filtered.cvg <- sapply(sel.idx, length)
    # determine results from tool:
    cvg_constraints(settings) <- list("primer_efficiency" = c("min" = eff.cutoff)) # filter
    constraint.df <- check_constraints(primer.df, template.df, settings,
                        active.constraints = "primer_coverage")
    expect_equal(filtered.cvg, constraint.df$primer_coverage) 
    expect_equal(filtered.cvd.seqs, constraint.df$Covered_Seqs)
})

test_that("3prime mismatches", {
    # check filtering for mismatches at the 3' ends
    data(Ippolito)
    primer.df <- primer.df[1:2,]
    # get baseline results if we don't filter
    cvg_constraints(settings) <- list("terminal_mismatch_pos" = c("min" = 0)) # no filtering
    conOptions(settings)$allowed_mismatches <- 3
    constraint.df <- check_constraints(primer.df, template.df, settings,
                        active.constraints = "primer_coverage")
    # decide for a cutoff:
    terminal.cutoff <- 6 # everything smalelr than 0.1 efficiency should be filtered
    # retrieve observed values
    terminal.pos <- lapply(strsplit(constraint.df$terminal_mismatch_pos, split = ","), as.numeric)
    sel.idx <- lapply(terminal.pos, function(x) which(x >= terminal.cutoff))
    filtered.cvd.seqs <- sapply(seq_along(sel.idx), function(x) paste0(strsplit(constraint.df$Covered_Seqs[x], split = ",")[[1]][sel.idx[[x]]], collapse = ","))
    filtered.cvg <- sapply(sel.idx, length)
    # determine results from tool:
    cvg_constraints(settings) <- list("terminal_mismatch_pos" = c("min" = terminal.cutoff)) # filter
    constraint.df <- check_constraints(primer.df, template.df, settings,
                        active.constraints = "primer_coverage")
    expect_equal(filtered.cvg, constraint.df$primer_coverage) 
    expect_equal(filtered.cvd.seqs, constraint.df$Covered_Seqs)
})

test_that("basic_coverage", {
    # check that basic reported coverage is correct:
        # position of binding
        # mismatch positions
        # number of mismatches
        # covered seqs
        # coverage count
    data(Ippolito)
    primer.df <- primer.df[1:2,]
    template.df <- template.df[template.df$Identifier %in% c(4,5,8,9)]
    cvg_constraints(settings) <- list() # no advanced cvg requirements
    conOptions(settings)$allowed_mismatches <- 0
    constraint.df <- check_constraints(primer.df, template.df, settings,
                        active.constraints = "primer_coverage")
    expect_equal(constraint.df$primer_coverage, c(2,2))
    expect_equal(constraint.df$Coverage_Ratio, c(0.5,0.5))
    expect_equal(constraint.df$Nbr_of_mismatches_fw, c("0,0", "0,0"))
    expect_equal(constraint.df$Binding_Position_Start_fw, c("58,58", "58,58"))
    expect_equal(constraint.df$Binding_Position_End_fw, c("80,80", "80,80"))
    cvg_constraints(settings) <- list("terminal_mismatch_pos" = c(min = 6))
    constraint.df <- check_constraints(primer.df, template.df, settings,
                        active.constraints = "primer_coverage")

})
