context("Primers") 

test_that("Primer CSV Input", {
    # tests whether csv input/output of primers works
    data(Ippolito)
    primer.csv <- tempfile("my_primers", fileext = ".xml")
    rownames(primer.df) <- NULL
    write.csv(primer.df, file = primer.csv, row.names = FALSE)
    csv.primers <- read_primers(primer.csv)
    expect_that(csv.primers, equals(primer.df))
    # check for error behavior on loading incorrect primers
    test.df <- data.frame(A = 1:3, B = 6:8)
    write.csv(test.df, file = primer.csv, row.names = FALSE)
    expect_error(read_primers(primer.csv))
})

test_that("Primer Set Significance", {
    # check that significance computations provide reasonable results
    data(Ippolito)
    primer.location <- system.file("extdata", "IMGT_data", "primers", "IGHV", 
                        "Ippolito2012.fasta", package = "openPrimeR")
    primers <- read_primers(primer.location, "_fw", "_rev")
    expect_warning(primer_significance(primers)) # nothing evaluated -> warning
    conOptions(settings)$allowed_mismatches <- 0
    constraint.df <- suppressWarnings(check_constraints(primers, template.df, settings,
                        active.constraints = "primer_coverage"))
    expect_warning(primer_significance(constraint.df))
    expect_warning(primer_significance(constraint.df, active.constraints = "gc_ratio"))
    p <- primer_significance(primer.df) # should be significant
    expect_lt(as.numeric(p), 1e-3)

})
