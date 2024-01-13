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
    expect_warning(primer_significance(primer.df, active.constraints = "primer_coverage"))
    p <- primer_significance(primer.df) # should not be significant
    expect_gt(as.numeric(p), 5e-2)
})

test_that("Primer fw+rev", {
    # check that primer sets with pairs of primers are analyzed correctly
    data(Ippolito)
    primer.location <- system.file("extdata", "IMGT_data", "primers", "Test", 
                        "test_both.fasta", package = "openPrimeR")
    primers <- read_primers(primer.location, "_fw", "_rev")
    constraint.df <- check_constraints(primers[1:2,], template.df, settings)
    temp.file <- tempfile(fileext = ".pdf")
    # check whether a report can be generated
    create_report(constraint.df, template.df, temp.file, settings)
})
