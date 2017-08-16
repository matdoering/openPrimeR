context("Repeats and Runs")

test_that("Number of repeats", {
    # test dinucleotide repeats
    expect_that(nbr.of.repeats("ATATATAT"), equals(4))
    expect_that(nbr.of.repeats("TCGCGCG"), equals(3))
    expect_that(nbr.of.repeats("T"), equals(0))
    expect_that(nbr.of.repeats("TACG"), equals(1))
    expect_that(nbr.of.repeats(5), equals(NULL))
}) 

test_that("Number of runs", {
    # test single nucleotide runs
    expect_that(nbr.of.runs("TAAAAGG"), equals(4))
    expect_that(nbr.of.runs("ACTG"), equals(1))
    expect_that(nbr.of.runs(numeric(0)), equals(NULL))
})
