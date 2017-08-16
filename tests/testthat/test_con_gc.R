context("GC clamp and GC ratio")

test_that("GC Clamp", {
    # test primers for GC clamp
    expect_that(evaluate.GC.clamp("CAGAC"), equals(1))
    expect_that(evaluate.GC.clamp("CAGACC"), equals(2))
    expect_that(evaluate.GC.clamp("CAGACCG"), equals(3))
    expect_that(evaluate.GC.clamp("CAGACCGA"), equals(0))
}) 

test_that("GC Ratio", {
    # test GC ratios of sequences
    expect_that(compute.gc.ratio("CGACTACG"), equals(5/8))
    expect_that(compute.gc.ratio("ACG"), equals(2/3))
    expect_that(compute.gc.ratio("ATT"), equals(0))
    # test ambiguities: mean of disambiguated seqs
    seq <- "CGAGTAGTMC" # M -> A/C
    expect <- (5/10 + 6/10) / 2
    expect_that(compute.gc.ratio(seq), equals(expect))
})
