context("Primer length")

test_that("Primer length", {
    # test primer lengths
    data(Ippolito)
    expect_that(nchar(primer.df$Forward), equals(primer.df$primer_length_fw))
    # no rev primers: length should be 0
    expect_that(nchar(primer.df$Reverse), equals(rep(0, nrow(primer.df))))
}) 

