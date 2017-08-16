context("Primer specificity")

test_that("primer_specificity", {
    data(Ippolito)
    template.df <- template.df[1:3,]
    # 1:13 binding, 14:... off
    template.df$InputSequence <- c("gatacacccccccccccccc",
                                   "gatacagattacagattaca",
                                   "gatacagattacagattaca")
    template.df <- assign_binding_regions(template.df, fw = c(1,13), rev = c(1,5))
    # create a new primer for testing
    seqs <- c("gataca", "gattaca", "ttaagaaaagga")
    fname <- tempfile("test.fasta")
    seqinr::write.fasta(as.list(seqs), as.list(c("testPrimer1_fw", "testPrimer2_fw", "testPrimer3_fw")), fname)
    primer.df <- read_primers(fname)
    conOptions(settings)$allowed_mismatches <- 0
    cvg_constraints(settings) <- list()
    constraint.df <- check_constraints(primer.df, template.df, settings,
                        active.constraints = "primer_coverage")
    expect_that(constraint.df$primer_coverage,
                equals(c(3, 2, 0), tolerance = 0.0))
    expect_that(constraint.df$primer_specificity,
                equals(c(1, 0.5, 0), tolerance = 0.01))
    # check coverage/specificity with mismatch binding
    seqs <- c("gaccca", "gacccca", "ttaagaaaagga") # 2 and 3 mismatches for binding required
    fname <- tempfile("test.fasta")
    seqinr::write.fasta(as.list(seqs), as.list(c("testPrimer1_fw", "testPrimer2_fw", "testPrimer3_fw")), fname)
    primer.df <- read_primers(fname)
    conOptions(settings)$allowed_mismatches <- 3
    constraint.df <- check_constraints(primer.df, template.df, settings,
                        active.constraints = "primer_coverage")
    expect_that(constraint.df$primer_coverage,
                equals(c(3, 3, 0), tolerance = 0.0))
    expect_that(constraint.df$primer_specificity,
                equals(c(3/(3+3), 3/(3+3), 0), tolerance = 0.01))
    # test filtering of cvg wrt specificity:
    # this part requires efficiency/annealing computations with oligoArrayAux
    if (!check.tool.function()["OligoArrayAux"]) {
        # cannot test without oligo
        skip("OligoArrayAux not available.")
    }
    data(Ippolito)
    conOptions(settings)$allowed_mismatches = 5
    seqs <- c("gatacaga")
    fname <- tempfile("test.fasta")
    seqinr::write.fasta(as.list(seqs), as.list(c("testPrimer1_fw")), fname)
    cvg_constraints(settings) <- list("annealing_DeltaG" = c("max" = (-10))) # select only bona-fide events
    primer.df <- read_primers(fname)
    constraint.df <- check_constraints(primer.df, template.df, settings,
                        active.constraints = "primer_coverage")
    low.spec <- constraint.df$primer_specificity
    expect_equal(low.spec, 1)
    expect_lt(length(unlist(strsplit(constraint.df$Off_Covered_Seqs, split = ","))), 
              length(unlist(strsplit(constraint.df$Basic_Off_Covered_Seqs, split = ","))))
}) 
