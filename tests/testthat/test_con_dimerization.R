context("Dimerization")

test_that("Self Dimerization", {
    # test self dimerization
    if (!check.tool.function()["OligoArrayAux"]) {
        # cannot test without oligo
        skip("OligoArrayAux not available.")
    }
    data(Ippolito)
    # select a primer for testing:
    # fix the annealing temp
    PCR(settings)$annealing_temp <- 57
    sel <- primer.df$Forward == "caggtgcagctgcaggagtcsg"
    constraint.df <- check_constraints(primer.df[sel,], template.df, 
                 settings, active.constraints = "self_dimerization")
    # worst-case binding mode: 12 maches, 6 mismatches
    #        5'-CAGGTGCAGCTGCAGGAGTCSG-3'
    #           |XXX||||||||||XXX|
    #    3'-GSCTGAGGACGTCGACGTGGAC-5'
    # free energy
    expect_that(constraint.df[,"Self_Dimer_DeltaG"], 
                equals(-6.996, tolerance = 0.5))
    # add a reverse primer and check worst-case binding again
    primer.df[sel, "Reverse"] <- "ctcctgcagctgcaggagtcsg"
    # new binding mode:
    #     5'-CTCCTGCAGCTGCAGGAGTCSG-3'</b><br>
    #        ||||||||||||||||||
    # 3'-GSCTGAGGACGTCGACGTCCTC-5'</b>
    constraint.df <- check_constraints(primer.df[sel,], template.df, 
                    settings, active.constraints = "self_dimerization")
    # match score should increase by 6 since we introduced 3 matches 
    # and prevented 3 mismatches
    # free energy
    expect_that(constraint.df[,"Self_Dimer_DeltaG"], 
                equals(-13.028, tolerance = 0.5))
}) 

test_that("Cross Dimerization", {
    # test cross dimerization
    if (!check.tool.function()["OligoArrayAux"]) {
        # cannot test without oligo
        skip("OligoArrayAux not available.")
    }
    data(Ippolito)
    constraint.df <- check_constraints(primer.df, template.df, 
                 settings, active.constraints = "cross_dimerization")
    sel <- "gaggtgcagctgktggagwcy"
    test.df <- constraint.df[which(constraint.df$Forward == sel),]
    #                 5'-GAGGTGCAGCTGKTGGAGWCY-3'
    #                    X||X||||X|||
    #       3'-GSCTGAGGACGTCGACGTGGAC-5'
    expect_that(test.df[, "Cross_Dimer_DeltaG"], 
                equals(-6.775, tolerance = 0.5))
    prev.G <- test.df[, "Cross_Dimer_DeltaG"]
    # modify cross-dimer partner reverse
    constraint.df[which(constraint.df$ID == test.df$Cross_Dimer_Idx1), "Reverse"] <-
                "cagctgcagctgcaggagtcsg"
    constraint.df <- check_constraints(constraint.df, template.df, 
                 settings, active.constraints = "cross_dimerization")
    test.df <- constraint.df[which(constraint.df$Forward == sel),]
    # energy should be smaller since this conformation has fewer mismatches
    expect_lte(test.df[, "Cross_Dimer_DeltaG"], prev.G)
})

