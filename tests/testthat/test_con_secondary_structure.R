context("Secondary Structures")

test_that("Secondary Structure", {
    # test self dimerization
    if (!check.tool.function()["ViennaRNA"]) {
        # cannot test without viennaRNA
        skip("Secondary structure tests require ViennaRNA.")
    }
    data(Ippolito)
    # select a primer for testing:
    # fix the annealing temp
    # test whether automatic setting of Ta works: structure energy should be different for the higher automatic Ta than the one we will set next:
    sel <- primer.df$Forward == "caggtgcagctgcaggagtcsg"
    constraint.df <- check_constraints(primer.df[sel,], template.df, 
                 settings, active.constraints = "secondary_structure")
    expect_that(constraint.df[, "Structure_deltaG"], 
                equals(0.0, tolerance = 0.01))
    # let's use a fixed Ta for the next computations: structure energy should be different
    PCR(settings)$annealing_temp <- 57
    sel <- primer.df$Forward == "caggtgcagctgcaggagtcsg"
    constraint.df <- check_constraints(primer.df[sel,], template.df, 
                 settings, active.constraints = "secondary_structure")
    expect_that(constraint.df[, "Structure_deltaG"], 
                equals(-0.12, tolerance = 0.01))
    expect_that(constraint.df[, "Structure_deltaG_fw"], 
                equals(constraint.df[, "Structure_deltaG"],
                tolerance = 0.01))

    # add a reverse primer with considerable structure formation
    primer.df[sel, "Reverse"] <- "cccccccaaaaagggggtttttttttt"
    constraint.df <- check_constraints(primer.df[sel,], template.df, 
                    settings, active.constraints = "secondary_structure")
    expect_that(constraint.df[,"Structure_deltaG"], 
                equals(-6.77, tolerance = 0.1))
    expect_that(constraint.df[,"Structure_deltaG"], 
                equals(constraint.df[, "Structure_deltaG_rev"], 
                tolerance = 0.1))
}) 
