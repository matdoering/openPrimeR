context("Settings") 

test_that("Settings validity", {
    # check whether settings validity is detected correctly
    data(Ippolito)
    # check that settings can be changed without errors
    constraints(settings)[["gc_clamp"]] <- c("min" = 0, "max" = 5)
    constraints(settings) <- constraints(settings)[1:4]
    constraintLimits(settings) <- constraintLimits(settings)[1:4]
    PCR(settings)$Na_concentration <- 0.00001
    conOptions(settings)$allowed_mismatches <- 2
    # check required fields -> should produce errors
    expect_error({PCR(settings) <- PCR(settings)[-which(names(PCR(settings)) == "Na_concentration")]})
    expect_error({conOptions(settings) <- conOptions(settings)[-which(names(conOptions(settings)) == "allowed_mismatches")]})
    expect_error({constraints(settings)$gc_clamp <- c("min" = 64, "max" = 30)})
    expect_error({constraints(settings)$gc_clamp <- c("max" = 30, "min" = 40)})
    test <- NULL
})

test_that("Settings concordance", {
    # check whether settings adjust when constraints/constraintLimit changes
    data(Ippolito)
    # widen the constraints beyond the limits:
    constraints(settings)$gc_clamp[2] <- 20
    new.setting <- constraints(settings)$gc_clamp[2]
    expect_that(constraintLimits(settings)$gc_clamp[2], equals(new.setting))
    # narrow down the limits beyond the constraint
    new.setting <- c("min" = 0, "max" = 0)
    constraintLimits(settings)$no_runs <- new.setting
    expect_that(constraints(settings)$no_runs, equals(new.setting))
    # check whether removal keeps constraints/limits concordant
    constraints(settings) <- constraints(settings)[names(constraints(settings)) %in% c("gc_clamp", "no_runs")]
    expect_that(names(constraints(settings)), equals(names(constraintLimits(settings))))
    # the same for addition of constraints
    constraintLimits(settings)$primer_length <- c("min" = 18, "max" = 30)
    expect_that(constraints(settings)$primer_length, equals(constraintLimits(settings)$primer_length))
})


test_that("Settings IO", {
    # check whether settings are stored/loaded correctly
    settings.xml <- system.file("extdata", "settings", 
                    "C_Taq_PCR_high_stringency.xml", package = "openPrimeR")
    settings <- read_settings(settings.xml)
    constraints(settings)$melting_temp_diff <- c("max" = 3)
    out.xml <- tempfile("my_settings", fileext = ".xml")
    write_settings(settings, out.xml)
    loaded.settings <- try(read_settings(out.xml), silent = TRUE)
    expect_that(class(loaded.settings), equals(class(settings)))
    expect_that(length(slotNames(loaded.settings)), equals(length(slotNames(settings))))
    for (i in seq_along(slotNames(loaded.settings))) {
        loaded.name <- slotNames(loaded.settings)[i]
        old.name <- slotNames(settings)[i]
        # ensure that slot names are identical
        expect_that(loaded.name, equals(old.name))
        loaded <- constraints(slot(loaded.settings, loaded.name))
        old <- constraints(slot(settings, loaded.name))
        # check that list entries are identical
        for (j in seq_along(loaded)) {
            con.name <- names(loaded)[j]
            loaded.val <- loaded[[j]]
            old.val <- old[[con.name]]
            expect_that(loaded.val, equals(old.val))
        }
    }
})
