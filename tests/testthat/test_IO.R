context("Input/Output")

test_that("Settings IO Consistency", {
    # load settings
    filename <- system.file("extdata", "settings", 
                  "C_Taq_PCR_high_stringency.xml", package = "openPrimeR")
    settings <- read_settings(filename)
	settings.file <- tempfile("test_settings", fileext = ".xml")
    write_settings(settings, settings.file)
    # verify that test_settings agree with settings:
    test.settings <- read_settings(settings.file)
    for (i in seq_along(settings)) {
        slot.name <- slotNames(settings)[i]
        setting.slots <- slot(settings, slot.name)
        test.slots <- slot(test.settings, slot.name)
        expect_that(length(setting.slots), equals(length(test.slots)))
        expect_that(names(setting.slots), equals(names(test.slots)))
        expect_that(setting.slots, equals(test.slots))
    }
}) 
 
