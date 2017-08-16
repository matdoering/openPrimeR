context("Templates") 

test_that("Template CSV Input", {
    # tests whether csv input/output of templates works
    data(Ippolito)
    # test whether correct template files are loaded
    template.files <- system.file("extdata", "IMGT_data", "comparison", "templates", 
                              "IGH_templates.csv", package = "openPrimeR")
    template.data <- read_templates(template.files)
    # test whether we can store personal templates
    template.csv <- tempfile("my_templates", fileext = ".csv")
    write.csv(template.df, file = template.csv, row.names = FALSE)
    # reset rownames for test, would be different due to subsetting otherwise
    rownames(template.df) <- NULL
    template.data <- read_templates(template.csv)
    expect_that(template.data, equals(template.df))
    # removed some tests checking for errors when reading since these make problems with testthat environment ..
})
 
