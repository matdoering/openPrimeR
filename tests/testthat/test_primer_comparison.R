context("Primer Comparison") 

test_that("primer_comparison", {
    # check that the filtering procedure excludes the correct primers
    data(Comparison)
    primer.data <- primer.data[1:3]
    template.data <- template.data[1:3]
    constraint.settings <- constraints(settings)
    constraints <- names(constraint.settings)
    # create basic plots etc.
    p <- plot_constraint_deviation(primer.data, settings)
    p <- plot_primer_binding_regions(primer.data, template.data)
    p <- plot_constraint(primer.data, settings, constraints)
    p <- plot_cvg_constraints(primer.data, settings)
    p <- plot_constraint_fulfillment(primer.data, settings)
    p <- plot_template_cvg(primer.data, template.data)
    tab <- get_comparison_table(template.data, primer.data)
    ######
    # test what happens if a template set is lacking required objects
    #######
    # missing coverage
    primer.data[[1]] <- primer.data[[1]][, !colnames(primer.data[[1]]) %in% c("Covered_Seqs", "primer_coverage")]
    p <- expect_warning(plot_template_cvg(primer.data, template.data)) # ok :-) 
    p <- expect_error(plot_primer_binding_regions(primer.data, template.data), ".*primer coverage.*")
    tab <- expect_warning(get_comparison_table(template.data, primer.data)) # ok 
    # missing constraint
    primer.data[[1]] <- primer.data[[1]][, !colnames(primer.data[[1]]) %in% c("melting_temp")]
    p <- expect_warning(plot_constraint_deviation(primer.data, settings)) # ok
    p <- expect_warning(plot_constraint(primer.data, settings, "melting_temp_range")) # ok
    p <- expect_warning(plot_constraint_fulfillment(primer.data, settings))
})
