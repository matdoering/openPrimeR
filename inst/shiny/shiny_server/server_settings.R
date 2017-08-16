##########
# Shiny server functionalities for tool settings
############
CoreObserver <- observeEvent(input$no_of_cores, {
    # the used number of cores by the tool.
    doParallel.available <- requireNamespace("doParallel", quietly = TRUE)
	is.win <- grepl("windows", .Platform$OS.type)
    if (doParallel.available && !is.win) { # no parallel support for m$
        doParallel::registerDoParallel(cores = min(input$no_of_cores, parallel::detectCores()))
    } else {
        # doParallel isn't available
        # disable slider and don't register the parallel backend
		if (is.win) {
			warning("parallel backend not available for windows.")
		} else {
			warning("doParallel not available: using 1 core only.")
		}
        shinyjs::disable("no_of_cores")
    }
}, ignoreNULL = TRUE)


mismatchPreventObserver1 <- observeEvent(c(input$allow_3prime_mismatch), {
    if (input$allow_3prime_mismatch == "active") {
        # activate 3 prime mismatch input and set to 0
        shinyjs::disable("disallowed_mismatch_pos") 
    } else {
        shinyjs::enable("disallowed_mismatch_pos") 
        updateSliderInput(session, "disallowed_mismatch_pos", value = 1)
    }
})
mismatchPreventObserver2 <- observeEvent(c(input$disallowed_mismatch_pos), {
    if (input$disallowed_mismatch_pos == 0) {
       updateRadioButtons(session, "allow_3prime_mismatch", selected="active")
    } else {
       updateRadioButtons(session, "allow_3prime_mismatch", selected="inactive")
    }
})
mismatchObserver1 <- observeEvent(c(input$are_mismatches_allowed), {
    if (input$are_mismatches_allowed != "active")  {
        # deactivate mismatch input and set to 0
        shinyjs::disable("allowed_mismatches") 
    } else {
        shinyjs::enable("allowed_mismatches") 
        updateSliderInput(session, "allowed_mismatches", value = 1)
    }
})
mismatchObserver2 <- observeEvent(input$allowed_mismatches, {
    if (input$allowed_mismatches == 0) {
       updateRadioButtons(session, "are_mismatches_allowed", selected="inactive")
    } else {
       updateRadioButtons(session, "are_mismatches_allowed", selected="active")
    }
})

allowed_nbr_of_mismatches <- reactive({
    return(ifelse(input$are_mismatches_allowed == "active", input$allowed_mismatches, 0))
})
stop.codons.allowed <- reactive({
    # Allow mismatches to create stop codons for coverage evaluation?
    return(ifelse(input$allowed_stop_codons == "active", TRUE, FALSE))
})
substitutions.allowed <- reactive({
    return(ifelse(input$allowed_substitutions == "active", TRUE, FALSE))
})

Na.concentration <- reactive({
    # Converts input Na salt concentration from mM to M
    validate(need(is.numeric(input$Na_concentration), "Concentration should be a numeric."))
    return(input$Na_concentration * 1e-3)
})
Mg.concentration <- reactive({
    # Converts the input Mg salt concentration from mM to M
    validate(need(is.numeric(input$Mg_concentration), "Concentration should be a numeric."))
    return(input$Mg_concentration * 1e-3)
})
K.concentration <- reactive({
    # Converts the input potassium concentration from mM to M
    validate(need(is.numeric(input$K_concentration), "Concentration should be a numeric."))
    return(input$K_concentration * 1e-3)
})
Tris.concentration <- reactive({ 
    # Converts the Tris buffer concentration concentration to its ion concentration, by halving the buffer concentration and converting from mM to M
    validate(need(is.numeric(input$Tris_concentration), "Concentration should be a numeric."))
    return((input$Tris_concentration/2) * 1e-3) # n.b.: tris buffer concentration is halved to receive [Tris+] concentration (assume Pk buffer)
})
use.taq.polymerase <- reactive({
    return(ifelse(input$use_taq_polymerase == "active", TRUE, FALSE))
})
annealing.temperature <- reactive({
    # The currently active annealing temperature (either automatically determined or input by the user)
    if (input$automatic_annealing_temp == "active") {
        # this should be computed 'at the right time' 
        #print("COMPUTING ANNEALING TEMPERATURE")
        # get Ta for currently active set:
        primer.data <- switch(input$set_meta_selector, 
            "all" = primer.data(), 
            "filtered" = current.filtered.primers(), 
            "optimized" = optimal.primers())
        template.data <- switch(input$set_meta_selector,
            "all" = rv_templates$cvg_all,
            "filtered" = rv_templates$cvg_filtered,
            "optimized" = rv_templates$cvg_optimized)
        annealing.temp <- try(min(openPrimeR:::compute_annealing_temp(primer.data, run.mode(), 
                                template.data, Na.concentration(), Mg.concentration(), K.concentration(), 
                                Tris.concentration(), primer.concentration()), na.rm=TRUE), silent = TRUE)
        #print("ANNEALING TEMPERATURE COMPUTED")
        if (class(annealing.temp) == "try-error") { # T_a can't be computed because coverage isnt' available yet
            Ta <- NULL # Ta will be computed by the backend when it's possible & necessary
        } else {
            Ta <- annealing.temp
        }
    } else {
        #print("USING INPUT ANNEALING TEMPERATURE")
        validate(need(is.numeric(input$annealing_temp), "Annealing temperature should be a numeric."))
        Ta <- input$annealing_temp
    }
    #message("Ta is: ", Ta)
    return(Ta)
})
primer.concentration <- reactive({
    # convert primer concentration from nM to M
    validate(need(is.numeric(input$primer_concentration), "Concentration should be a numeric."))
    return(input$primer_concentration * 1e-9)
})
template.concentration <- reactive({
    # convert template concentration from nM to M
    validate(need(is.numeric(input$template_concentration), "Concentration should be a numeric."))
    return(input$template_concentration * 1e-9)
})

constraintToolsObsever <- observeEvent(loaded.constraint.settings(), {
    # Disable constraints that are not supported by the available tools of the user
    # new function that can supplant MELTING software is available
    #if (!AVAILABLE.TOOLS()["MELTING"]) {
        #updateRadioButtons(session, "constraint_melting_temp_range", selected="inactive")
        #shinyjs::disable("constraint_melting_temp_range")
        #updateRadioButtons(session, "constraint_melting_temp_diff", selected="inactive")
        #shinyjs::disable("constraint_melting_temp_diff")
    #}
    if (!AVAILABLE.TOOLS()["ViennaRNA"]) {
        updateRadioButtons(session, "constraint_secondary_structure", selected="inactive", inline=TRUE)
        shinyjs::disable("constraint_secondary_structure")
        # disable optimization of template structures button
        shinyjs::disable("modify_binding_regions_secondary_structures")
    }
    if (!AVAILABLE.TOOLS()["OligoArrayAux"]) {
        # disable primer efficiency
        updateRadioButtons(session, "constraint_primer_efficiency", selected="inactive", inline=TRUE)
        shinyjs::disable("constraint_primer_efficiency")
        # disable self dimerization
        updateRadioButtons(session, "constraint_self_dimerization", selected="inactive")
        shinyjs::disable("constraint_self_dimerization")
        # disable cross dimerization
        updateRadioButtons(session, "constraint_cross_dimerization", selected="inactive")
        shinyjs::disable("constraint_cross_dimerization")
        # disable coverage model
        updateRadioButtons(session, "constraint_coverage_model", selected="inactive", inline=TRUE)
        shinyjs::disable("constraint_coverage_model")
        # disable annealing DeltaG
        updateRadioButtons(session, "constraint_annealing_DeltaG", selected="inactive", inline=TRUE)
        shinyjs::disable("constraint_annealing_DeltaG")
    }
    if (!AVAILABLE.TOOLS()["MAFFT"]) {
        # disable tree-init
        updateRadioButtons(session, "init_algo", selected = "naive")
        shinyjs::disable("init_algo")
        # disable conservation computation:
        shinyjs::disable("modify_binding_regions_conservation")
    }
})

# doesn't work consistenly ... therefore disabled
#constraintLimitObserver <- observeEvent(c(input$primer_analysis_type, loaded.constraint.settings()), {
    # disable constraint limit selection if we're not designing primers
    #   loaded.constraint.settings(): dependency to trigger the observer when primer_analysis_type isn't changed by the user
    # n.b.: sometimes shinyjs doesn't update correctly 
    # even though the observer is called.
#    print("CONSTRAINT LIMIT OBSERVER")
#    if (input$primer_analysis_type != "design") {
#        #print("not designing: disabling")
#        shinyjs::disable("limit_allowed_gc_clamp")
#        shinyjs::disable("limit_allowed_gc_ratio")
#        shinyjs::disable("limit_allowed_no_runs")
#        shinyjs::disable("limit_allowed_no_repeats")
#        shinyjs::disable("limit_allowed_primer_specificity")
#        shinyjs::disable("limit_allowed_melting_temp_range")
#        shinyjs::disable("limit_allowed_melting_temp_diff")
#        shinyjs::disable("limit_allowed_secondary_structure")
#        shinyjs::disable("limit_allowed_self_dimerization")
#        shinyjs::disable("limit_allowed_cross_dimerization")
#    } else {
#        #print("designing: enabling")
#        shinyjs::enable("limit_allowed_primer_coverage")
#        shinyjs::enable("limit_allowed_gc_clamp")
#        shinyjs::enable("limit_allowed_gc_ratio")
#        shinyjs::enable("limit_allowed_no_runs")
#        shinyjs::enable("limit_allowed_no_repeats")
#        shinyjs::enable("limit_allowed_primer_specificity")
#        shinyjs::enable("limit_allowed_melting_temp_range")
#        shinyjs::enable("limit_allowed_melting_temp_diff")
#        shinyjs::enable("limit_allowed_secondary_structure")
#        shinyjs::enable("limit_allowed_self_dimerization")
#        shinyjs::enable("limit_allowed_cross_dimerization")
#    }
#})

plotSelectorObserver <- observeEvent(current.settings(), {
    # decide which constraints can be plotted depending on loaded settings
    settings <- current.settings()
    # update constraint detail plot:
    choices <- names(openPrimeR::constraints(settings))
    names(choices) <- openPrimeR:::constraints_to_unit(names(openPrimeR::constraints(settings)), FALSE)
    updateSelectizeInput(session, "selected_other_result",
        choices = choices)
    updateSelectizeInput(session, "selected_other_plot",
        choices = choices)
    # update cvg constraint detail plot:
    choices <- names(openPrimeR::cvg_constraints(settings))
    names(choices) <- openPrimeR:::constraints_to_unit(names(openPrimeR::cvg_constraints(settings)), FALSE)
    updateSelectizeInput(session, "selected_cvg_constraints",
        choices = choices)
    updateSelectizeInput(session, "selected_cvg_comp_constraints",
        choices = choices)



})
AVAILABLE.TOOLS <- reactive({
    return(openPrimeR:::check.tool.function(frontend = TRUE))
})
annealingObserver <- observeEvent(input$constraint_annealing_DeltaG, {
    # Enables/disables constraint sliders
    if (input$constraint_annealing_DeltaG == "inactive") {
        shinyjs::disable("allowed_annealing_DeltaG")
    } else {
        if (AVAILABLE.TOOLS()["OligoArrayAux"]) {
            shinyjs::enable("allowed_annealing_DeltaG")
        } else {
            updateRadioButtons(session, "constraint_annealing_DeltaG", selected = "inactive")
            shinyjs::disable("allowed_annealing_DeltaG")
            shinyjs::disable("constraint_annealing_DeltaG")
        }
    }
})

constraintEfficiencyObserver <- observeEvent(input$constraint_primer_efficiency, {
    # Enables/disables constraint sliders
    if (input$constraint_primer_efficiency == "inactive") {
        shinyjs::disable("allowed_primer_efficiency")
    } else {
        if (AVAILABLE.TOOLS()["OligoArrayAux"]) {
            shinyjs::enable("allowed_primer_efficiency")
        } else {
            updateRadioButtons(session, "constraint_primer_efficiency", selected = "inactive")
            shinyjs::disable("allowed_primer_efficiency")
            shinyjs::disable("constraint_primer_efficiency")
        }
    }
})
constraintCoverageModelObserver <- observeEvent(input$constraint_coverage_model, {
    # Enables/disables constraint sliders
    if (input$constraint_coverage_model == "inactive") {
        shinyjs::disable("allowed_coverage_model")
    } else {
        if (AVAILABLE.TOOLS()["OligoArrayAux"]) {
            shinyjs::enable("allowed_coverage_model")
        } else {
            updateRadioButtons(session, "constraint_coverage_model", selected = "inactive")
            shinyjs::disable("allowed_coverage_model")
            shinyjs::disable("constraint_coverage_model")
        }
    }
})
constraintCoverageObserver <- observeEvent(input$constraint_primer_coverage, {
    # Enables/disables constraint sliders

    if (input$constraint_primer_coverage == "inactive") {
        shinyjs::disable("allowed_primer_coverage")
    } else {
        shinyjs::enable("allowed_primer_coverage")
    }
})
constraintLengthObserver <- observeEvent(input$constraint_primer_length, {
    # Enables/disables constraint sliders

    if (input$constraint_primer_length == "inactive") {
        shinyjs::disable("allowed_primer_length")
    } else {
        shinyjs::enable("allowed_primer_length")
    }
})
constraintGC_ClampObserver <- observeEvent(input$constraint_gc_clamp, {
    # Enables/disables constraint sliders

    if (input$constraint_gc_clamp == "inactive") {
        shinyjs::disable("allowed_gc_clamp")
        shinyjs::disable("limit_allowed_gc_clamp")
    } else {
        shinyjs::enable("allowed_gc_clamp")
        if (input$primer_analysis_type == "design") {
            shinyjs::enable("limit_allowed_gc_clamp")
        }
    }
})
constraintGC_ratio_Observer <- observeEvent(input$constraint_gc_ratio, {
    # Enables/disables constraint sliders

    if (input$constraint_gc_ratio == "inactive") {
        shinyjs::disable("allowed_gc_ratio")
        shinyjs::disable("limit_allowed_gc_ratio")
    } else {
        shinyjs::enable("allowed_gc_ratio")
        if (input$primer_analysis_type == "design") {
            shinyjs::enable("limit_allowed_gc_ratio")
        }
    }
})
constraintNoRunsObserver <- observeEvent(input$constraint_no_runs, {
    # Enables/disables constraint sliders

    if (input$constraint_no_runs == "inactive") {
        shinyjs::disable("allowed_no_runs")
        shinyjs::disable("limit_allowed_no_runs")
    } else {
        shinyjs::enable("allowed_no_runs")
        if (input$primer_analysis_type == "design") {
            shinyjs::enable("limit_allowed_no_runs")
        }
    }
})

constraintNoRepeatsObserver <- observeEvent(input$constraint_no_repeats, {
    # Enables/disables constraint sliders

    if (input$constraint_no_repeats == "inactive") {
        shinyjs::disable("allowed_no_repeats")
        shinyjs::disable("limit_allowed_no_repeats")
    } else {
        shinyjs::enable("allowed_no_repeats")
        if (input$primer_analysis_type == "design") {
            shinyjs::enable("limit_allowed_no_repeats")
        }
    }
})
constraintPrimerSpecObserver <- observeEvent(input$constraint_primer_specificity, {
    # Enables/disables constraint sliders

   if (input$constraint_primer_specificity == "inactive") {
        shinyjs::disable("allowed_primer_specificity")
        shinyjs::disable("limit_allowed_primer_specificity")
    } else {
        shinyjs::enable("allowed_primer_specificity")
        if (input$primer_analysis_type == "design") {
            shinyjs::enable("limit_allowed_primer_specificity")
        }
    }
})
constraintMeltingTempDiffObserver <- observeEvent(input$constraint_melting_temp_diff, {
    # Enables/disables constraint sliders
    if (input$constraint_melting_temp_diff == "inactive") {
        shinyjs::disable("allowed_melting_temp_diff")
        shinyjs::disable("limit_allowed_melting_temp_diff")
    } else {
        #if (AVAILABLE.TOOLS()["MELTING"]) { # only enable when tool is available
        shinyjs::enable("allowed_melting_temp_diff")
        if (input$primer_analysis_type == "design") {
            shinyjs::enable("limit_allowed_melting_temp_diff")
        }
        #} else {
        #    updateRadioButtons(session, "constraint_melting_temp_diff", selected = "inactive")
        #    shinyjs::disable("allowed_melting_temp_diff")
        #    shinyjs::disable("limit_allowed_melting_temp_diff")
        #    shinyjs::disable("constraint_melting_temp_diff")
        #}
    }
})
constraintMeltingTempRangeObserver <- observeEvent(input$constraint_melting_temp_range, {
    # Enables/disables constraint sliders

    if (input$constraint_melting_temp_range == "inactive") {
        shinyjs::disable("allowed_melting_temp_range")
        shinyjs::disable("limit_allowed_melting_temp_range")
    } else {
        #if (AVAILABLE.TOOLS()["MELTING"]) { # only enable when tool is available
        shinyjs::enable("allowed_melting_temp_range")
        if (input$primer_analysis_type == "design") {
            shinyjs::enable("limit_allowed_melting_temp_range")
        }
        #} else {
            #updateRadioButtons(session, "constraint_melting_temp_range", selected = "inactive")
            #shinyjs::disable("allowed_melting_temp_range")
            #shinyjs::disable("limit_allowed_melting_temp_range")
            #shinyjs::disable("constraint_melting_temp_range")
        #}
    }
})
constraintSecondaryStructObserver <- observeEvent(input$constraint_secondary_structure, {
    # Enables/disables constraint sliders
    if (input$constraint_secondary_structure == "inactive") {
        shinyjs::disable("allowed_secondary_structure")
        shinyjs::disable("limit_allowed_secondary_structure")
    } else {
        if (AVAILABLE.TOOLS()["ViennaRNA"]) {
            shinyjs::enable("allowed_secondary_structure")
            if (input$primer_analysis_type == "design") {
                shinyjs::enable("limit_allowed_secondary_structure")
            }
        } else {
            updateRadioButtons(session, "constraint_secondary_structure", selected = "inactive")
            shinyjs::disable("allowed_secondary_structure")
            shinyjs::disable("limit_allowed_secondary_structure")
            shinyjs::disable("constraint_secondary_structure")
        }
    }
})
constraintSelfDimerObserver <- observeEvent(input$constraint_self_dimerization, {
    # Enables/disables constraint sliders
    if (input$constraint_self_dimerization == "inactive") {
        shinyjs::disable("allowed_self_dimerization")
        shinyjs::disable("limit_allowed_self_dimerization")
    } else {
        if (AVAILABLE.TOOLS()["OligoArrayAux"]) {
            shinyjs::enable("allowed_self_dimerization")
            if (input$primer_analysis_type == "design") {
                shinyjs::enable("limit_allowed_self_dimerization")
            }
        } else {
            shinyjs::disable("allowed_self_dimerization")
            shinyjs::disable("limit_allowed_self_dimerization")
            # disable button itself to be sure
            updateRadioButtons(session, "constraint_self_dimerization", selected = "inactive")
            shinyjs::disable("constraint_self_dimerization")
        }
    }
})
constraintCrossDimerObserver <- observeEvent(input$constraint_cross_dimerization, {
    # Enables/disables constraint sliders

    #message("cross dimerization was changed to: ")
    #message(input$constraint_cross_dimerization)
    if (input$constraint_cross_dimerization == "inactive") {
        shinyjs::disable("allowed_cross_dimerization")
        shinyjs::disable("limit_allowed_cross_dimerization")
    } else {
        if (AVAILABLE.TOOLS()["OligoArrayAux"]) {
            shinyjs::enable("allowed_cross_dimerization")
            if (input$primer_analysis_type == "design") {
                shinyjs::enable("limit_allowed_cross_dimerization")
            }
        } else {
            updateRadioButtons(session, "constraint_cross_dimerization", selected = "inactive")
            shinyjs::disable("allowed_cross_dimerization")
            shinyjs::disable("limit_allowed_cross_dimerization")
            shinyjs::disable("constraint_cross_dimerization")
        }
    }
})
constraint_settings_cvg <- reactive({
    eff.setting <- NULL
    if (input$constraint_primer_efficiency == "active") {
        # primer efficiency is active
        eff.setting <- input$allowed_primer_efficiency
        names(eff.setting) <- c("min", "max")
    }
    term.mm.setting <- NULL
    if (input$allow_3prime_mismatch == "inactive") {
        # 3 prime mismatches are forbidden (inactive here)
        term.mm.setting <- input$disallowed_mismatch_pos
        names(term.mm.setting) <- "min"
    }
    annealing.deltaG.setting <- NULL
    if (input$constraint_annealing_DeltaG == "active") {
        # 3 prime mismatches are forbidden
        annealing.deltaG.setting <- input$allowed_annealing_DeltaG
        names(annealing.deltaG.setting) <- "max"
    }
    cvg.model.setting <- NULL
    if (input$constraint_coverage_model == "active") {
        # FPR limit for model:
        cvg.model.setting <- input$allowed_coverage_model 
        names(cvg.model.setting) <- "max"
    }
    settings <- openPrimeR:::get.cvg.constraint.settings(input$allowed_stop_codons == "active",
                                eff.setting, term.mm.setting, annealing.deltaG.setting, input$allowed_substitutions == "active", cvg.model.setting)
    return(settings)
})
constraint_settings_other <- reactive({ 
    # xml output for constraints that are not part of the 'filtering procedure'
    settings <- openPrimeR:::get.other.constraint.settings(allowed_nbr_of_mismatches(), 
        input$allowed_other_binding_ratio,input$allowed_binding_region_definition)
    return(settings)
})

analysisTypeObserver <- observeEvent(input$primer_analysis_type, {
    # ensures that primer coverage & length constraint are active for primer design

    if (input$primer_analysis_type == "design") {
        # deactivate button control: we need a couple of constraints
        updateRadioButtons(session, "constraint_primer_coverage", selected = "active")
        shinyjs::disable("constraint_primer_coverage")
        updateRadioButtons(session, "constraint_primer_length", selected = "active")
        shinyjs::disable("constraint_primer_length")
        # disable other binding ratio and set slider to 0
        updateSliderInput(session, "allowed_other_binding_ratio", value = 0)
        shinyjs::disable("allowed_other_binding_ratio") 

    } else {
        # enable other binding ratio
        shinyjs::enable("allowed_other_binding_ratio") 
        updateSliderInput(session, "allowed_other_binding_ratio", value = 1)
        # enable constraints that were deactivated before
        shinyjs::enable("constraint_primer_coverage")
        shinyjs::enable("constraint_primer_length")
        shinyjs::enable("allowed_other_binding_ratio") 
    }
    # change bg color depending on class
    modes <- c("evaluate", "design", "compare")
    for (mode in modes) {
        bg.class <- paste0("bg_", mode)
        shinyjs::toggleClass("headerPanel", class = bg.class, input$primer_analysis_type == mode)
        #shinyjs::addClass("headerPanel", class = "grad")

    }
})

ConstraintFileObserver <- observeEvent(input$load_constraints, { 
    # update current input settings file on user upload of settings
    rv_cur.input.data$settings <- input$load_constraints
})
SettingsChoiceObserver <- observeEvent(input$primer_analysis_type, {  
    # change available settings when modifying the analysis mode
    setting.options <- get.available.settings.view(
                        system.file("extdata", "settings", 
                                    package = "openPrimeR"), 
                        #input$use_taq_polymerase == "active"
                        taq.PCR = NULL, analysis.mode = input$primer_analysis_type)
    ## update selection:
    updateSelectizeInput(session, "load_available_constraints", 
                        choices = setting.options)
})

ConstraintAvailableObserver <- observeEvent(input$load_settings_button, {
    # converts the selected provided constraint settings by the app into a file choice
    if (length(input$load_settings_button) == 0) { # initial settings load:
        selection <- isolate(input$load_available_constraints)
    } else {
        selection <- input$load_available_constraints
    }
    if (length(selection) == 0 || selection == "") {
        return()
    }
    app.settings.folder <- system.file("extdata", "settings", 
                        package = "openPrimeR") 
    available.settings <- get.available.settings(app.settings.folder, 
                                                #input$use_taq_polymerase == "active"
                                                taq.PCR = NULL, analysis.mode = input$primer_analysis_type)
    path <- available.settings[grep(selection, available.settings)]
    out <- list("datapath" = path, "name" = selection)
    rv_cur.input.data$settings <- out
}, ignoreNULL = FALSE)

current.settings <- reactive({
    # use the backend units here:
    annealing.temp <- NULL
    if (input$automatic_annealing_temp != "active") {
        # set the input annealing temperature if provided
        annealing.temp <- annealing.temperature()
    }
    PCR.settings <- openPrimeR:::get.PCR.settings(use.taq.polymerase(), annealing.temp, Na.concentration(), 
                                                  Mg.concentration(), K.concentration(), 
                                                  Tris.concentration(), primer.concentration(), 
                                                  template.concentration(), input$cycles)
    # construction of current settings: uses the molar concentration for PCR concentrations i.e. the backend value. conversion of units takes place once again when we output to xml.
    #print("current pcr conditions:")
    #print(PCR.settings)
    # convert 'active' to TRUE/FALSE
    other.constraint.settings <- constraint_settings_other()
    cvg.constraint.settings <- constraint_settings_cvg()
    # get constraints:
    con.setting <- constraints()$active_settings
    constraint.limits <- constraint.limits()
    # input is ok but initialize screws up
    #print(con.setting)
    #print(constraint.limits)
    #print(cvg.constraint.settings)
    #print(PCR.settings)
    #print(other.constraint.settings)
    settings <- openPrimeR:::DesignSettings(con.setting,
        constraint.limits, cvg.constraint.settings,
        PCR.settings, other.constraint.settings)
    return(settings)
})
loaded.constraint.settings <- reactive({
    # loaded constraint settings from provided settings or file upload

    if (is.null(rv_cur.input.data$settings)) {
        return(NULL)
    }
    # read settings with `frontend = TRUE` such that we keep the original xml input data (PCR units are not changed). for the backend calls we use always current.settings() where we have the appropriate units from the current slider settings. 
    withProgress(message = 'Loading settings XML file ...', value = 0, {
        con.data <- withWarnings(openPrimeR:::read_settings(rv_cur.input.data$settings$datapath, frontend = TRUE))
    })
    # error handling
    for (i in seq_along(con.data$errors)) {
        error <- con.data$errors[[i]]
        print(error)
        if (inherits(error, "XML_Parsing_Error")) {
            toggleModal(session, "XML_Parsing_Error")
        } else {
            # unknown parsing error
            toggleModal(session, "XML_Parsing_Error")
        }
    }
    if (length(con.data$errors) != 0) {
        con.data <- NULL
    } else {
        con.data <- con.data$value
    }
    # activate confirm settings button when constraints are available
    if (length(con.data) == 0) {
        shinyjs::disable("confirm_settings_choice")
    } else {
        shinyjs::enable("confirm_settings_choice")
    }
    validate(need(con.data, "Could not read constraint data. Please check your input!"))
    isolate({
        if (length(current.settings()) != 0 && openPrimeR::PCR(current.settings())$cycles != 100) {
            # this is a bit hacky: assume that if the maximum of 100 cycles are used, the default settings are loaded -> don't switch to settings view initially
            updateTabsetPanel(session, "main", selected = "settings_view") # update view to settings tab
        }
    })
    return(con.data)
})

constraintsFromXMLObserver <- observeEvent(c(rv_cur.input.data$settings, input$reset_constraints), { 
    # loads constraints from input xml file when settings are uploaded  -> modify the UI elements accordingly
    # or should be restored (input$reset_constraints)
    if (is.null(loaded.constraint.settings())) {
      return(NULL)
    }
    withProgress(message = 'Implementing settings ...', value = 0, {
    # store all constraints here for deactivation/activation:
    UI.CONSTRAINT.MAPPING.FILTERS <- list(primer_coverage = "allowed_primer_coverage", 
        primer_length = "allowed_primer_length", 
        gc_clamp = "allowed_gc_clamp", gc_ratio = "allowed_gc_ratio", no_runs = "allowed_no_runs", 
        no_repeats = "allowed_no_repeats", melting_temp_range = "allowed_melting_temp_range", 
        primer_specificity = "allowed_primer_specificity", 
        self_dimerization = "allowed_self_dimerization", cross_dimerization = "allowed_cross_dimerization", 
        secondary_structure = "allowed_secondary_structure",
        melting_temp_diff = "allowed_melting_temp_diff")

    con.data <- loaded.constraint.settings()
    con.f <- openPrimeR:::constraints(con.data)
    for (i in seq_along(con.f)) {
        id <- paste("allowed_", names(con.f)[i], sep = "") # sliders should always be called allowed_<con_name> 
        data <- con.f[[i]]
        # update radio button: on/off?
        active.id <- paste("constraint_", names(con.f)[i], sep = "")
        updateRadioButtons(session, active.id, selected = "active")
        updateSliderInput(session, id,  value = unname(data))
    }
    # deactivate all constraints that were not selected in the xml
    possible.filters <- names(UI.CONSTRAINT.MAPPING.FILTERS)
    inactive.filters <- setdiff(possible.filters, names(con.f))
    #message(inactive.filters)
    for (i in seq_along(inactive.filters)) {
        inactive.id <- paste("constraint_", inactive.filters[i], sep = "")
        updateRadioButtons(session, inactive.id, selected = "inactive")
    }
    # update boundaries 
    con.f <- openPrimeR:::constraintLimits(con.data)
    for (i in seq_along(con.f)) {
        id <- names(con.f)[i]
        id <- paste("limit_allowed_", id, sep = "")
        data <- con.f[[i]]
        updateSliderInput(session, id,  value = unname(data))
    }
    #######################
    # Coverage constraints
    #######################
    cvg.conditions <- openPrimeR:::cvg_constraints(con.data)
    if ("primer_efficiency" %in% names(cvg.conditions)) {
        updateRadioButtons(session, "constraint_primer_efficiency", selected = "active")
        shinyjs::enable("allowed_primer_efficiency")
        updateSliderInput(session, "allowed_primer_efficiency", value = unname(cvg.conditions$primer_efficiency))
    } else {
        # allowed range can't be disabled upon startup, dont know why TODO
        updateRadioButtons(session, "constraint_primer_efficiency", selected = "inactive")
        shinyjs::disable("allowed_primer_efficiency")
    }
    if ("annealing_DeltaG" %in% names(cvg.conditions)) {
        updateRadioButtons(session, "constraint_annealing_DeltaG", selected = "active")
        shinyjs::enable("allowed_annealing_DeltaG")
        updateSliderInput(session, "allowed_annealing_DeltaG", value = unname(cvg.conditions$annealing_DeltaG["max"]))
    } else {
        shinyjs::disable("allowed_annealing_DeltaG")
        updateRadioButtons(session, "constraint_annealing_DeltaG", selected = "inactive")
    }
     if ("coverage_model" %in% names(cvg.conditions)) {
        updateRadioButtons(session, "constraint_coverage_model", selected = "active")
        shinyjs::enable("allowed_coverage_model")
        updateSliderInput(session, "allowed_coverage_model", value = unname(cvg.conditions$coverage_model["max"]))
    } else {
        shinyjs::disable("allowed_coverage_model")
        updateRadioButtons(session, "constraint_coverage_model", selected = "inactive")
    }

    if ("terminal_mismatch_pos" %in% names(cvg.conditions)) {
        # prevent mismatch binding
        updateRadioButtons(session, "allow_3prime_mismatch", selected = "inactive")
        shinyjs::enable("disallowed_mismatch_pos")
        updateSliderInput(session, "disallowed_mismatch_pos", value = cvg.conditions$terminal_mismatch_pos["min"] - 1)
    } else {
        # allow mismatch binding
        updateSliderInput(session, "disallowed_mismatch_pos", value = 0)
        shinyjs::disable("disallowed_mismatch_pos")
        updateRadioButtons(session, "allow_3prime_mismatch", selected = "active")
    }
    if ("stop_codon" %in% names(cvg.conditions) && all(cvg.conditions$stop_codon <= 0)) {
        # stop codon check
        updateRadioButtons(session, "allowed_stop_codons", selected = "inactive")
    } else {
        # no stop codon check
        updateRadioButtons(session, "allowed_stop_codons", selected = "active")
    }
    if ("substitution" %in% names(cvg.conditions) && all(cvg.conditions$substitution <= 0)) {
        # substitution check
        updateRadioButtons(session, "allowed_substitutions", selected = "inactive")
    } else {
        # stop codon check
        updateRadioButtons(session, "allowed_substitutions", selected = "active")
    }
    ##################
    # PCR conditions:
    ###################
    PCR.conditions <- openPrimeR:::PCR(con.data)
    for (i in seq_along(PCR.conditions)) {
        id <- names(PCR.conditions)[i] # constraint xml entries should correspond to slider IDs in UI
        data <- PCR.conditions[[i]]
        updateSliderInput(session, id,  value = unname(data))
    }
    # Constraint settings (e.g. binding conditions)
    constraint.settings <- openPrimeR:::conOptions(con.data)
    for (i in seq_along(constraint.settings)) {
        id <- names(constraint.settings)[i]
        if (id == "allowed_other_binding_ratio" && input$primer_analysis_type != "design") {
            # overwrite the input other binding ratio for eval/comparison mode
            #message("Enabling other binding ratio")
            shinyjs::enable("allowed_other_binding_ratio") 
            updateSliderInput(session, "allowed_other_binding_ratio", value = 1)
        } else {
            # all other options: simply activate according to input
            data <- constraint.settings[[i]]
            if (is.numeric(data)) {
                updateSliderInput(session, id,  value = unname(data))
            } else if (is.logical(data)) {
                # convert from boolean to shiny UI identifier for buttons
                updateRadioButtons(session, id,  selected = ifelse(data, "active", "inactive"))
            } else if (is.character(data)) {
                updateRadioButtons(session, id,  selected = unname(data))
            }
        }
    }
    }) # progress end
})

output$current_constraints <- DT::renderDataTable({
    # renders a data frame showing the current constraints
    tab <- switch(input$selected_settings_table,
                "constraints" = constraints.view(),
                "cvg_constraints" = cvg.constraints.view(),
                "opts" = constraint.options.view(),
                "PCR_options" = PCR.view())
    DT::datatable(
        tab,
        rownames = FALSE, escape = FALSE,
        caption="Overview of current settings.",
        options = list("dom" = "t", pageLength = 25)
    )
})
PCR.view <- reactive({
    validate(need(current.settings(), "No settings available yet."))
    table <- openPrimeR:::create.PCR.table(openPrimeR::PCR(current.settings()),
                            format.type = "shiny")
    return(table)
})
constraint.options.view <- reactive({
    validate(need(current.settings(), "No settings available yet."))
    table <- openPrimeR:::create.options.table(openPrimeR::conOptions(current.settings()), format.type = "shiny")
    return(table)
})
cvg.constraints.view <- reactive({
    validate(need(current.settings(), "No active constraints here."))
    cvg.table <- openPrimeR:::create.constraint.table(openPrimeR::cvg_constraints(current.settings()), 
                            format.type = "shiny")
    return(cvg.table)
})
constraints.view <- reactive({
    # the constraints table
    validate(need(constraints(), "No active constraints here."))
    used.constraints.fw <- NULL # relaxed constraints from opti
    used.constraints.rev <- NULL
    if (length(rv_primers$optimal_data) != 0) {
        if (!is.null(rv_primers$optimal_data$used_constraints[["fw"]])) {
            used.constraints.fw <- openPrimeR::constraints(rv_primers$optimal_data$used_constraints[["fw"]])
        }
        if (!is.null(rv_primers$optimal_data$used_constraints[["rev"]])) {
            used.constraints.rev <- openPrimeR::constraints(rv_primers$optimal_data$used_constraints[["rev"]])
       }
    }
    constraints <- constraints()$active_settings
    filter.table <- openPrimeR:::create.constraint.table(constraints, 
                        constraint.limits(), used.constraints.fw,
                        used.constraints.rev, format.type = "shiny")
    return(filter.table)
})

active.constraints <- reactive({
    # names of the currently active constraints
    con <- constraints()
    return(con[["active"]])
})

constraint.limits <- reactive({
    # list with all constraint limits for filtering

    # primer length: can't be relaxed -> set to allowed value
    constraint.values <- list(
        # no limit for primer coverage: not relaxed
        "primer_coverage" = c("min" =  input$allowed_primer_coverage[1]), 
        "gc_clamp" = c("min" = input$limit_allowed_gc_clamp[1], 
                       "max" = input$limit_allowed_gc_clamp[2]), 
        "gc_ratio" = c("min" = input$limit_allowed_gc_ratio[1], 
                       "max" = input$limit_allowed_gc_ratio[2]), 
        "no_runs" = c("min" = input$limit_allowed_no_runs[1], 
                      "max" = input$limit_allowed_no_runs[2]), 
        "no_repeats" = c("min" = input$limit_allowed_no_repeats[1], 
                         "max" = input$limit_allowed_no_repeats[2]), 
        "self_dimerization" = c("min" = input$limit_allowed_self_dimerization[1]),
        "cross_dimerization" = c("min" = input$limit_allowed_cross_dimerization[1]),
        "secondary_structure" = c("min" = input$limit_allowed_secondary_structure[1]),
        "melting_temp_range" = c("min" = input$limit_allowed_melting_temp_range[1], 
                                "max" = input$limit_allowed_melting_temp_range[2]), 
        "melting_temp_diff" = c("min" = input$limit_allowed_melting_temp_diff[1],
                                "max" = input$limit_allowed_melting_temp_diff[2]),
        # no limits for primer length: not relaxed
        "primer_length" = c("min" = input$allowed_primer_length[1], 
                            "max" = input$allowed_primer_length[2]), 
        "primer_specificity" = c("min" = input$limit_allowed_primer_specificity[1], 
                                 "max" = input$limit_allowed_primer_specificity[2]))
    return(constraint.values)
})

constraints <- reactive({
    # list with the constraint settings
    # info about all constraints: use radio buttons
    # order of constraint evaluations was determined here before
    # we implemented this in DesignSettings itself

    # primer coverage
    primer.coverage <- input$constraint_primer_coverage
    active.constraints <- NULL
    if (primer.coverage == "active") {
        active.constraints <- c(active.constraints, "primer_coverage") # need primer coverage first for relaxation ..
    }
    # primer length
    primer.length <- input$constraint_primer_length
    if (primer.length == "active") {
        active.constraints <- c(active.constraints, "primer_length")
    }
    # primer specificity
    primer.specificity <- input$constraint_primer_specificity
    if (primer.specificity == "active") {
        active.constraints <- c(active.constraints, "primer_specificity")
    }
    # gc clamp
    gc.clamp <- input$constraint_gc_clamp
    if (gc.clamp == "active") {
        active.constraints <- c(active.constraints, "gc_clamp")
    }
    # gc ratio
    gc.ratio <- input$constraint_gc_ratio
    if (gc.ratio == "active") {
        active.constraints <- c(active.constraints, "gc_ratio")
    }
    # runs: repetition of the same base
    no.of.runs <- input$constraint_no_runs
    if (no.of.runs == "active") {
        active.constraints <- c(active.constraints, "no_runs")
    }
    # repeat: repetition of a dinucleotide
    no.of.repeats <- input$constraint_no_repeats
    if (no.of.repeats == "active") {
        active.constraints <- c(active.constraints, "no_repeats")
    } 
    # self-dimerization
    self.complementary.end <- input$constraint_self_dimerization
    if (self.complementary.end == "active") {
        active.constraints <- c(active.constraints, "self_dimerization")
    }
    # melting temperature
    melting.temp <- input$constraint_melting_temp_range
    if (melting.temp == "active") {
        active.constraints <- c(active.constraints, "melting_temp_range")
    }
    melting.temp.diff <- input$constraint_melting_temp_diff
    if (melting.temp.diff == "active") {
        active.constraints <- c(active.constraints, "melting_temp_diff")
    }
    # secondary structures
    secondary.structure <- input$constraint_secondary_structure
    if (secondary.structure == "active") {
        active.constraints <- c(active.constraints, "secondary_structure")
    }
    # cross-dimerization
    cross.complementary.end <- input$constraint_cross_dimerization
    if (cross.complementary.end == "active") {
        active.constraints <- c(active.constraints, "cross_dimerization")
    }
    #message(active.constraints)
    #########
    constraint.values <- list(
        "primer_coverage" = c("min" = input$allowed_primer_coverage[1]), 
        "gc_clamp" = c("min" = input$allowed_gc_clamp[1], 
                       "max" = input$allowed_gc_clamp[2]), 
        "gc_ratio" = c("min" = input$allowed_gc_ratio[1], 
                       "max" = input$allowed_gc_ratio[2]), 
        "no_runs" = c("min" = input$allowed_no_runs[1], 
                      "max" = input$allowed_no_runs[2]), 
        "no_repeats" = c("min" = input$allowed_no_repeats[1], 
                         "max" = input$allowed_no_repeats[2]), 
        "self_dimerization" = c("min" = input$allowed_self_dimerization[1]),
        "cross_dimerization" = c("min" = input$allowed_cross_dimerization[1]),
        "secondary_structure" = c("min" = input$allowed_secondary_structure[1]),
        "melting_temp_range" = c("min" = input$allowed_melting_temp_range[1], 
                                 "max" = input$allowed_melting_temp_range[2]), 
        "melting_temp_diff" = c("min" = input$allowed_melting_temp_diff[1], 
                                "max" = input$allowed_melting_temp_diff[2]),
        "primer_length" = c("min" = input$allowed_primer_length[1], 
                            "max" = input$allowed_primer_length[2]), 
        "primer_specificity" = c("min" = input$allowed_primer_specificity[1],
                                 "max" = input$allowed_primer_specificity[2]))
    m <- match(active.constraints, names(constraint.values))
    active.settings <- constraint.values[m] # actually used settings atm
    result <- list("active" = active.constraints, "values" = constraint.values, "active_settings" = active.settings)
    return(result)
})


activateFilterObserver <- observeEvent(input$activate_all_filters, { 
    # activate all filters on click
    s <- "active"
    updateRadioButtons(session, "constraint_primer_coverage", selected=s, inline=TRUE)
    updateRadioButtons(session, "constraint_primer_length", selected=s, inline=TRUE)
    updateRadioButtons(session, "constraint_gc_clamp", selected=s, inline = TRUE)
    updateRadioButtons(session, "constraint_gc_ratio", selected=s, inline = TRUE)
    updateRadioButtons(session, "constraint_no_runs", selected=s, inline=TRUE)
    updateRadioButtons(session, "constraint_no_repeats", selected=s, inline=TRUE)
    updateRadioButtons(session, "constraint_melting_temp_range", selected=s, inline=TRUE)
    updateRadioButtons(session, "constraint_secondary_structure", selected=s, inline=TRUE)
    updateRadioButtons(session, "constraint_primer_specificity", selected=s, inline=TRUE)
    updateRadioButtons(session, "constraint_self_dimerization", selected=s, inline=TRUE)
    updateRadioButtons(session, "constraint_cross_dimerization", selected=s, inline=TRUE)
    updateRadioButtons(session, "constraint_melting_temp_diff", selected=s, inline=TRUE)
})

deactivateConstraintsObserver <- observeEvent(input$deactivate_all_filters, { 
    # deactivate all constraints on clicking the button
    s <- "inactive"
    if (input$primer_analysis_type != "design") { # never deactivate these when designing
        updateRadioButtons(session, "constraint_primer_coverage", selected=s, inline=TRUE)
        updateRadioButtons(session, "constraint_primer_length", selected=s, inline=TRUE)
    }
    updateRadioButtons(session, "constraint_gc_clamp", selected=s, inline = TRUE)
    updateRadioButtons(session, "constraint_gc_ratio", selected=s, inline = TRUE)
    updateRadioButtons(session, "constraint_no_runs", selected=s, inline=TRUE)
    updateRadioButtons(session, "constraint_no_repeats", selected=s, inline=TRUE)
    updateRadioButtons(session, "constraint_melting_temp_range", selected=s, inline=TRUE)
    updateRadioButtons(session, "constraint_secondary_structure", selected=s, inline=TRUE)
    updateRadioButtons(session, "constraint_primer_specificity", selected=s, inline=TRUE)
    updateRadioButtons(session, "constraint_self_dimerization", selected=s, inline=TRUE)
    updateRadioButtons(session, "constraint_cross_dimerization", selected=s, inline=TRUE)
    updateRadioButtons(session, "constraint_melting_temp_diff", selected=s, inline=TRUE)
})

ToolInfoObserver <- observeEvent(input$third_party_tools, {
    toggleModal(session, "MissingTools") 
})

output$ToolOverview <- DT::renderDataTable({
    tool.df <- openPrimeR:::build.tool.overview(AVAILABLE.TOOLS(), TRUE)
    DT::formatStyle(DT::datatable(tool.df, rownames=FALSE, escape = FALSE,
                    caption="Overview of the installation status of required third-party tools. If third-party software is missing, the indicated features will not be available."),
        "Status", backgroundColor = DT::styleEqual(c("Unavailable", "Available"), 
        c("#ff9999", "#99d6ff")))
})
StartObserver <- observe ({
    # run only on startup
    if (!AVAILABLE.TOOLS()["MAFFT"]) {
        # disable tree-init
        updateRadioButtons(session, "init_algo", selected = "naive")
        shinyjs::disable("init_algo")
    }
    if (!all(AVAILABLE.TOOLS()) ) {
        # inform user about missing tools only
        toggleModal(session, "MissingTools") 
        updateActionButton(session, "third_party_tools", icon = icon("exclamation-triangle"))

   }
})
output$CoverageBox <- renderUI({
    # provides a box summarizing the current coverage settings
    cvg.constraints <- openPrimeR::cvg_constraints(current.settings())
    con.options <- openPrimeR:::rename.constraint.options(openPrimeR::conOptions(current.settings()))
    basic.format <- paste0("<b>Basic Coverage</b><ul style = 'list-style-type:none; padding-left:10px;'>",
                    paste0(sapply(seq_along(con.options), function(x) 
                        paste0("<li><i>", names(con.options)[x], "</i>", ": ", unname(con.options)[[x]], "</li>")),
                        collapse = ""),
                "</ul>")
    if (length(cvg.constraints) != 0) { # at least one coverage constraint is present
        # show a listing of active constraints
        con.strings <- sapply(seq_along(cvg.constraints), function(x) {
            values <- cvg.constraints[[x]]
            if (length(values) == 1 && names(values) == "min") {
                out <- paste0("&ge; ", unname(values))
            } else if (length(values) == 1 && names(values) == "max") {
                out <- paste0("&le; ", unname(values))
            } else if (length(values) == 2) { # length 2
                out <- paste0("[", paste0(unname(cvg.constraints[[x]]), collapse = ","), "]")
            } else {
                out <- ""
            }
        })
        extended.format <- paste0("<b>Coverage</b><ul style = 'list-style-type:none; padding-left:10px;'>",
                        paste0(sapply(seq_along(cvg.constraints), function(x) 
                            paste0("<li><i>", 
                                openPrimeR:::constraints_to_unit(names(cvg.constraints)[x], FALSE, "HTML"),
                                "</i>", ": ", con.strings[x], "</li>"
                            )
                        ),
                        collapse = ""), "</ul>")
    } else {
        # show a warning that no coverage constraint was active
        extended.format <- paste0("<b>Coverage</b><ul style = 'list-style-type:none; padding-left:10px;'>",
                                "<li style = 'color:red'>",
                                "Warning: Estimated coverage may be inaccurate since no coverage constraint was active.",
                                "</li>",
                                "</ul>")
    }
    out <- HTML(paste0(h4("Coverage Conditions"), br(),
                        basic.format, extended.format))
    return(out)
})
