####### Server UI navigation functionalities rv_navigation.data states: vector: can we
####### go to the next step?, page: current page in the sequence
rv_navigation.data <- reactiveValues(states = NULL, page = 1)


loadedConstraintsObserver <- observeEvent(loaded.constraint.settings(), {
    # enable/disable constraint settings panels depending on settings availability
    if (length(loaded.constraint.settings()) == 0) {
        # disable binding/other conditions using jqery
        shinyjs::disable(selector = "#input_settings_collapse div[value=binding_conditions_panel]")
        shinyjs::disable(selector = "#input_settings_collapse div[value=customize_settings_panel]")
        shinyjs::disable(selector = "#input_settings_collapse div[value=customize_general_settings_panel]")
    } else {
        # settings loaded -> enable constraint panels
        shinyjs::enable(selector = "#input_settings_collapse div[value=binding_conditions_panel]")
        shinyjs::enable(selector = "#input_settings_collapse div[value=customize_settings_panel]")
        shinyjs::enable(selector = "#input_settings_collapse div[value=customize_general_settings_panel]")
    }
})
loadedTemplatesObserver <- observeEvent(seq.data(), {
    # enable/disable panels after template data has been uploaded
    if (length(seq.data()) == 0) {
        shinyjs::disable(selector = "#customize_allowed_regions_collapse div[value=customize_allowed_regions_panel]")
        shinyjs::disable(selector = "#other_template_data_collapse div[value=allowed_template_panel]")
    } else {
        shinyjs::enable(selector = "#customize_allowed_regions_collapse div[value=customize_allowed_regions_panel]")
        shinyjs::enable(selector = "#other_template_data_collapse div[value=allowed_template_panel]")
    }
    
})
######## EXIT
exitPanelObserver <- observeEvent(input$quit_tool, {
    # create the exit bsmodal when user presses the exit button
    session$sendCustomMessage(type = "jsCode", list(value = "$('#ExitInfo').modal('show')"))
})

exitButtonObserver <- observeEvent(input$exitButton, {
    session$sendCustomMessage(type = "jsCode", list(value = "$('#ExitInfo').modal('hide')"))  # hide exit dialog
    session$sendCustomMessage(type = "jsCode", list(value = "$('#ExitScreen').modal('show')"))  # show exit screen
    stopApp()
})

######## RESET
resetPanelObserver <- observeEvent(input$reset_tool, {
    # create the reset bsmodal when user pressets the reset tab
    session$sendCustomMessage(type = "jsCode", list(value = "$('#ResetInfo').modal('show')"))
})

observeEvent(input$reset_button, {
    # reset the shiny session
    shinyjs::js$reset()
})

######## BASIC NAVIGATION
cur.event.sequence <- reactive({
    # sets the current event sequence on change to the analysis type and resets
    # rv_navigation.data$stages according to the selection
    sequence <- NULL
    default.event.sequence <- c("templates", "primers", "settings", "analyze")
    # all analysis modes have the same event sequence:
    if (input$primer_analysis_type == "evaluate") {
        sequence <- default.event.sequence
    } else if (input$primer_analysis_type == "design") {
        sequence <- default.event.sequence
    } else if (input$primer_analysis_type == "compare") {
        sequence <- default.event.sequence
    }
    # reset navigation process when we start a new primer analysis type:
    rv_navigation.data$stages <- rep(FALSE, length(sequence))  # reset navigation progress
    rv_navigation.data$page <- 1 # jump to the template page
    update.following.navigation(session, "templates", "disable") # disable all navigation tabs following 'templates'
    return(sequence)
})
pageObserver <- observeEvent(input$settingsPanel, {
    # update selected page by user according to the current selected panel (important
    # in case users doesn't use the navigation buttons but the panels directly)
    if (input$settingsPanel == "input_data_panel") {
        rv_navigation.data$page <- 1
    } else if (input$settingsPanel == "primer_input_tab") {
        rv_navigation.data$page <- 2
    } else if (input$settingsPanel == "constraint_panel") {
        rv_navigation.data$page <- 3
    } else if (input$settingsPanel == "analyze_panel") {
        rv_navigation.data$page <- 4
    }
    # message(paste('current navigation page is: ', rv_navigation.data$page, sep = ''))
})
observeEvent(input$prevBtn, {
    # move back one page when prevBtn was pressed
    rv_navigation.data$page <- rv_navigation.data$page - 1  # logic for phases
})

observeEvent(input$nextBtn, {
    # advance one page when nextBtn was pressed
    rv_navigation.data$page <- rv_navigation.data$page + 1  # logic for phases
})

navigationObserver <- observeEvent(rv_navigation.data$page, {
    # update currently active tab in data input when the page changes (next/prev
    # button) ui response
    if (length(cur.event.sequence()) == 0) {
        return()
    }
    cur.phase <- cur.event.sequence()[rv_navigation.data$page]
    if (cur.phase == "templates") {
        updateTabsetPanel(session, "settingsPanel", selected = "input_data_panel")
    } else if (cur.phase == "primers") {
        updateTabsetPanel(session, "settingsPanel", selected = "primer_input_tab")
    } else if (cur.phase == "settings") {
        updateTabsetPanel(session, "settingsPanel", selected = "constraint_panel")
    } else if (cur.phase == "analyze") {
        updateTabsetPanel(session, "settingsPanel", selected = "analyze_panel")
    } else {
        warning(paste("Non-supported navigation phase: ", cur.phase, sep = ""))
    }
})
NavigationStateObserver <- observe({
    # determine if we have fulfilled the criterion for a certain phase to determine
    # whether we can use the next button or not
    if (length(cur.event.sequence()) == 0) {
        return()
    }
    cur.phase <- cur.event.sequence()[rv_navigation.data$page]
    if (tail(cur.event.sequence(), n = 1) == cur.phase) {
        # we're already in the last phase -> never activate the 'next' button
        return()
    }
    if (cur.phase == "templates") {
        # user is in the template input phase check whether templates have been uploaded
        # seq data or rv_comparison.data uploaded
        if (input$primer_analysis_type == "compare") {
            # template comparison
            pass <- length(rv_comparison.data$seqs) != 0
        } else  { # for design/evaluate
            pass <- length(seq.data()) != 0 && !is.na(seq.data())
        }
        rv_navigation.data$state[rv_navigation.data$page] <- pass  # next button should activate/disable now
        action <- ifelse(pass, "enable", "disable")
        update.following.navigation(session, cur.phase, action)
        # message('passed the template phase!')
        # deactivate look-ahead (template choice can invalidate primers)
        if ((input$primer_analysis_type == "compare" && length(rv_comparison.data$primer_fnames) == 0) ||
            (input$primer_analysis_type != "compare" && length(primer.data()) == 0)) {
                update.following.navigation(session, "primers", "disable")
        }
    } else if (cur.phase == "primers") {
        # user should enter primers now check whether primers are available for design:
        # no primer input necessary -> activate next button directly
        if (input$primer_analysis_type == "design") {
            # no required input for primer design :-)
            pass <- TRUE
            # message('passed the primer phase!')
        } else if (input$primer_analysis_type == "compare") { 
            #  comparison:
            pass <- length(rv_comparison.data$primer_fnames) != 0
        } else { # evaluate
            pass <- length(primer.data()) != 0
        }
        rv_navigation.data$state[rv_navigation.data$page] <- pass
        action <- ifelse(pass, "enable", "disable")
        update.following.navigation(session, cur.phase, action)
    } else if (cur.phase == "settings") {
        # update settings selector with checkmark
        pass <- length(loaded.constraint.settings()) != 0
        action <- ifelse(pass, "enable", "disable")
        rv_navigation.data$state[rv_navigation.data$page] <- pass  # next button should activate now
        update.following.navigation(session, cur.phase, action)
    } else {
        message(paste("Unknown current phase: ", cur.phase, sep = ""))
    }
})
navigationBarObserver <- observe({
    show.panels <- c("input_data_panel", "primer_input_tab", "constraint_panel", 
        "analyze_panel")  # show navigation buttons only on these pages
    if (input$settingsPanel %in% show.panels) {
        # display navigation buttons message('we are in a show panel for navigation ..')
        # last page
        if (rv_navigation.data$page >= length(rv_navigation.data$stages)) {
            shinyjs::hide("nextBtn")
        } else {
            shinyjs::show("nextBtn")
        }
        if (rv_navigation.data$page == 1) {
            # no prev button on first page
            shinyjs::hide("prevBtn")
        } else {
            shinyjs::show("prevBtn")
        }
    } else {
        # hide navigation buttons
        shinyjs::hide("nextBtn")
        shinyjs::hide("prevBtn")
    }
    shinyjs::toggleState(id = "prevBtn", condition = rv_navigation.data$page > 1)  # we can always go back :-)
    shinyjs::toggleState(id = "nextBtn", condition = rv_navigation.data$state[rv_navigation.data$page])  # we can only go to the next page if we are allowed to in the current stored state
})

######### STYLES: highlight currently active bscollapse

DesignPrimersActionsObserver <- observeEvent(input$design_options_algorithms_collapse, 
    {
        # Update style of open design bsPanel
        panels <- c("design_algo_panel", "design_algo_options_panel")
        styles <- rep("default", length(panels))
        names(styles) <- panels
        if (length(input$design_options_algorithms_collapse) != 0) {
            styles[input$design_options_algorithms_collapse] <- "primary"  # set open tab to primary
        }
        updateCollapse(session, "design_options_algorithms_collapse", style = styles)
})
DesignStyleObserver <- observeEvent(input$design_options_algorithms_collapse, {
    # Update the style of primer upload panel according to active panel
    panels <- c("design_algo_panel", "design_optimize_binding")
    styles <- rep("default", length(panels))
    names(styles) <- panels
    if (length(input$design_options_algorithms_collapse) != 0) {
        # a panel is selected
        active.tab <- tail(input$design_options_algorithms_collapse, n = 1)
        styles[active.tab] <- "primary"  # set open tab to primary
    }
    updateCollapse(session, "design_options_algorithms_collapse", style = styles)
})

PrimerUploadStyleObserver <- observeEvent(input$primer_upload_collapse, {
    # Update the style of primer upload panel according to active panel
    panels <- c("primer_header_structure_panel", "primer_upload_panel")
    styles <- rep("default", length(panels))
    names(styles) <- panels
    if (length(input$primer_upload_collapse) != 0) {
        # a panel is selected
        active.tab <- tail(input$primer_upload_collapse, n = 1)
        styles[active.tab] <- "primary"  # set open tab to primary
    }
    updateCollapse(session, "primer_upload_collapse", style = styles)
})
CoverageConditionsStyleObserver <- observeEvent(input$coverage_conditions_collapse, {
    # Update the style of the coverage conditions 
    panels <- c("basic_cvg_conditions", "ext_cvg_conditions")
    styles <- rep("default", length(panels))
    names(styles) <- panels
    if (length(input$coverage_conditions_collapse) != 0) {
        # a panel is selected
        active.tab <- tail(input$coverage_conditions_collapse, n = 1)
        styles[active.tab] <- "primary"  # set open tab to primary
    }
    updateCollapse(session, "coverage_conditions_collapse", style = styles)
})
CoverageConditionsExtendedStyleObserver <- observeEvent(input$extended_cvg_collapse, {
    # Update the style of the extended coverage conditions 
    panels <- c("ext_cvg_conditions_binding", "ext_cvg_conditions_mismatches")
    styles <- rep("default", length(panels))
    names(styles) <- panels
    if (length(input$extended_cvg_collapse) != 0) {
        # a panel is selected
        active.tab <- tail(input$extended_cvg_collapse, n = 1)
        styles[active.tab] <- "primary"  # set open tab to primary
    }
    updateCollapse(session, "extended_cvg_collapse", style = styles)
})


SettingsCollapseStyleObserver <- observeEvent(input$input_settings_collapse, {
    # Update style of settings panel according to user selection
    settings.panels <- c("load_settings_panel", "binding_conditions_panel", "customize_settings_panel", 
        "customize_general_settings_panel")  # values of the collapsePanels
    styles <- rep("default", length(settings.panels))  # one style for every collapse panel
    names(styles) <- settings.panels
    if (length(input$input_settings_collapse) != 0) {
        active.tab <- tail(input$input_settings_collapse, n = 1)
        styles[active.tab] <- "primary"  # set open tab to primary
        updateCollapse(session, "input_settings_collapse", style = styles)
    }
})
TemplateCollapseStyleObserver <- observeEvent(input$template_collapse_analysis, 
    {
        # clean up open tabs
        tab.names <- c("template_input_panel", "allowed_template_panel")
        styles <- rep("default", length(tab.names))
        names(styles) <- tab.names
        if (length(input$template_collapse_analysis) != 0) {
            # select only the last chosen bscollapse panel
            # tail() is necessary since multiple may be selected when
            # input elements change (the last active one counts)
            active.tab <- tail(input$template_collapse_analysis, n = 1)
            styles[active.tab] <- "primary"  # set open tab to primary
            updateCollapse(session, "template_collapse_analysis", style = styles)
        } else {
            # non-selected -> all are set to default
            updateCollapse(session, "template_collapse_analysis", style = styles)
        }
        #print("Template input observer:")
        #print(input$personal_template_options)
        #print(input$template_collapse_analysis)
})
TemplatePersonalCollapseStyleObserver <- observeEvent(input$personal_template_options, {
    # change the style of personal template options (upload) collapse 
    tab.names <- c("basic_personal_template_options", "template_expert_panel")
    styles <- rep("default", length(tab.names))
    names(styles) <- tab.names
    if (length(input$personal_template_options) != 0) {
        # select only the last chosen bscollapse panel
        # tail() is necessary since multiple may be selected when
        # input elements change (the last active one counts)
        active.tab <- tail(input$personal_template_options, n = 1)
        styles[active.tab] <- "primary"  # set open tab to primary
        updateCollapse(session, "personal_template_options", style = styles)
    } else {
        # non-selected -> all are set to default
        updateCollapse(session, "personal_template_options", style = styles)
    }
})
observeEvent(input$IMGT_template_confirm_button, {
    # on confirmation of IMGT templates, move to binding region customization
        # allowed panel isn't active yet
        # trigger opening, but then remove the duplicated entry?
        updateCollapse(session, "template_collapse_analysis", close = c("template_input_panel"), 
            open = c("allowed_template_panel"))
})

session$onSessionEnded(function() {
    # actions to perform at the end of a shiny session: cleanu
    session$sendCustomMessage(type = "jsCode", list(value = "$('#ExitScreen').modal('show')"))  # show exit screen
    #cat("Ending the Shiny session ...")
    # reset to default ggplot theme
    ggplot2::theme_set(OLD.GG.THEME)
    #cat("\n\to Removing all observers ...\n")
    # suspend all observers:
    IMGT_EvaluatedObserver$suspend()
    CoreObserver$suspend()
    StartObserver$suspend()
    pageObserver$suspend()
    resetPanelObserver$suspend()
    loadedConstraintsObserver$suspend()
    loadedTemplatesObserver$suspend()
    ComparisonPrimerSuppliedObserver$suspend()
    ComparisonPrimerInputObserver$suspend()
    availablePrimerComparisonUpdater$suspend()
    comparisonTemplateOtherObserver$suspend()
    comparisonTemplateIMGTObserver$suspend()
    analysisTypeObserver$suspend()
    ConstraintFileObserver$suspend()
    ConstraintAvailableObserver$suspend()
    SettingsCollapseStyleObserver$suspend()
    NavigationStateObserver$suspend()
    InputPrimerObserver$suspend()
    activateFilterObserver$suspend()
    deactivateConstraintsObserver$suspend()
    comparisonResetObserver$suspend()
    comparisonFilterObserver$suspend()
    IMGT_TemplateDataObserver$suspend()
    templateInvalidationObserver$suspend()
    primerViewObserver$suspend()
    primerViewObserverGroup$suspend()
    PrimerTabObserver$suspend()
    SeqTabObserver$suspend()
    primerComparisonObserver$suspend()
    exitPanelObserver$suspend()
    exitButtonObserver$suspend()
    # constraint observers:
    constraintsFromXMLObserver$suspend()
    constraintCrossDimerObserver$suspend()
    constraintSelfDimerObserver$suspend()
    constraintSecondaryStructObserver$suspend()
    constraintMeltingTempDiffObserver$suspend()
    constraintMeltingTempRangeObserver$suspend()
    constraintPrimerSpecObserver$suspend()
    constraintNoRepeatsObserver$suspend()
    constraintNoRunsObserver$suspend()
    constraintGC_ratio_Observer$suspend()
    constraintGC_ClampObserver$suspend()
    constraintLengthObserver$suspend()
    constraintCoverageObserver$suspend()
    constraintEfficiencyObserver$suspend()
    #constraintLimitObserver$suspend()
    constraintToolsObsever$suspend()
    #constraintSecondaryStructOptiObserver$suspend()
    #constraintEfficiencyOptiObserver$suspend()
    #constraintMeltingTempOptiObserver$suspend()
    #comparisonConstraintObserver$suspend()
    comparisonFileObserver$suspend()
    comparisonTemplateObserver$suspend()
    notifyRelaxation$suspend()  # pop-up when relaxation occured
    notifyNotAllowedBinding$suspend()
    #message("o Cleaning temporary files ...")
    tmp.files <- list.files(tempdir(), full.names = TRUE)
    tmp.files <- tmp.files[!dir.exists(tmp.files)]  # only select files -> widgetbinding is required on restart
    sapply(tmp.files, function(x) unlink(x, recursive = TRUE))
    #cat("\to Killing remaining children ...\n")  # after closing browser!
    #children <- parallel:::children()  # any child processes still remaining
    #parallel:::mckill(children, KILL.METHOD.ALT)  # forcefully (because shiny is still running) kill all remaining children
    
    #cat("\to Triggering garbage collection ...\n")
    
    gc()
    #detach("package:shinyjs") # causes problems with print() otherwise
})

primerInvalidationObserver <- observeEvent(input.primers(), {
    # Invalidate previously computed rv_values when primers have changed
    reset.reactive.values(rv_values)
    reset.reactive.values(rv_templates, keep = c("SeqTab", "raw_seqs"))  # templates not valid anymore
    reset.reactive.values(rv_primers, keep = c("PrimerTab", "all"))  # retain reactive values that were just loaded on input, clean all others
    if ("primer_coverage" %in% colnames(input.primers())) {
        # special case for csv: don't reset evaluated constraints
        # update cvg in templates:
        template.df <- try(openPrimeR::update_template_cvg(current.seqs(), input.primers()))
        if (class(template.df) != "try-error") {
            # primers & templates seem to correspond to one another -> use input coverage values
            rv_templates$cvg_all <- template.df
            primers <- input.primers()
            # adjust binding region (if template binding region was customized)
            # important for binding region plots (relative positions)
            old.template.df <- template.df
            old.template.df$Allowed_Start_fw <- template.df$Allowed_Start_fw_initial
            old.template.df$Allowed_End_fw <- template.df$Allowed_End_fw_initial
            old.template.df$Allowed_Start_rev <- template.df$Allowed_Start_rev_initial
            old.template.df$Allowed_End_rev <- template.df$Allowed_End_rev_initial
            # update the binding positions of the templates relative to the current selected binding region
            primers <- openPrimeR:::update_primer_binding_regions(primers, template.df, old.template.df)
            rv_primers$evaluated_primers <- primers
        } else { # otherwise: ignore coverage of input primers & warn user
            toggleModal(session, "TemplateCoverageUpdateFailed")
        }
    } 
})
templateInvalidationObserver <- observeEvent(c(seq.data.input()), {
    # Invalidate previously computed rv_values when templates or binding regions
    reset.reactive.values(rv_values)
    reset.reactive.values(rv_templates, keep = c("SeqTab", "raw_seqs"))  # templates not valid anymore
    reset.reactive.values(rv_primers) # invalidate loaded primer data when templates change
    rv_cur.input.data$primers <- NULL # don't load old primer data in primer.data() function

})
########
##### Specific navigation events
######
############# TEMPLATE NAVIGATION
observeEvent(input$confirm_binding_conditions, {
    # close binding conditions panel and open constraint customization panel
    updateCollapse(session, "input_settings_collapse", close = c("binding_conditions_panel"), 
        open = "customize_settings_panel")
})

observeEvent(input$confirm_uploaded_allowed_regions, {
    # navigates from uploading allowed regions to customizing allowed regions update
    # collapse selection after allowed regions have been confirmed
    updateCollapse(session, "other_template_data_collapse", close = "allowed_template_panel")
    # open the customization of allowed regions panel
    updateCollapse(session, "customize_allowed_regions_collapse", open = "customize_allowed_regions_panel")
})

observeEvent(input$confirm_uploaded_templates, {
    # close template input panel when settings have been confirmed by user and go to
    # allowed regions
    updateCollapse(session, "template_collapse_analysis", close = c("template_input_panel"), 
        open = c("allowed_template_panel"))
})
######### PRIMER NAVIGATION
#observeEvent(input$confirm_primer_header_structure, {
    ## when primer header structure has been confirmed, move to the upload of the
    ## primer data
    #updateCollapse(session, "primer_upload_collapse", close = c("primer_header_structure_panel"), 
        #open = c("primer_upload_panel"))
#})
############# SETTINGS NAVIGATION
observeEvent(input$load_settings_choice, {
    # show panel for default (provided) settings or personal settings?
    open <- NULL
    closed <- NULL
    if (input$load_settings_choice == "default") {
        # open default tab
        open <- "load_default_settings_panel"
        closed <- "load_personal_settings"
    } else if (input$load_settings_choice == "personal") {
        closed <- "load_default_settings_panel"
        open <- "load_personal_settings"
    } else {
        message(paste("Not implemented settings load choice: ", input$load_settings_choice, 
            sep = ""))
    }
    updateCollapse(session, "load_settings_collapse", open = open, close = closed)
})

observeEvent(input$confirm_settings_choice, {
    # close xml setting choice panel when settings have been confirmed by user
    updateCollapse(session, "input_settings_collapse", close = c("load_settings_panel"), 
        open = "binding_conditions_panel")
})
observeEvent(input$confirm_PCR_settings, {
    # open the Other settings bsCollapse
    updateCollapse(session, "input_settings_collapse", close = c("customize_general_settings_panel"), 
        open = "customize_other_settings_panel")
})


observeEvent(input$confirm_constraints, {
    # close constraint customization panel when user confirms
    updateCollapse(session, "input_settings_collapse", close = c("customize_settings_panel"), 
        open = "customize_general_settings_panel")
})
###### DESIGN NAVIGATION
algorithmChoiceObserver <- observeEvent(input$confirm_algorithm_choice, {
    # move from algo selection to algo options panel
    updateCollapse(session, "design_options_algorithms_collapse",
        close = "design_algo_panel", open = "design_optimize_binding") 
})
