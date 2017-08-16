######### Shiny server functionalities for dimerization

# rv_dimer.data: selected_idx: selected idx from dimerization table to be shown
# in detail using popup
rv_dimer.data <- reactiveValues(selected_idx = NULL)

#dimerDetailObserver <- observeEvent(input$dimer_data_rows_selected, {
    ## if something is clicked, show details of dimerization (all conformations)
    #sel.ID <- input$dimer_data_row_last_clicked  # last clicked row: since ID is the first column, we need to match to ID
    ## store in reactive rv_values to access by reactive function
    #rv_dimer.data$selected_idx <- as.numeric(sel.ID)  # only works if rownames are reset to 1:N
    #toggleModal(session, "DimerDetail")  # show modal when selection changes
#})
#cur_dimer_detail <- reactive({
    ## show the current dimer details for selection
    #
    ## Need to differentiate self-comp and cross-comp because of two indices/one index
    #if (length(rv_dimer.data$selected_idx) == 0 || length(input$selected_dimerization_data) == 
        #0) {
        #data <- NULL
    #} else {
        ## show info on selected dimer idx
        #data <- NULL
        ## check input> self/cross
        #if (input$selected_dimerization_data == "Self-Dimerization") {
            #dimer.ID <- cur.dimer.data()$ID[rv_dimer.data$selected_idx]
            #idx <- which(all.self.dimer.data()$ID == dimer.ID)
            #data <- all.self.dimer.data()[idx, ]
            #data <- openPrimeR:::view.self.dimer.table(data)
        #} else if (input$selected_dimerization_data == "Cross-Dimerization") {
            #dimer.ID <- cur.dimer.data()[rv_dimer.data$selected_idx, "ID 1"]
            #dimer.ID2 <- cur.dimer.data()[rv_dimer.data$selected_idx, "ID 2"]
            #idx <- which(all.cross.dimer.data()$ID_1 == dimer.ID & all.cross.dimer.data()$ID_2 == 
                #dimer.ID2)
            #data <- all.cross.dimer.data()[idx, ]
            #data <- openPrimeR:::view.cross.dimer.table(data)
        #}
        #o <- order(data$DeltaG)
        #data <- data[o, ]
        ## remove columns due to reduced with of modal
        #data <- openPrimeR:::exclude.cols(c("Match Length", 
                                            #"Mismatch Length"), data)
    #}
    #return(data)
#})
#output$dimer_detail <- DT::renderDataTable({
    ## show a table for the currently selected primers dimerization status
    #validate(need(input$selected_dimerization_data, "Please select a type of dimerization."))
    #validate(need(cur_dimer_detail(), "No dimerizing primers available."))
    #cutoff <- ifelse(input$selected_dimerization_data == "Self-Dimerization", input$allowed_self_dimerization, 
        #allowed.cross.dimerization.setting())
    #DT::formatStyle(DT::datatable(cur_dimer_detail(), caption = "", 
                        #escape = FALSE, rownames = FALSE),
                    #"DeltaG", 
                    #backgroundColor = DT::styleInterval(cutoff, 
                                      #c("#ff9999", "#99d6ff")))
    #})
output$dimer_text <- renderUI({
    # output an overview text for dimerization
    text <- HTML(paste("<h3>", 
                dimer.text.info(cur.dimer.data(), primer.data(),
                                cur.dimer.cutoff()), 
                "</h3>", sep = ""))
    return(text)
})

#cur.all.dimer.data <- reactive({
    ## all dimer data for current selection (self-dimerization/cross-dimerization)
    #if (input$selected_dimerization_data == "Self-Dimerization") {
        #data <- all.self.dimer.data()
    #} else if (input$selected_dimerization_data == "Cross-Dimerization") {
        #data <- all.cross.dimer.data()
    #} else {
        #data <- NULL
    #}
    #if (length(data) == 0) {
        #return(NULL)
    #}
    #return(data)
#})
cur.dimer.data <- reactive({
    # formatted current primer data
    if (input$selected_dimerization_data == "Self-Dimerization") {
        data <- self.dimer.data()
    } else if (input$selected_dimerization_data == "Cross-Dimerization") {
        data <- cross.dimer.data()
    } else {
        data <- NULL
    }
    if (length(data) == 0) {
        return(NULL)
    }
    o <- order(data$DeltaG)  # order by DeltaG by default
    data <- data[o, ]
    rownames(data) <- NULL  # to reset the row names
    # reorder: ID var in the front
    idx <- grep("Primer", colnames(data))
    out <- data
    if (length(idx) != 0) {
        out <- cbind(data[,idx], data[, -idx])
        colnames(out) <- c(colnames(data)[idx], colnames(data)[-idx])
    }
    return(out)
})
output$dimer_table <- DT::renderDataTable({
    # show a table of dimerization info
    data <- NULL
    if (length(input$selected_dimerization_data) != 0 &&
        length(cur.dimer.data()) != 0) {
        data <- openPrimeR:::dimerization.table(cur.dimer.data(), 
                cur.dimer.cutoff(), input$selected_dimerization_data)
    } 
    DT::datatable(data, caption = "Dimerization frequencies of individual primers with their mean free energies.", escape = FALSE)
})
cur.dimer.cutoff <- reactive({
    # cutoff for identifying dimerization events
    cutoff <- ifelse(input$selected_dimerization_data == "Self-Dimerization", input$allowed_self_dimerization, 
        input$allowed_cross_dimerization)
    return(cutoff)
})

output$dimer_data <- DT::renderDataTable({
    # table of individual dimerization events
    validate(need(input$selected_dimerization_data, "Please select a type of dimerization."))
    validate(need(cur.dimer.data(), "No dimerizing primers available."))
    # rownames are active for datatable output for unique row identification
    # (otherwise we can't identify pairs!)
    cutoff <- cur.dimer.cutoff()
    dimer.data <- cur.dimer.data()
    excl.cols <- c("Idx1", "Idx2")
    rm.idx <- sapply(excl.cols, function(x) grep(x, colnames(dimer.data)))
    if (length(rm.idx) != 0) {
        dimer.data <- dimer.data[, -rm.idx]
    }
    DT::formatStyle(DT::datatable(dimer.data, caption = "Overview of possible primer dimers. Primers with dimerizing conformations (exceeding the DeltaG threshold) are highlighted in red, while non-dimerizing conformations are shown in blue.", escape = FALSE, selection = "single"),
        "DeltaG", backgroundColor = DT::styleInterval(cutoff, 
        c("#ff9999", "#99d6ff")))
})

self.dimer.data <- reactive({
    # show worst-case structures for each primer from all.self.dimer.data
    primer.data <- switch(input$set_meta_selector, 
        "all" = primer.data(), 
        "filtered" = current.filtered.primers(), 
        "optimized" = optimal.primers())
    annealing.temp <- isolate(annealing.temperature())
    data <- openPrimeR:::compute.all.self.dimers.frontend(primer.data, primer.concentration(), Na.concentration(), 
        Mg.concentration(), K.concentration(), Tris.concentration(), annealing.temp, for.shiny = TRUE)
    if (length(data) != 0) {
        # annotate with primerID / rename Self_Dimer_DeltaG to DeltaG 
        data$Primer <- primer.data$ID[data$Self_Dimer_Idx1]
        colnames(data)[colnames(data) == "Self_Dimer_DeltaG"] <- "DeltaG"
        data <- na.omit(data)
    }
    return(data)
})
output$dimer_distribution <- renderPlot({
    # histogram of dimerization deltaG values
    validate(need(input$selected_dimerization_data, "Please select a type of dimerization."))
    validate(need(cur.dimer.data(), "No dimerizing primers available."))
    openPrimeR:::plot.dimer.dist(cur.dimer.data(), cur.dimer.cutoff())
})
all.cross.dimer.data <- reactive({
    # all cross dimer data, possibly multiple per primer (database)
    primer.data <- switch(input$set_meta_selector, 
        "all" = primer.data(), 
        "filtered" = current.filtered.primers(), 
        "optimized" = optimal.primers())
    annealing.temp <- isolate(annealing.temperature())
    data <- openPrimeR:::compute.all.cross.dimers.unfiltered(primer.data, primer.concentration(), 
        Na.concentration(), Mg.concentration(), K.concentration(), Tris.concentration(), annealing.temp, for.shiny = TRUE)
    return(data)
})
cross.dimer.data <- reactive({
    # worst-case cross dimer data per primer
    primer.data <- switch(input$set_meta_selector, 
        "all" = primer.data(), 
        "filtered" = current.filtered.primers(), 
        "optimized" = optimal.primers())
    annealing.temp <- isolate(annealing.temperature())
    data <- openPrimeR:::compute.all.cross.dimers(primer.data, primer.concentration(), Na.concentration(), 
        Mg.concentration(), K.concentration(), Tris.concentration(), annealing.temp, 
        all.cross.dimer.data(), for.shiny = TRUE)
    data <- na.omit(data)
    return(data)
})

