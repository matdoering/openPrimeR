################
# Shiny server functionalities for evaluation
###################

evaluated.primers <- observeEvent(input$evaluateButton, {
    # data frame of evaluated primers with annotated constraints updates reactive
    # values
    if (input$evaluateButton == 0) {
        return(NULL)
    }
    if (length(primer.data()) == 0 || length(current.seqs()) == 0) {
        # show bsmodal
        session$sendCustomMessage(type = "jsCode", list(value = "$('#NotifyNoDataAvailable').modal('show')"))
        return(NULL)
    }
    primer.df <- primer.data()
    if (length(rv_primers$evaluated_primers) != 0) {
        # if primers were evaluated -> use the evaluated data set to save some
        # computations
        primer.df <- rv_primers$evaluated_primers
    }
    con <- constraints()
    active.constraints <- con[["active"]]
    ###### Create a callback function to update progress.
    progress <- shiny::Progress$new()
    progress$set(message = "Evaluating", value = 0)
    on.exit(progress$close())
    updateProgress <- function(value, detail, option) {
        if (is.null(value)) {
            value <- progress$getValue()
            value <- value + (progress$getMax() - value)/5
        }
        if (option == "inc") {
            progress$inc(amount = value, detail = detail)
        } else {
            progress$set(value = value, detail = detail)
        }
    }
    to.compute.constraints <- active.constraints
    annealing.temp <- isolate(annealing.temperature()) # don't trigger
    settings <- current.settings() 
    eval.data <- openPrimeR:::check_constraints(primer.df, current.seqs(),
        settings, active.constraints, 
        to.compute.constraints, for.shiny = TRUE, updateProgress = updateProgress)

    rv_primers$available_constraints <- unique(c(rv_primers$last_eval_constraints, active.constraints)) # add evaluated constraints as 'available'
    rv_primers$evaluated_primers <- eval.data
    rv_primers$all <- eval.data
    if ("primer_coverage" %in% colnames(eval.data)) {
        rv_templates$cvg_all <- openPrimeR:::update_template_cvg(current.seqs(), 
            eval.data, run.mode())  # update templates with cvg values 
    }
    switch.view.selection("all", input$main, session)  # view evaluated primers
})
