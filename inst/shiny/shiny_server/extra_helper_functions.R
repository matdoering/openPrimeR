###
# Shiny helper functions
###
# HELP FUNCTIONS
view.help <- function(session) {
    updateTabsetPanel(session, "main", selected = "help")
}
# input: template help
view.template.help <- function(session) {
    view.help(session)
    updateNavlistPanel(session, "help_panel", selected = "help_templates")
}
view.template.help.element <- function(session, id) {
    view.template.help(session)
    updateTabsetPanel(session, "help_input_tabset", selected = id)
}
# input: primer help
view.primer.help <- function(session) {
    view.help(session)
    updateNavlistPanel(session, "help_panel", selected = "help_input_primers")
}
view.primer.help.element <- function(session, id) {
    view.primer.help(session)
    updateTabsetPanel(session, "help_input_primer_tabset", selected = id)
}
##############
# constraints
view.constraints.help <- function(session) {
    view.help(session)
    updateNavlistPanel(session, "help_panel", selected = "help_constraints")
}
view.PCR.help <- function(session) {
    view.constraints.help(session)
    updateTabsetPanel(session, "help_constraints_tab", selected="help_settings_PCR")
}
view.constraint.general.help <- function(id, session) {
    view.constraints.help(session)
    updateTabsetPanel(session, "help_constraints_tab", selected=id)

}
# constraints: filter 
view.filter.help <- function(id, session) {
    view.constraints.help(session)
    updateTabsetPanel(session, "help_constraints_tab", selected="help_filtering_constraints")
    # print(paste0("viewing filter help: ", id))
    updateTabsetPanel(session, "help_filters_tab", selected=id)
}
# view cvg constraint help
view.cvg.help <- function(id, session) {
    view.constraints.help(session)
    updateTabsetPanel(session, "help_constraints_tab", selected= "help_coverage_conditions")
    updateTabsetPanel(session, "help_cvg_constraints_tab", selected=id)
}

# constraints: optimization
view.opti.constraints.help <- function(session, id) {
   view.constraints.help(session)
    updateTabsetPanel(session, "help_constraints_tab", selected="help_opti_constraints")
    updateTabsetPanel(session, "help_opti_tab", selected=id)
}
##########
# evaluation help
view.eval.help <- function(session) {
    view.help(session)
    updateNavlistPanel(session, "help_panel", selected = "help_evaluation")
}
view.eval.help.entry <- function(session, id) {
    view.eval.help(session)
    updateTabsetPanel(session, "help_evaluation_tab", selected=id)
}
###
# compare help
###
view.compare.help <- function(session) {
    view.eval.help(session)
    updateTabsetPanel(session, "help_evaluation_tab",   selected="help_comparison")
}

#view.compare.help.entry <- function(session, id) {
#    view.compare.help(session)
#    updateTabsetPanel(session, "help_compare_tab", selected = id)
#}
#################
# optimization help
###
view.opti.help <- function(session) {
    view.eval.help(session)
    updateTabsetPanel(session, "help_evaluation_tab",   selected="help_optimization")
}
view.opti.help.element <- function(session, id) {
    view.opti.help(session)
    updateTabsetPanel(session, "help_optimization_operations", selected = id)
}
view.opti.init.help <- function(session) {
    view.opti.help(session)
    updateTabsetPanel(session, "help_optimization_operations", selected="opti_init")

}
view.opti.init.help.element <- function(session, id) {
    view.opti.init.help(session)

    updateTabsetPanel(session, "help_optimization_init",  selected = id)
}
############# # end of help view functions####
create.help.button <- function(id) {
    button.name <- paste("help_", id, sep = "")
    a <- actionButton(button.name, "", icon = icon("info-sign", lib = "glyphicon"),
        style='padding:0px; font-size:85%; border-width:0px; opacity:0%; background-color:#f5f5f5')
    return(a) 
}
create.anchor <- function(id) {
    # creates a reference to an id using anchors
    #url <- paste("127.0.0.1:", port, sep = "")
    url <- ""
    link <- paste("<a href='#", id, "'>[?]</a>", sep = "")
    return(HTML(link))
}
switch.primer.view <- function(current.tab) {
    print("switching primer view function..")
    # switches to the Primers tab if current.tab is not a 'display tab' (PRIMER.TABS)
    #   current.tab: currently selected main tab in UI
    if (current.tab %in% PRIMER.TABS) {
        return(current.tab)
    } else {
        return("Primers")
    }
}

unset.subprocess.busy <- function(session) {
    # unset busy status
    # remove the busy class (modal can be closed by user, blue background (style.css) disappears.
    shinyjs::removeClass(id = "BusyInfo", class = "busy") 
    toggleModal(session, "BusyInfo", toggle = "toggle")
}
set.subprocess.busy <- function(session) {
    # set busy status
    shinyjs::addClass(id = "BusyInfo", class = "busy")
    #session$sendCustomMessage(type='jsCode', list(value = "$('#BusyInfo').modal('toggle')")) # show busy modal.
    toggleModal(session, "BusyInfo", toggle = "toggle")
}
update.following.navigation <- function(session, cur.phase, action) {
    # enable/disable following tabs depending on user action
    # cur.phase: templates/primers/settings
    # action: enable disable
    debug <- FALSE # if debug is on, navigation is not restricted here
    if (debug) {
        message("DEBUG mode active: update.following.navigation")
    } else {
        panels <- c("settingsPanel" = "primer_input_tab", 
                    "settingsPanel" = "constraint_panel",
                    "settingsPanel" = "analyze_panel")
         #message("updating tab navigation ...")
         if (action == "disable") { # disable all following navigation elements
             if (cur.phase == "templates") {
                panels <- panels
             } else if (cur.phase == "primers") {
                panels <- panels[2:length(panels)]
             } else if (cur.phase == "settings") {
                panels <- panels[3:length(panels)]
             } else {
                message(paste("nothing to change for phase: ", cur.phase, sep = ""))
             }
        } else { # action is "enable" -> activate only the following navigation element
            if (cur.phase == "templates") {
                panels <- panels[1]
             } else if (cur.phase == "primers") {
                panels <- panels[2]
             } else if (cur.phase == "settings") {
                panels <- panels[3]
             } else {
                message(paste("nothing to change for phase: ", cur.phase, sep = ""))
             }
        }
         for (i in seq_along(panels)) {
            panel_selector = paste("#", names(panels)[i], " li a[data-value=", panels[i], "]", sep = "")
            if (action == "disable") {
                #message("disabling:")
                shinyjs::disable(selector = panel_selector) 
                shinyjs::addClass(class = "disabled", selector = panel_selector) # add disabled class to grey out/change cursor (not done by shinyjs for selectors!!)
            } else {
                #message("enabling")
                shinyjs::enable(selector = panel_selector) 
                shinyjs::removeClass(class = "disabled", selector = panel_selector) # remove greyed-out color
            }
            #message(panel_selector)
         }
    }
}
create.constraint.table.row.custom <- function(radio.button, slider.setting, slider.limit, help.element = NULL) {
    na.element <- "" # could do some other style (icon?)
    if (length(radio.button) == 0) {
        radio.button <- na.element 
    }
    if (length(slider.setting) == 0) {
        slider.setting <- na.element    
    }
    nbr.cols <- 3
    if (length(slider.limit) == 0) {
        # don't show slider limit column
        nbr.cols <- 2
    }
    if (length(help.element) == 0) {
        help.element <- na.element
    }
    if (nbr.cols == 3) { 
        out <- tagList(HTML("<tr>
                    <td>"), 
                    radio.button,
                HTML("</td>",
                    "<td>"),
                    slider.setting,
                HTML("</td>
                    <td>"
                ),
                slider.limit,
                HTML("</td></tr>"))
    } else if (nbr.cols == 2) {
         out <- tagList(HTML("<tr>
                    <td>"), 
                    radio.button,
                HTML("</td>",
                    "<td>"),
                    slider.setting,
                HTML("</td></tr>"))
    } else {
        out <- NULL
    }
    return(out)
}

create.Tm.text <- function(Tm.range, Tm.diff) {
    Tm.range <- round(Tm.range, 2)
    Tm.diff <- round(Tm.diff, 2)
    Tm.r <- paste("from ", Tm.range[1], "&#8451; to ", Tm.range[2], "&#8451;", sep = "")
    text <- paste("The melting temperature ranges ", Tm.r, " (&Delta;T<sub>m</sub> = ", Tm.diff, "&#8451;).", sep = "")
    return(text)
}

reset.reactive.values <- function(values, keep = NULL) {
    # reset an input reactive values list to default values.
    # args:
    #   values: reactive values
    #   keep: vector of names in values to keep
    for (i in seq_len(length(names(values)))) {
        name <- names(values)[i]
        if (name %in% keep) { # don't remove
            next
        }
        c <- class(values[[name]])
        if (c != "NULL" && length(attributes(values[[name]])) == 0) {
            exp <- paste(class(values[[name]]), "(1)", sep = "")
            default.val <- eval(parse(text=exp))
        } else { # NULL or other non-basic class
            default.val <- NULL
        }
        values[[name]] <- default.val
        #print(paste0("reset.reactive.values():", name, " set to: ", default.val))
    }
}

switch.view.selection <- function(selected, tab, session) {
    # when frontend computed something (eval/filtering/optimization), change all selection sliders to the corresponding results
    #   selected: the results to be shown (all/filtered/optimized)
    #   tab: current main tab selected (input$main)
    #   session: current shiny session object
    #updateTabsetPanel(session, "main", selected = switch.primer.view(tab))
    updateSelectInput(session, "set_meta_selector", selected = selected)
}

#' Error/Warning Handler.
#'
#' Evaluates an expression without throwing warning/errors.
#'
#' Instead of throwing warnings/errors, the warnings/errors are
#' returned in a list object, alongside with the results.
#' 
#' @param expr Expression to evaluate
#' @return List containing values, warnings, and errors.
#' @keywords internal
withWarnings <- function(expr) {
    # call a function and store all warning in a list item "warnings". if we used tryCatch, execution would stop in case of warnings which we don't want
    myWarnings <- NULL
    myErrors <- NULL
    res <- withCallingHandlers(tryCatch(expr, error = function(e) {
      myErrors <<- append(myErrors, list(e))
      NULL
    }), warning = function(w) {
      myWarnings <<- append(myWarnings, list(w))
      invokeRestart("muffleWarning")
    })
	return(list(value = res, warnings = myWarnings, errors = myErrors))
} 
#' Custom Catch Function
#'
#' Returns a data frame if the input expression causes an error, 
#' otherwise returns the value of the evaluated expression.
#' 
#' @param expr Expression to evaluate.
#'
#' @return The evaluated expression.
#' @keywords internal
myCatch <- function(expr) {
    data <- try(expr())
    if (class(data)[1] == "try-error") {
        return(data.frame())
    } else {
        return(expr)
    }
}
combine.settings <- function(filter.settings, opti.settings, analysis.type) {
    # combines filtering and opti constraint settings for output
    constraints <- filter.settings
    overlap <- intersect(names(filter.settings), names(opti.settings))
    remainder <- setdiff(names(opti.settings), overlap)
    if (analysis.type == "design" && length(overlap) != 0) {
        # overwrite the filtering settings with the opti settings when designing
        constraints[overlap] <- opti.settings[overlap]
   } 
   if (length(remainder) != 0) {
        # only add those opti settings that aren't in the filters yet.
        constraints <- c(constraints, opti.settings[remainder])
   }
   return(constraints)
}
my_column <- function (width, ..., offset = 0)  # TODO modify
{
    if (!is.numeric(width) || (width < 1) || (width > 12)) 
        stop("column width must be between 1 and 12")
    colClass <- paste0("col-sm-", width)
    if (offset > 0) 
        colClass <- paste0(colClass, " col-sm-offset-", offset)
    div(class = colClass, ...)
}
traffic_light <- function() {
    code <- '<div id="light" style = display: inline-block;">
        <span id="red"></span>
        <span id="orange"></span>
        <span id="green"></span>
    </div>'
    return(HTML(code))
}
