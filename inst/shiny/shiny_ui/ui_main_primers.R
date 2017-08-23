########
# Primer main UI panel
########
tabPanel("Primers", 
    icon = icon("tag", lib = "glyphicon"),
    br(),
    # individual primer cvg view
    radioButtons("view_cvg_individual",
    tagList(icon = icon("television"),
        "Show individual primer coverage"), 
    c("On" = "active", "Off" = "inactive"), 
    selected = "inactive", inline = TRUE
    ),
    bsTooltip("view_cvg_individual", 
    "Should the primer table show the coverage of individual templates or groups of templates?",
    "right", options = list(container = "body")
    ),
    DT::dataTableOutput("PrimerTab"), 
    br(),
    div(class="rightAligned",
        style = "color:blue; font-weight:bold", 
        htmlOutput("runModeText")
    )
    #### 
    # old idea: show additional primer info on click?
    #bsModal("PrimerDetail", "Primer Details", "",  # add more info than just the table...
            #size = "large", DT::dataTableOutput("primer_detail_table"))
)
