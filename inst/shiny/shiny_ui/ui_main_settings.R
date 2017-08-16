########
# UI panel for the current analysis settings
#######
tabPanel("Settings",
    value = "settings_view",
    icon  = icon("wrench"),
    br(),
    #####################
    # CURRENT SETTINGS
    ######################
    selectInput("selected_settings_table", "Selected data", 
                c("Constraints" = "constraints", 
                "Coverage Constraints" = "cvg_constraints", 
                "Options" = "opts", 
                "PCR" = "PCR_options")),
    DT::dataTableOutput("current_constraints")
)
