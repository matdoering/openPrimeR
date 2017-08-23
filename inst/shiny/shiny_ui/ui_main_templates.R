######
# UI template main panel
######
tabPanel("Templates", 
    value = "template_view_panel",
    icon = icon("book", lib = "glyphicon"),
    br(),
    DT::dataTableOutput('SeqTab')
)
 
