##########
# Side panel for UI downloads
#######
tabPanel("Download", value = "DownloadTab",
    icon = icon("download-alt", lib = "glyphicon"),
    h2("Download"),
    div(p("Download the currently selected data.
           Raw primer/template sets are provided in CSV format, while formatted sets are provided as FASTA files. Both types of files can be used as  
           input to the tool.")),
    # input of sample name
    textInput(inputId = "sample_name", 
        label = tagList(icon("briefcase", lib = "glyphicon"), 
        "Identifier"), 
         value = "Sample"), 
    bsTooltip("sample_name", 
        "The identifier for the primer analysis.",
        "right", options = list(container = "body")),
    selectInput("downloadDataSet", 
        tagList(icon("database"), "Data set"), 
        choices = list("Summary" = "all", 
					"Reports" = c(
						"Evaluation Report" = "eval_report",
						"Comparison Report" = "comparison_report",
                        "Coverage Spreadsheet" = "coverage_spreadsheet"
					),
					"Core data" = c(
						"Templates" = "sequences", 
						"Primers" = "primers",
						"Primer subset" = "primers_subset",
						"Settings" = "settings"
					),
					"Mismatches" = c(
						"Mismatch table (fw)" = "mismatch_table_fw",
						"Mismatch table (rev)" = "mismatch_table_rev"
					),
					"Dimerization" = c(
						"Self dimerization" = "self_dimerization", 
						"Cross dimerization" = "cross_dimerization"
					)
        )
    ),
    radioButtons("download_style",  # download raw/user-formatted data
            tagList(icon("television"), "Output format"), 
            choices = c("Raw" = "raw", "Formatted" = "formatted"), inline=TRUE),
    bsTooltip("download_style", 
        "Download data wihout any special formatting (machine-redable, e.g. CSV) or with formatting (e.g. FASTA).",
        "right", 
        options = list(container = "body")
    ),
    downloadButton('downloadData', 'Download') # download button
)

