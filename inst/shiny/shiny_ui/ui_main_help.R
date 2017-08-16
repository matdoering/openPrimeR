######
# Help main panel in UI
######
tabPanel("Help",value="help",
    icon = icon("info-sign", lib = "glyphicon"),
    br(),
    navlistPanel( # help navigation panel on the left-hand side
        paste("openPrimeR", " v", packageVersion("openPrimeR"), sep = ""),
        tabPanel("FAQ", value="FAQ", 
            # FAQ entries
            icon = icon("question-sign", lib = "glyphicon"),
            includeHTML(file.path(FAQ.folder, "screencast.html")),
            includeHTML(file.path(FAQ.folder, "tools.html")),
            includeHTML(file.path(FAQ.folder, "ambiguities.html")),
            includeHTML(file.path(FAQ.folder, "design.html"))

        ),
        tabPanel("Templates", value = "help_templates",
            # input help
            icon = icon("book", lib = "glyphicon"),
            tabsetPanel(id = "help_input_tabset",
                tabPanel("Overview", value="help_input_templates_overview",
                          icon = icon("fullscreen", lib="glyphicon"),
                        includeHTML(file.path(help.input.template.folder, "overview.html"))
                        ),
                        tabPanel("Binding regions", value="help_input_templates_allowed",
                                 icon = icon("bookmark", lib = "glyphicon"),
                                 includeHTML(file.path(help.input.template.folder, "allowed.html"))
                        )
                    )
                 ),
                tabPanel("Primers", value = "help_input_primers",
                    # primer help panel
                    icon = icon("tag", lib = "glyphicon"),
                        tabsetPanel(id = "help_input_primer_tabset", 
                        tabPanel("FASTA Input", value = "help_input_primers_overview",
                            includeHTML(file.path(help.input.primer.folder, "FASTA_input.html"))
                        ),
                        tabPanel("CSV Input", value = "help_input_primers_comparison",
                            includeHTML(file.path(help.input.primer.folder, "CSV_input.html"))
                        )
                    )
        ),
        tabPanel("Settings", value = "help_constraints",
            # constraint help panel
            #icon  = icon("resize-small", lib = "glyphicon"),
            icon = icon("menu-hamburger", lib = "glyphicon"),
            tabsetPanel(id = "help_constraints_tab",
                tabPanel("Overview", value = "help_constraints_overview",
                    icon = icon("fullscreen", lib="glyphicon"),
                    includeHTML(file.path(help.constraint.folder, "overview.html"))
                ),
                ###
                # START CONSTRAINT SETTINGS HELP
                ####
                tabPanel("Constraints", value = "help_filtering_constraints", # tabpanels have values, tabsets have ids!
                    icon = icon("filter", lib = "glyphicon"),
                    br(),
                    tabsetPanel(id = "help_filters_tab",
                        tabPanel("Overview", value = "help_tab_overview_filters",
                          icon = icon("fullscreen", lib="glyphicon"),
                            includeHTML(file.path(help.filter.folder, "overview.html"))
                        ),
                        tabPanel("Primer coverage", value="help_tab_primer_coverage",
                            includeHTML(file.path(help.filter.folder, "primer_coverage.html"))
                        ),
                        tabPanel("Primer length",value = "help_tab_primer_length",
                            includeHTML(file.path(help.filter.folder, "primer_length.html"))
                        ),
                        tabPanel("GC clamp", value = "help_tab_gc_clamp",
                            includeHTML(file.path(help.filter.folder, "GC_clamp.html"))
                        ),
                        tabPanel("GC ratio", value = "help_tab_gc_ratio",
                            includeHTML(file.path(help.filter.folder, "GC_content.html"))
                        ),
                        tabPanel("Run length", value = "help_tab_run_length",
                            includeHTML(file.path(help.filter.folder, "run_length.html"))
                        ),
                        tabPanel("Repeat length", value = "help_tab_repeat_length",
                            includeHTML(file.path(help.filter.folder, "repeat_length.html"))
                        ),
                        tabPanel("Melting temperature", value = "help_tab_melting_temperature",
                            includeHTML(file.path(help.filter.folder, "melting_temperature.html"))
                        ),
                        tabPanel("Melting temperature deviation", value = "help_tab_melting_temperature_diff",
                            includeHTML(file.path(help.filter.folder, "melting_temp_opti.html"))
                        ),
                        tabPanel("Secondary structure", value = "help_tab_secondary_structure",
                            includeHTML(file.path(help.filter.folder, "secondary_structure.html"))
                        ),
                        tabPanel("Primer specificity", value = "help_tab_primer_specificity", 
                            includeHTML(file.path(help.filter.folder, "primer_specificity.html"))
                        ),
                        tabPanel("Cross dimerization", value = "help_tab_cross_complementarity", 
                            includeHTML(file.path(help.filter.folder, "cross_complementarity.html"))
                        ),
                        tabPanel("Self dimerization", value = "help_tab_self_complementarity", 
                            includeHTML(file.path(help.filter.folder, "self_complementarity.html"))
                        )
                    )
                ),
                #########
                # COVERAGE CONSTRAINTS
                #########
                tabPanel("Coverage Constraints", value = "help_coverage_conditions", 
                    icon = icon("star"),
                    br(),
                    tabsetPanel(id = "help_cvg_constraints_tab",
                        tabPanel("Basic conditions", value="help_tab_coverage_basic",
                            includeHTML(file.path(help.cvg.con.folder, "basic_primer_coverage.html"))
                        ),
                        tabPanel("Coverage model", value="help_tab_coverage_model",
                            includeHTML(file.path(help.cvg.con.folder, "coverage_model.html"))
                        ),
                         tabPanel("Primer efficiency", value="help_tab_primer_efficiency",
                            includeHTML(file.path(help.cvg.con.folder, "primer_efficiency.html"))
                        ),
                         tabPanel("Annealing energy", value="help_tab_annealing_DeltaG",
                            includeHTML(file.path(help.cvg.con.folder, "annealing_DeltaG.html"))
                        ),
                        tabPanel("Codon design", value="help_tab_codon_design",
                            includeHTML(file.path(help.cvg.con.folder, "codon_design.html"))
                        )
                    )
                ),
                # PCR SETTINGS HELP
                tabPanel("PCR", value = "help_settings_PCR",
                    icon = icon("flask"),
                    includeHTML(file.path(help.constraint.folder, "PCR_settings.html"))
                )
            )
        ),
        tabPanel("Analyze", value = "help_evaluation",
            # evaluation help
            icon = icon("cogs"),
            tabsetPanel(id = "help_evaluation_tab",
                tabPanel("Overview", value = "help_eval_actions",
                          icon = icon("fullscreen", lib="glyphicon"),
                    includeHTML(file.path(help.eval.folder, "overview.html"))
                ),
                tabPanel("Evaluation", value = "help_eval_eval",
                    icon = icon("bar-chart"),
                    includeHTML(file.path(help.eval.folder, "evaluation.html"))
                ),
                tabPanel("Design", value="help_optimization",
                    icon = icon("sort-amount-desc"),
                    br(),
                    tabsetPanel(id = "help_optimization_operations",
                        tabPanel("Design",value="help_opti_optimization",
                            icon = icon("sort-amount-desc"),
                            includeHTML(file.path(help.eval.folder, "optimization.html"))
                        ),
                        tabPanel("Binding regions",value="opti_templates_help",
                            icon = icon("book", lib = "glyphicon"),
                            includeHTML(file.path(help.eval.folder, "optimization_templates.html"))
                        )
                    )
                ),
                tabPanel("Comparison", value="help_comparison",
                    icon = icon("scale", lib = "glyphicon"),
                    includeHTML(file.path(help.eval.folder, "comparison.html"))
                )
            ) 
        ),
        tabPanel("Contact",
            # contact panel
            icon = icon("envelope", lib = "glyphicon"),
            br(),
            includeHTML(file.path(help.folder, "contact.html"))
        ),
        tabPanel("References", value="reference_papers",
            # reference panel
            icon = icon("search", lib = "glyphicon"),
            includeHTML(file.path(help.folder, "references.html"))

        ),
        tabPanel("License", value = "license_info",
            # license information
            icon = icon("exclamation-circle"),
            includeHTML(file.path(help.folder, "license.html"))
        ),
        tabPanel("Disclaimer", value = "disclaimer",
            # Disclaimer
            icon = icon("balance-scale"),
            includeHTML(file.path(help.folder, "disclaimer.html"))
        )
        , widths = c(3,9), id = "help_panel" # width of navigation and main part
    )
)

