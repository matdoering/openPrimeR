#############
# Set paths for the shiny app 
###############
# frontend directory 
# important: if we don't use devtools, then the path MAY be set to the R library path instead of the local path -> tools/ won't be found!
# avoid using system.file!
if (exists("SHINY.PATH")) {
    appDir <- SHINY.PATH
} else {
    appDir <- system.file("shiny", package = "openPrimeR")
}
base.dir <- normalizePath(file.path(appDir, "..", "..", "..", ".."))
#message("Base dir is: ", base.dir)
tool.folder <- file.path(base.dir, "tools")
tool.src.folder <- file.path(base.dir, "tools_src")
app.data.folder <- file.path(appDir, "..", "extdata") # data shipped with the app
src.ui.folder <- file.path(appDir, "shiny_ui")
######## IMGT folders
TOOL.SOURCE.FILES <- list.files(tool.src.folder, pattern = "\\.tar", full.names = TRUE)
# MELT
MELT.PATH <- file.path(tool.folder, "MELTING")
MELT.SRC <- TOOL.SOURCE.FILES[grep("MELTING", TOOL.SOURCE.FILES)[1]]
MELT.LOCATION <- file.path(MELT.PATH, "executable", "melting-batch")
MELT.TANDEM.LOCATION <- file.path(MELT.PATH, "Data", "AllawiSantaluciaPeyret1997_1998_1999tanmm.xml")  # original tandem mismatch parameters
MELT.TANDEM.LOCATION.MOD <- file.path(MELT.PATH, "Data", "AllawiSantaluciaPeyret1997_1998_1999tanmm_mod.xml")  # modified tandem mismatch parameters (assume mean)
# VIENNARNA
VIENNA.PATH <- file.path(tool.folder, "ViennaRNA")
VIENNA.SRC <- TOOL.SOURCE.FILES[grep("ViennaRNA", TOOL.SOURCE.FILES)[1]]
VIENNA.SRC.WIN <- list.files(tool.src.folder, full.names = TRUE)
VIENNA.SRC.WIN <- VIENNA.SRC.WIN[grep("ViennaRNAWin", VIENNA.SRC.WIN)[1]]
VIENNA.LOCATION <- file.path(VIENNA.PATH, "bin", "RNAfold")
if (Sys.info()["sysname"] == "Windows") {
	VIENNA.LOCATION <- VIENNA.PATH
}
VIENNA.LOCATION.PARAM <- file.path(VIENNA.PATH, "share", "ViennaRNA")
if (Sys.info()["sysname"] == "Windows") {
	VIENNA.LOCATION.PARAM <- file.path(VIENNA.PATH, "Misc")
}
# OLIGOARRAYAUX
OLIGOARRAY.INSTALL.PATH <- file.path(tool.folder, "oligoarrayaux")
OLIGOARRAY.PATH <- file.path(OLIGOARRAY.INSTALL.PATH, "bin")
OLIGOARRAY.SRC <- TOOL.SOURCE.FILES[grep("oligoarrayaux", TOOL.SOURCE.FILES)[1]]
OLIGOARRAY.SRC.WIN <- list.files(tool.src.folder, full.names = TRUE)
OLIGOARRAY.SRC.WIN <- OLIGOARRAY.SRC.WIN[grep("OligoArrayAuxWin",  OLIGOARRAY.SRC.WIN)[1]]
# MAFFT:
MAFFT.INSTALL.PATH <- file.path(tool.folder, "MAFFT")
MAFFT.PATH <- file.path(MAFFT.INSTALL.PATH, "bin", "mafft")
if (Sys.info()["sysname"] == "Windows") {
	MAFFT.PATH <- file.path(MAFFT.INSTALL.PATH, "ms", "bin")
}
MAFFT.SRC <- TOOL.SOURCE.FILES[grep("mafft", TOOL.SOURCE.FILES)[1]]
MAFFT.SRC.WIN <- list.files(tool.src.folder, pattern = "\\.zip", full.names = TRUE)
MAFFT.SRC.WIN <- MAFFT.SRC.WIN[grep("mafft", MAFFT.SRC.WIN)[1]]
# PANDOC:
PANDOC.INSTALL.PATH <- file.path(tool.folder, "Pandoc")
PANDOC.SRC.WIN <- file.path(tool.src.folder, "pandocWin")
PANDOC.SRC.MAC <- list.files(tool.src.folder, pattern = "\\.pkg", full.names = TRUE)
PANDOC.SRC.MAC <- PANDOC.SRC.MAC[grep("pandoc",PANDOC.SRC.MAC)[1]]
PANDOC.SRC <- list.files(tool.src.folder, pattern = "\\.deb", full.names = TRUE)
PANDOC.SRC <- PANDOC.SRC[grep("pandoc",PANDOC.SRC)[1]]
if (Sys.info()["sysname"] == "Darwin") {
    # mac
    PHANTOMJS.LOCATION <- file.path(tool.src.folder, "phantomjs-2.1.1-macosx")
} else if (Sys.info()["sysname"] == "Windows") {
	PHANTOMJS.LOCATION <- file.path(tool.src.folder, "phantomjs-2.1.1-windows")
} else {
    if (.Machine$sizeof.pointer == 8) {
        # 64 bit
        PHANTOMJS.LOCATION <- file.path(tool.src.folder, "phantomjs-2.1.1-linux-x86_64")
    } else {
        # 32 bit
        PHANTOMJS.LOCATION <- file.path(tool.src.folder, "phantomjs-2.1.1-linux-i686")
    }
}
########### help folders
help.folder <- file.path(appDir, "help")
www.folder <- file.path(appDir, "www")
help.eval.folder <- file.path(help.folder, "evaluation")
help.constraint.folder <- file.path(help.folder, "constraints")
help.filter.folder <- file.path(help.constraint.folder, "filters")
help.cvg.con.folder <- file.path(help.constraint.folder, "coverage")
help.opti.folder <- file.path(help.folder, "optimization")
help.opti.filter.folder <- file.path(help.constraint.folder, "optimization")
help.opti.options.folder <- file.path(help.opti.folder, "options")
help.opti.init.folder <- file.path(help.opti.options.folder, "initialization")
help.opti.opti.folder <- file.path(help.opti.options.folder, "optimization")
FAQ.folder <- file.path(help.folder, "FAQ")
help.input.folder <- file.path(help.folder, "input")
help.input.template.folder <- file.path(help.input.folder, "templates")
help.input.primer.folder <- file.path(help.input.folder, "primers")
