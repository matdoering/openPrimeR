############
# Helper functions for shiny app installation script
##############
myip <- function(...) {
    return(readLines("http://ip-api.com/line"))
}
code.to.country <- function(code, countries) {
    idx <- match(code, countries$Code)
    return(countries[idx, "CountryShort"])
}

get.close.countries <- function(country.code, K = NULL, tool.data.folder = NULL) {
    # get the K closest country codes load country dists
	if (length(tool.data.folder) == 0) {
		tool.data.folder <- app.data.folder # use global var
	} 
	load(file.path(tool.data.folder, "country_dists.Rdata"))
    m <- match(country.code, rownames(country.dists))
    if (is.na(m)) {
        return(NULL)
    }
    dists <- country.dists[m, -m]
    o <- order(dists)
    if (length(K) == 0) {
        K <- length(dists)
    }
    out <- colnames(country.dists)[-m][o[1:K]]
    return(out)
}
set.CRAN.repository <- function(rep) {
    # sets the cran repository for installing libraries rep: selected repository
    if (length(rep) != 0) {
        local({
            r <- getOption("repos")
            r["CRAN"] <- rep
            options(repos = r)
        })
    } else {
        return(NULL)
    }
}
cranometer <- function(ms = getCRANmirrors(all = FALSE, local.only = FALSE), country = tolower(myip()[3]), tool.data.folder = NULL) {
    options(timeout = 2)
    dest = tempfile()
	# remove https mirror so we don't consider duplicates and don't have a problem if not supported by user machine
	ms <- ms[!grepl("https", ms$Name),]
    idx <- which(ms$CountryCode == country)
    # message(country)
    if (length(idx) != 0) {
        ms <- ms[idx, ]
    } else {
        # choose a close country that is contained
        countries <- get.close.countries(country, tool.data.folder = tool.data.folder)
        if (length(countries) == 0) {
            # no country found
            return(NULL)
        }
        m <- lapply(countries, function(x) which(ms$CountryCode == x))
        idx <- which(sapply(m, length) != 0)
        # use mirrors from the closest neighboring country
        if (length(idx) == 0) {
            warning("No countries close by found ...")
            return(NULL)
        }
        nbr.countries <- 1
        sel.mirror.idx <- unlist(m[idx][1:(1 + nbr.countries - 1)])
        ms <- ms[sel.mirror.idx, ]
    }
    nms = dim(ms)[1]
    ms$t = rep(NA, nms)
    for (i in 1:nms) {
        m = ms[i, ]
        url = paste(m$URL, "/src/base/NEWS", sep = "")
        t = try(system.time(try(download.file(url, dest), silent = TRUE), gcFirst = TRUE), 
            silent = TRUE)
        if (file.exists(dest)) {
            file.remove(dest)
            ms$t[i] = t["elapsed"]
        } else {
            ms$t[i] = NA
        }
    }
    ms$t <- as.numeric(ms$t)
    o <- order(ms$t)
    ms <- ms[o, ]
    return(ms)
}
my_trimws <- function (x, which = c("both", "left", "right")) {
    which <- match.arg(which)
    mysub <- function(re, x) sub(re, "", x, perl = TRUE)
    if (which == "left")
        return(mysub("^[ \t\r\n]+", x))
    if (which == "right")
        return(mysub("[ \t\r\n]+$", x)) 
    mysub("[ \t\r\n]+$", mysub("^[ \t\r\n]+", x))
}                    
get_deps <- function(path, dependencies = NA) {
    dcf <- read.dcf(file.path(path, "DESCRIPTION"))
	if (is.na(dependencies[1])) {
		dependencies <- c("Depends", "Imports", "Suggests")
	} else {
		dependencies <- dependencies
	}
    jj <- intersect(dependencies, colnames(dcf))
    val <- unlist(strsplit(dcf[, jj], ","), use.names=FALSE)
    val <- gsub("\\s.*", "", my_trimws(val))
    result <- val[val != "R"]
	return(result)
}

get.OS <- function() {
    if (grepl("^darwin", R.version$os)) {
        return("mac")
    } else if (grepl("windows", .Platform$OS.type)) {
        return("win")
    } else if (grepl("unix", .Platform$OS.type)) {
        return("linux")
    } else {
        return(.Platform$OS.type)
    }
}
check.libPaths <- function() {
    # checks whether user has a personal library or not. if not, personal library is
    # added (should also work for mac) determine OS
    #cat("o Determining personal library path ...\n")
    is.MAC <- get.OS() == "mac"
    maj.version <- R.version$major
    min.version <- strsplit(R.version$minor, split = "\\.")[[1]][1]
    version <- paste(maj.version, ".", min.version, sep = "")
    libPath <- NULL
    if (is.MAC) {
        libPath <- file.path("Library", "R", version, "library")
    } else {
        libPath <- file.path("R", paste(R.version$platform, "-library", sep = ""), 
            version)
    }
    libPath <- file.path(path.expand("~"), libPath)
    pers.lib.exists <- any(grepl(libPath, .libPaths()))
    if (!pers.lib.exists) {
        message(paste("Creating personal library: ", libPath, sep = ""))
        dir.create(libPath, recursive = TRUE)
        .libPaths(libPath)  # add new personal library to libpaths
    }
}
select.mirror <- function(x = NULL, tool.data.folder= NULL) {
    # x: pkg name to be installed
    mirror.speeds <- try(cranometer(tool.data.folder = tool.data.folder))
    if (length(mirror.speeds) == 0 || class(mirror.speeds) == "try-error") {
        # no internet connection -> continue
        return(NULL)
    }
    sel <- NULL
    if (!is.null(x)) {
        for (i in 1:nrow(mirror.speeds)) {
            # try if mirror really works by installing a package
            cur.repos <- mirror.speeds$URL[i]
            #message(paste("mirror select install: ", x))
            install.msg <- tryCatch(install.packages(x, repos = cur.repos, dep = TRUE), 
                warning = function(w) {
                  msg <- w$message
                  if (grepl("unable to connect", msg)) {
                    # retry with another mirror, remove this mirror
                    "error"
                  }
                })
            if (length(install.msg) == 0) {
                # no problem connecting
                sel <- i
                break
            } else {
                message("ERROR in selecting mirror: could not install/connect")
            }
        }
    } else {
		# just output the selected mirror without testing if it works for installing
        sel <- which.min(mirror.speeds$t)
    }
    return(mirror.speeds$URL[sel])
}
set.CRAN.mirror <- function(DEFAULT.MIRROR, x, tool.data.folder = NULL) {
    # choose a mirror if necessary DEFAULT.MIRROR: repository to use, if NULL
    # determine fastest x: pkg name
    used.repository <- NULL  # selected repo
	if (getOption("repos")["CRAN"] == "@CRAN@") {
		message("Selecting the fastest CRAN mirror ...")
		used.repository <- try(select.mirror(x, tool.data.folder = tool.data.folder))
		if (length(used.repository) == 0 || class(used.repository) == "try-error" && !interactive) {
            used.repository <- DEFAULT.MIRROR  # go to default
			set.CRAN.repository(used.repository)
        } else if (length(used.repository) == 0 || class(used.repository) == "try-error") {
			# let the user choose
            chooseCRANmirror()
			used.repository <- getOption("repos")["CRAN"]
		} else {
			# mirror was selected automatically correctly
			set.CRAN.repository(used.repository)
        }
		message(paste("Selected CRAN mirror: ", used.repository, sep = ""))
    } else {
		used.repository <- getOption("repos")["CRAN"]
	}
    return(used.repository)
}
pkgTest <- function(x, repository = "http://cran.uni-muenster.de/", biocLite = NULL, 
					load_namespace_only = FALSE, CRAN.pkgs = NULL, dependencies = NA,
					force_install = FALSE) {
	# force_install: install pkg also if already present
    # checks if required package is installed and loads it
    # load_namespace_only: check if the namespace exists or doesn't
    # -> determine repository automatically x: pkg name
    used.repository <- repository
	if (length(CRAN.pkgs) == 0) {
		CRAN.pkgs <- available.packages()
	}
    message("Loading package: ", x)
    if (!load_namespace_only && suppressMessages(suppressWarnings(!require(x, character.only = TRUE)))) {
        used.repository <- set.CRAN.mirror(repository, x)
        # message(used.repository) install packages
		if (!suppressWarnings(require(x, character.only = TRUE))) {
			# check for pkg CRAN availabiliity
			CRAN.pkg <- any(grepl(x, CRAN.pkgs[, "Package"]))
			if (CRAN.pkg) {
				# we have a CRAN pkg 
				install.packages(x, dependencies = dependencies)
			} else {
				# assume it's a Bioconductor pkg
				source("http://bioconductor.org/biocLite.R")
				biocLite(x)
			}
			if (!require(x, character.only = TRUE)) {
                warning(paste("Package ", x, " could not be attached!", sep = ""))
            }
        } 
    } else if (load_namespace_only && (!requireNamespace(x, quietly = TRUE) || force_install)) {
        used.repository <- set.CRAN.mirror(repository, x)
		CRAN.pkg <- any(grepl(x, CRAN.pkgs[, "Package"]))
		if (CRAN.pkg) {
			message("Installing `", x, "` from CRAN.")
			install.packages(x, dependencies = dependencies)
		} else {
			message("Installing `", x, "` from Bioconductor.")
            # biocLite pkg
            source("http://bioconductor.org/biocLite.R")
            biocLite(x)
        } 
		if (!requireNamespace(x, quietly = TRUE)) {
			warning(paste("Package ", x, " could not be loaded!", sep = ""))
		}
    }
    return(used.repository)
}
update.required.packages <- function(required.pkgs, USED.MIRROR, dependencies = NA) {
	message("Updating dependencies ...")
    if (length(USED.MIRROR) == 0) {
        # mirror hasn't been set
        USED.MIRROR <- set.CRAN.mirror(USED.MIRROR, NULL)
    }
    suppressWarnings(update.packages(oldPkgs = my_deps, repos = USED.MIRROR))
	#pkgs <- old.packages(repos = USED.MIRROR)[, "Package"]  # pkgs that can be updated
    #up.pkgs <- intersect(pkgs, required.pkgs)  # pkgs to update
    #for (i in seq_along(up.pkgs)) {
        #message("Updating: ", up.pkgs[i])
		#install.packages(pkgs = up.pkgs[i], dependencies = dependencies)
    #}
    return(USED.MIRROR)
}


# TOOL INSTALLATION
dir.copy <- function(src.dir, dest.dir, overwrite) {
    message(paste("Copying data from '", src.dir, "' to '", dest.dir, "'", sep = ""))
    file.names <- dir(src.dir)
    suppressWarnings(dir.create(dest.dir, recursive = TRUE))
    file.copy(from = file.path(src.dir, file.names), to = dest.dir, overwrite = overwrite, 
        recursive = TRUE)
}
install.RAxML <- function(RAXML.INSTALL.PATH, RAXML.SRC, tool.src.folder) {
    dir.create(RAXML.INSTALL.PATH)
    untar(RAXML.SRC, exdir = tool.src.folder)
    d <- list.dirs(tool.src.folder, recursive = FALSE)
    src.dir <- d[grep("RAxML", d)[1]]
    prievo.dir <- getwd()
    setwd(src.dir)
    make.cmd <- "make -f Makefile.SSE3.PTHREADS.gcc"
    res1 <- system(make.cmd)
    d <- list.files()
    setwd(prievo.dir)
    # copy raxml executable to tools/ folder
    raxml.path <- file.path(src.dir, d[grep("raxml", d)])
    cpy.call <- paste("cp ", raxml.path, " ", file.path(RAXML.INSTALL.PATH, d[grep("raxml", 
        d)]), sep = "")
    res2 <- system(cpy.call)
    return(res1 == 0 && res2 == 0)
}
install.viennaRNA <- function(VIENNA.PATH, VIENNA.SRC, VIENNA.SRC.WIN, tool.src.folder, VIENNA.LOCATION) {
    message("Installing viennaRNA.")
    is.win <- get.OS() == "win"
    if (is.win) {
        VIENNA.SRC <- VIENNA.SRC.WIN
    }
    if (!file.exists(VIENNA.SRC)) {
        message(paste("Could not find ", VIENNA.SRC, ". Please make sure the viennaRNA tool data is available to use DNA folding features.", 
            sep = ""))
        return(FALSE)
    }
    if (is.win) {
        # win: we have extracted the files from the .exe installer (7zip).
		# move binaries to destination
        status <- dir.copy(VIENNA.SRC, VIENNA.PATH, overwrite = TRUE)
        return(all(status))
    }
    # linux/mac: from source
    prievo.dir <- getwd()
    message("Installing ViennaRNA ...")
    # extract the tarball
    message("Extracting tarball ...")
    untar(VIENNA.SRC, exdir = tool.src.folder)
    # install the package
    d <- list.dirs(tool.src.folder, recursive = FALSE)
    install.folder <- d[grep("ViennaRNA", d)[1]]
    setwd(install.folder)
    message("Cleaning configuration ...")
    clean.cmd <- "make distclean"
    system(clean.cmd)
    message("Configuring src files ...")
    # viennaRNA can't install to paths with spaces (due to libtool/makefile) ->
    # install to tmp and then move to destination
    destination <- file.path(prievo.dir, VIENNA.PATH)
    prefix <- file.path(tempdir(), "ViennaRNA")
    # viennaRNA requires python2 if we use swig -> set python version to be sure
    # config.call <- paste('PYTHON_VERSION=2 ./configure --prefix=', prefix, sep =
    # '')
    config.call <- paste("./configure --prefix=", prefix, " --without-python --without-perl --without-kinfold --without-forester", 
        sep = "")  # without-swig: deactivate scripting interfaces -> no python dependence
    message(config.call)
    res <- system(config.call)
    if (res != 0) {
        setwd(prievo.dir)
        warning("Error during viennaRNA configuration. Please view the console output!")
        return(FALSE)
    }
    message("Compiling ...")
    res <- system("make")
    if (res != 0) {
        warning("Error during viennaRNA compilation.")
        setwd(prievo.dir)
        return(FALSE)
    }
    message("Installing to tmp ...")
    res <- system("make install")
    if (res != 0) {
        warning("Error during viennaRNA installation.")
        setwd(prievo.dir)
        return(FALSE)
    }
    setwd(prievo.dir)
    if (file.exists(VIENNA.PATH)) {
        unlink(VIENNA.PATH, recursive = TRUE)
    }
    mv.cmd <- paste("mv ", prefix, " ", tool.folder, sep = "")
    message("Moving to destination ...")
    res <- system(mv.cmd)
    if (res != 0) {
        warning("Error while moving viennaRNA to install folder.")
        return(FALSE)
    }
    # switch directory back to my dir
    return(file.exists(VIENNA.LOCATION))
}

install.melt <- function(MELT.PATH, MELT.SRC, tool.src.folder) {
    if (!file.exists(MELT.SRC)) {
        message(paste("Could not find ", MELT.SRC, ". Please place the .tar.gz into the required folder to use MELTING features.", 
            sep = ""))
        return(FALSE)
    }
    message(paste("Installing MELTING into: '", MELT.PATH, "'", sep = ""))
    # extract the tarball
    untar(MELT.SRC, exdir = tool.src.folder)
    # install the package get folder name containing extracted tarball
    d <- list.dirs(tool.src.folder, recursive = FALSE)
    install.folder <- d[grep("MELTING", d)[1]]
    dir.copy(install.folder, MELT.PATH, overwrite = TRUE)
	# install custom melt config:
	#openPrimeR:::copy.melt.config(file.path(MELT.PATH, "executable", "melting-batch"))
    return(TRUE)
}
install.oligoArrayAux <- function(OLIGOARRAY.PATH, OLIGOARRAY.INSTALL.PATH,
                        OLIGOARRAY.SRC, OLIGOARRAY.SRC.WIN, tool.src.folder) {
    is.win <- get.OS() == "win"
    if (is.win) {
        OLIGOARRAY.SRC <- OLIGOARRAY.SRC.WIN
    }
    if (!file.exists(OLIGOARRAY.SRC)) {
        message(paste("Could not find ", OLIGOARRAY.SRC, ". Please ensure that the oligoArrayAux source file is in the required folder to use DECIPHER features.", 
            sep = ""))
        return(FALSE)
    }
    message(paste("Installing oligoarrayaux into: ", OLIGOARRAY.INSTALL.PATH, sep = ""))
    if (is.win) {
		# just copy to destination, we've already extracted the binaries from the setup using the 'innounp' tool.
	    dir.copy(OLIGOARRAY.SRC, OLIGOARRAY.INSTALL.PATH, overwrite = TRUE)
        return(TRUE)
    } 
    ###
    # for linux/mac:
    ####
    # extract the tarball
    untar(OLIGOARRAY.SRC, exdir = tool.src.folder, compressed = "bzip2")
    # install the package get folder name containing extracted tarball
    prievo.dir <- getwd()
    on.exit(setwd(prievo.dir))
    d <- list.dirs(tool.src.folder, recursive = FALSE)
    install.folder <- d[grep("oligoarrayaux", d)[1]]
    message("install folder:")
    message(install.folder)
    setwd(install.folder)
    message("Configuring src files ...")
    prefix <- file.path(tempdir(), "oligoarrayaux")
    config.call <- paste("./configure --prefix=", prefix, sep = "")
    res <- system(config.call)
    if (res != 0) {
        warning("configuration of oligoarrayaux failed.")
        return(FALSE)
    }
    message("Compiling ...")
    ret <- system("make")
    if (ret != 0) {
        warning("Compilation of oligoarrayaux failed.")
        return(FALSE)
    }
    message("Installing ...")
    ret <- system("make install")
    if (ret != 0) {
        warning("Installation of oligoarrayaux failed.")
        return(FALSE)
    }
    # switch directory back to my dir
    setwd(prievo.dir)
    if (file.exists(OLIGOARRAY.INSTALL.PATH)) {
        unlink(OLIGOARRAY.INSTALL.PATH, recursive = TRUE)
    }
    mv.cmd <- paste("mv ", prefix, " ", tool.folder, sep = "")  # move is required for paths containing spaces (not supported by make!)
    message("Moving to destination ...")
    res <- system(mv.cmd)
    if (res != 0) {
        warning("Moving oligoarrayaux to destination failed.")
        return(FALSE)
    }
    return(TRUE)
}
install.pandoc <- function(PANDOC.INSTALL.PATH, PANDOC.SRC, PANDOC.SRC.WIN, PANDOC.SRC.MAC, tool.src.folder) {
    # this is not recommended -> doesn't really work in most cases due to Haskell's large number of uninstallable packages
	OS <- get.OS()
    #########
	if (OS == "win") {
		# copy pandoc binaries to tools folder
		dir.copy(PANDOC.SRC.WIN, PANDOC.INSTALL.PATH, overwrite= TRUE)
		return(TRUE)
	} else if (OS == "linux") {
        # unix: only 64 bit debian supported
        message("Only the 64 bit debian pandoc version is available for Linux.")
        dir.create(PANDOC.INSTALL.PATH, recursive = TRUE, showWarnings = FALSE)
        cmd <- paste0("ar p ", PANDOC.SRC, " data.tar.gz | tar xvz --strip-components 3 -C ", PANDOC.INSTALL.PATH)
        check <- system(cmd)
        if (check != 0) {
            warning("The following command failed: ", cmd, 
                    "Installation of pandoc failed.")
            return(FALSE)
        } else {
            return(TRUE)
        }
    } else if (OS == "mac") {
        # macOS
        old.wd <- getwd()
        on.exit(setwd(old.wd))
        setwd(tool.src.folder)
        dir.create("pandoc-extract", showWarnings = FALSE)
        setwd("pandoc-extract")
        cmd <- paste0("xar -xf ", PANDOC.SRC.MAC)
        system(cmd)
        # create usr/local files locally
        cmd <- "cat pandoc.pkg/Payload | gunzip -dc | cpio -i"
        check <- system(cmd)
        if (check != 0) {
            warning("Could not run the following command: ", cmd,
                 "pandoc could not be installed.")
            return(FALSE)
        }
        # executables are now in ./usr/bin/, man pages in ./usr/share/man
        src.dir <- file.path("usr", "local")
		dir.copy(src.dir, PANDOC.INSTALL.PATH, overwrite= TRUE)
        return(TRUE)
    } else {
        warning("Your OS (", OS, ") is not supported!")
        return(FALSE)
    }
}
install.pandoc.from.src <- function(PANDOC.INSTALL.PATH, PANDOC.SRC, PANDOC.SRC.WIN, tool.src.folder) {
    # this is not recommended -> doesn't really work in most cases due to Haskell's large number of uninstallable packages
	is.win <- get.OS() == "win"
	if (is.win) {
		# copy pandoc binaries to tools folder
		dir.copy(PANDOC.SRC.WIN, PANDOC.INSTALL.PATH, overwrite= TRUE)
		return(TRUE)
	}
    # linux/macOS: install from sources
    # 1) untar the archive
    untar(PANDOC.SRC, exdir = tool.src.folder, compressed = "gzip")
    # unnecessary
    ## move sources to the tools destination
    pandoc.dir <- list.dirs(dirname(PANDOC.SRC), full.names = TRUE)
    pandoc.dir <- pandoc.dir[grep("Pandoc", pandoc.dir)[1]]
    # 2) use cabal to install pandoc
    if (Sys.which("cabal") == "")  {
        message("`cabal` was not available on your system. ",
            "Please install the Haskell platform on your system or install ",
            "`pandoc` manually using your system's package repository.")
        return(FALSE)
    }
    old.wd <- getwd()
    # need to set working dir to pandoc dir for the next steps
    setwd(pandoc.dir)
    on.exit(setwd(old.wd))
    # a) Install dependencies for pandoc
    cmd <- "cabal install --only-dependencies"
    check <- system(cmd)
    if (check != 0) {
        warning("Could not install pandoc dependencies.")
    }
    # b) configure pandoc installation
    cmd <- paste0("cabal configure")  # for --prefix: necessary if we use copy with destdir?
    check <- system(cmd)
    if (check != 0) {
        warning("Could not configure pandoc with `cabal`.")
        return(FALSE)
    }
    # c) build & test
    cmd <- "cabal build && cabal test"
    system(cmd)
    if (check != 0) {
        warning("Could not build pandoc.")
        return(FALSE)
    }
    cmd <- paste0("cabal copy --destdir=", PANDOC.INSTALL.PATH)
    if (cmd != 0) {
        warning("Could not copy pandoc to target directory.")
        return(FALSE)
    }
	return(TRUE)
}
install.MAFFT <- function(MAFFT.INSTALL.PATH, MAFFT.SRC, MAFFT.SRC.WIN, tool.src.folder) {
    # modify MAFFT Makefile for local installation
	is.win <- get.OS() == "win"
	if (is.win) {
		MAFFT.SRC <- MAFFT.SRC.WIN
		temp_dir <- tempdir()
		unzip(MAFFT.SRC, overwrite = TRUE, exdir = temp_dir)
		dir.copy(file.path(temp_dir, "mafft-win"), MAFFT.INSTALL.PATH, overwrite = TRUE)
		return(TRUE)
	}
    first.line <- paste("PREFIX = ", MAFFT.INSTALL.PATH, sep = "")
    untar(MAFFT.SRC, exdir = tool.src.folder, compressed = "gzip")
    d <- list.dirs(tool.src.folder, recursive = FALSE)
    src.dir <- d[grep("mafft", d)[1]]
    src.folder <- file.path(src.dir, "core")
    makefile.loc <- file.path(src.folder, "Makefile")
    tmp.loc <- file.path(tempdir(), "Makefile")
    # cmd <- paste('sed -i '1s#.*#', first.line, ''# ', makefile.loc, sep = '') need
    # to use sed -e for mac os, not sed -i ... (posix compatibility)
    cmd <- paste("sed -e \"1s#.*#", first.line, "\"# ", makefile.loc, " > ", tmp.loc, 
        " && mv ", tmp.loc, " ", makefile.loc, sep = "")
    ret <- system(cmd)
    if (ret != 0) {
        warning("sed command for MAFFT didn't work. Failed installation!")
        return(FALSE)
    }
    prievo.dir <- getwd()
    setwd(src.folder)
    make.cmd <- "make clean && make && make install"
    res <- system(make.cmd)
    setwd(prievo.dir)
    return(res == 0)
}
install.tools <- function(AVAILABLE.TOOLS = NULL) {
    # install third-party tools
    if (length(AVAILABLE.TOOLS) == 0) {
        message("Checking tools for installation status ...")
        AVAILABLE.TOOLS <- openPrimeR:::check.tool.installation(frontend = TRUE)
    } else {
        # don't try to install again, if we already tried to install (AVAILABLE.TOOLS is set in this case)
        return(AVAILABLE.TOOLS)
    }
    for (i in seq_along(AVAILABLE.TOOLS)) {
        tool <- names(AVAILABLE.TOOLS)[i]
        available <- AVAILABLE.TOOLS[i]
        if (available) {
            next
        }
        message(paste0("Installing: ", tool))
        if (tool == "MELTING") {
            if (file.exists(MELT.LOCATION)) {
                AVAILABLE.TOOLS["MELTING"] <- FALSE
                warning("MELT is installed, but path isn't set correctly.")
            } else if (install.melt(MELT.PATH, MELT.SRC, tool.src.folder)) {
                AVAILABLE.TOOLS["MELTING"] <- TRUE
                message("-> Melt was installed successfully.")
            } else {
                AVAILABLE.TOOLS["MELTING"] <- FALSE
                warning(paste("-> MELT could not be installed. Error during installation or missing source files in : ", 
                    MELT.SRC, ". Disabling melting temperature analysis.", sep = ""))
            }
        } else if (tool == "ViennaRNA") {
            if (file.exists(VIENNA.LOCATION)) {
                AVAILABLE.TOOLS["ViennaRNA"] <- FALSE
                warning("ViennaRNA is installed, but path isn't set correctly or checks failed (UNAFOLDDAT set?).")
            } else if (install.viennaRNA(VIENNA.PATH, VIENNA.SRC, VIENNA.SRC.WIN, tool.src.folder, VIENNA.LOCATION)) {
                AVAILABLE.TOOLS["ViennaRNA"] <- TRUE
                message("-> Successfully installed ViennaRNA.")
            } else {
                AVAILABLE.TOOLS["ViennaRNA"] <- FALSE
                warning(paste("ViennaRNA could not be installed: disabled secondary structure analysis. Check location: ", 
                    VIENNA.LOCATION, ".", sep = ""))
            }
        } else if (tool == "OligoArrayAux") {
            if (file.exists(OLIGOARRAY.PATH)) {
                warning("Tool is installed but not in the path.")
                AVAILABLE.TOOLS["OligoArrayAux"] <- TRUE
            } else if (install.oligoArrayAux(OLIGOARRAY.PATH, OLIGOARRAY.INSTALL.PATH, 
                     OLIGOARRAY.SRC, OLIGOARRAY.SRC.WIN, tool.src.folder)) {
                message("Installation of oligoarrayaux successful!")
                AVAILABLE.TOOLS["OligoArrayAux"] <- TRUE
            } else {
                AVAILABLE.TOOLS["OligoArrayAux"] <- FALSE
                warning("Installation of oligoarrayaux failed!")
            }
        } else if(tool == "Selenium") {
            if (Sys.which("Python") != "") {
                cmd <- "python -m pip install --user selenium"  # local install with --user
                installed <- system(cmd)
                if (installed != 0) {
                    warning(paste("-> Selenium for python couldn't be installed.", 
                    "Possibly pip is not available on your system.",
                    "Please refer to: https://packaging.python.org/installing/#requirements-for-installing-packages"))
                } else {
                    message("-> Successfully installed selenium for python!")
                }
            } else {
                warning(paste("Did not install selenium for python as python or pip weren't available. 
                Please install python to ensure if you want to use the IMGT database feature."))
		    }
        } else if (tool == "PhantomJS") {
            s <- dir.copy(PHANTOMJS.LOCATION, file.path(tool.folder, "phantomjs"), overwrite = TRUE)
            AVAILABLE.TOOLS["PhantomJS"] <- all(s)
        } else if (tool == "MAFFT") {
            if (file.exists(MAFFT.PATH)) {
                AVAILABLE.TOOLS["MAFFT"] <- FALSE
                warning("MAFFT is installed, but path isn't set correctly.")
            } else if (install.MAFFT(MAFFT.INSTALL.PATH, MAFFT.SRC, MAFFT.SRC.WIN, tool.src.folder)) {
                AVAILABLE.TOOLS["MAFFT"] <- TRUE
                message("-> Successfully installed MAFFT")
            } else {
                AVAILABLE.TOOLS["MAFFT"] <- FALSE
                warning(paste("-> MAFFT could not be installed. Check location: ", 
                    MAFFT.SRC, ".\n", sep = ""))
            }
        } else if (tool == "Pandoc") {
            if (install.pandoc(PANDOC.INSTALL.PATH, PANDOC.SRC, PANDOC.SRC.WIN,PANDOC.SRC.MAC, tool.src.folder)) {
                message(paste("Pandoc was installed successfully.\n",
                        "Please also ensure you also have a LateX distribution installed",
                        "to use PDF reporting features."))
                AVAILABLE.TOOLS["Pandoc"] <- TRUE
            } else {
                warning(paste("Pandoc could not be installed.",
                    "openPrimeR won't be able to generate PDF reports."))
                AVAILABLE.TOOLS["Pandoc"] <- FALSE
            }
        } else {
            warning(paste("Unknown tool specified:", tool))
        }    
    }
    return(AVAILABLE.TOOLS)
}
# create icon for the package
create_tool_icon <- function(base.path) {
    # base.path: the primer design base folder
	is.unix <- get.OS() == "linux"
	if (!is.unix) {
		# icon only makes sense for unix systems.
		return(NULL)
	}
    message("Creating desktop icon ...")
    icon.loc <- system.file("shiny", "www", "images", 
                    "logo.png", package = "openPrimeR")
    exec.loc <- file.path(normalizePath(base.path), "start.sh")
    if (!file.exists(icon.loc)) {
        warning(paste("Icon not found at:", icon.loc))
        return()
    }
    if (!file.exists(exec.loc)) {
        warning(paste("Start script not found at:", exec.loc))
        return()
    }
    link.loc <- file.path(base.path, "openPrimeR")
    code <- paste(
        "[Desktop Entry]",
        "Name=openPrimeR",
        "Comment=Start openPrimeR",
        paste0("Icon=", icon.loc),
        paste0("Exec=", exec.loc),
        "Type=Application", 
        "Terminal=false",
        "StartupNotify=true",
        "GenericName=Start openPrimeR",
        sep = "\n")
    writeLines(code, link.loc)
}

