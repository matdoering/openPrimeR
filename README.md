[![Build Status](https://travis-ci.org/matdoering/openPrimeR.svg?branch=master)](https://travis-ci.org/matdoering/openPrimeR)[![codecov](https://codecov.io/gh/matdoering/openPrimeR/branch/master/graph/badge.svg)](https://codecov.io/gh/matdoering/openPrimeR) 

# openPrimeR

## Synopsis
openPrimeR is an R package providing methods for designing, evaluating,and
comparing primer sets for multiplex polymerase chain reaction (PCR). The package provides a primer design
function that generates novel primer setes by solving a
set cover problem such that the number of covered template sequences is
maximized with the smallest possible set of primers. Moreover, existing primer sets can be evaluated
according to their coverage and their fulfillment of constraints on the
PCR-relevant physicochemical properties. For PCR tasks for which multiple
possible primer sets exist, openPrimeR can facilitate the selection of the
most suitable set by performing comparative analyses. The R package includes a Shiny application that
provides a comprehensive and intuitive user interface for the core functionalites of the package.

## More information
There is a [Conda repository for openPrimeR available](https://anaconda.org/bioconda/bioconductor-openprimer/badges).
For more information on how to install the third-party dependencies for openPrimeR, we refer to the corresponding [user-space repository](https://github.com/matdoering/openPrimeR-User), which provides installation routines.

### Docker Instructions

Our docker image is available at [dockerhub](https://hub.docker.com/r/mdoering88/openprimer/). In order to use the docker image, you need to [install docker](https://www.docker.com/) on your system and activate the docker daemon. 

After logging in with your docker account, in a console, enter 

>**docker pull mdoering88/openprimer**

to retrieve the latest docker image of the tool. Since the image is quite large the download (~4 GB) may take some time, especially if your internet connection is not very fast.

To run the image, enter

>**docker run -p 3838:3838 --rm mdoering88/openprimer**

After this, the tool is available by accessing **localhost:3838** in your web browser.

In case you want to have more control of the image you are running or you want to study the output of the tool, you can execute

>**docker run --rm -p 3838:3838 -v /tmp/logs/:/var/log/shiny-server/ mdoering88/openprimer:latest**

With this call, the tool's log file is stored in the */tmp/logs/* folder on your system with a filename starting with *shiny-shiny*). During the session, you can use 

>**tail -f shiny-shiny-X.txt**

to retrieve the current status messages. 

Moreover, in the above call to Docker, we have specified a tag, namely *:latest*, which means that we have started the most recent version of the tool. In case that no tag is provided, the latest available Docker image is used.

You can also work with the openPrimeR API through the container by executing

```
docker run -rm -it mdoering88/openprimer bash
```

and then starting an R session and loading the `openPrimeR` package.

## Changelog

Take a look at the [CHANGELOG](inst/NEWS) to view recent changes to the project.

