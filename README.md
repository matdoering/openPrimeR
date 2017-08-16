# openPrimeR
## Synopsis

openPrimeR is a web-based tool for designing, evaluating, and comparing primers for multiplex polymerase chain reaction (PCR). When designing primers, the tool selects the minimal set of primers that provides the maximal coverage of template sequences. Designed primers are of high quality as several molecular properties can be considered. Moreover, users can specify the required binding conditions in detail: template-specific binding regions can be defined and the maximal number of allowed mismatches for primer binding can be selected. Moreover, hybridization and elongation efficiencies of individual primer binding events are considered to determine which templates can actually be amplified. The tool enables the comparison of the properties (e.g. binding regions and coverage) of existing primer sets. 

## Motivation

Existing approaches for multiplex primer design did not fulfill our requirements: many tools are not very customizable (allowed regions, primer properties, mismatch binding, primer directionality), do not optimize the intended target function (primer set is not always minimized, MPPrimer), are not user friendly (DECIPHER), or are not interpretable (no evaluation/comparison of primers). 

## Using openPrimeR

We provide two distributions for openPrimeR. The first way of distribution is the openPrimeR R package with the Shiny app, which can be retrieved via GitHub. 
The second type of distribution is the Shiny app as a self-contained Docker container. The openPrimeR Docker image contains all of the dependencies of the tool, which makes it easily usable on any type of system, independent of the installed operating system, R distribution, etc. We recommend that you use the Docker image if you belong to the following groups of users:

1. You do not have any experience with R whatsoever and you do not intend to use the functionalities of the R package.
2. You only want to use the openPrimeR frontend in terms of the Shiny app.
3. You do not want to invest time in installing the package with all of its dependencies.

## Usage of the tool from GitHub

### Introduction to GitHub
In a console, enter

>**git clone https://github.molgen.mpg.de/mdoering/primer_design**

and then enter your username and password.

To update your local version of the tool to the curreent GitHub version at a later point in time, just run

>**git pull origin master**

in the project's base folder. 

In case that you have performed changes to the local files and you would like to revert to the last GitHub version, execute

>**git checkout -- .**

in the project folder in order to discard all local changes.

### Installing the tool: MacOS and Unix

In a console, enter the project's base folder and execute

>**./install.sh**

### Installing the tool: Windows

If you are using Windows, please execute the batch script

>**./install.cmd**

Please ensure that you have added R into your system's path before using the script. For more information, on changing the path, please refer to the [R FAQ](https://cran.r-project.org/bin/windows/base/rw-FAQ.html#Rcmd-is-not-found-in-my-PATH_0021).h

### Installation advice

The installation command installs all of the package's dependencies, that is, other R packages as well as the third-party tools that are required by openPrimeR. 
If there are problems during the installation, please consider the console output. For example, if you should get an error such as

> there is no package called 'gtable'

this can probably be fixed by running

> install.packages('gtable')

in case of CRAN packages such as *gtable* or. For Bioconductor packages, such as *XVector*, you can run

> biocLite('XVector')

fater loading *biocLite*. After this, you restart the installation script. Once the installation is finished, a browser will open the tool.
Again, please note the console output here for the case that there were problems with installing the third-party tools.

### Starting the tool

If you are using Unix/MacOS, simply execute the bash script

>**./start.sh**

from the project's base folder in a console to start a browser with the app.

If you are using Windows, you can run the batch file

>**start.cmd**

located in the project's base folder.

## Usage with Docker

Our docker image is available at [dockerhub](https://hub.docker.com/r/mdoering88/primer_design/). In order to use the docker image, you need to [install docker](https://www.docker.com/) on your system and activate the docker daemon. 

### Short instructions

If you are using Unix/MacOS, please open a console, enter the project's base directory and run the bash script

>**./start_docker.sh**

This procedure will

1. Download the Docker image if it isn't available yet.
2. Run the Docker image (see details below).
3. Start a web browser for using the app.

### Detailed instructions

After logging in with your docker account, in a console, enter 

>**docker pull mdoering88/primer_design**

to retrieve the latest docker image of the tool. Since the image is quite large the download (~4 GB) may take some time, especially if your internet connection is not very fast.

To run the image, enter

>**docker run -p 3838:3838 --rm mdoering88/primer_design**

After this, the tool is available by accessing **localhost:3838** in your web browser.

In case you want to have more control of the image you are running or you want to study the output of the tool, you can execute

>**docker run --rm -p 3838:3838 -v /tmp/logs/:/var/log/shiny-server/ mdoering88/primer_design:0.99.0**

With this call, the tool's log file is stored in the */tmp/logs/* folder on your system with a filename starting with *shiny-shiny*). During the session, you can use 

>**tail -f shiny-shiny-X.txt**

to retrieve the current status messages. 

Moreover, in the above call to Docker, we have specified a tag, namely *:0.99.0*, which means that we have started version *0.99.0* of the tool of the tool. In case that no tag is provided, the latest available Docker image is used.

## Contributors

openPrimeR is being developed at the Max Planck Institute for Informatics and the University of Cologne.

## License

See the [LICENSE](LICENSE.txt) file for license rights and limitations (GNU General Public License, Version 2.0).

## Changelog

Take a look at the [CHANGELOG](CHANGELOG.md) to view recent changes to the project.

## Requirements for the R package
- R >= version 3.1.0 
- OS: Linux, MacOS, Windows
- Third-party software: MAFFT, oligoarrayaux, MELTING, phantomjs, selenium for python, ViennaRNA

