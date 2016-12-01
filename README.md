# About

This github is home to a bunch of perl modules and scripts that were written during the course of my PhD. Most modules are written using Moose, the perl library for object-orientation. The main bulk of the modules are used to manipulate PDB files and in particular to send or process data to and from [bioptools](https://github.com/ACRMGroup/bioptools) programs.

# Prerequisites

These prerequisites are required for TCNUtils

* Get submodules
TCNUtils is currently dependent on a handful of modules in the [SAAP package](https://github.com/ACRMGroup/SAAP). SAAP is a submodule of TCNUtils, so can be obtained by running `git submodule update --init --recursive` on the command line.

* Install bioptools
Available via [github](https://github.com/ACRMGroup/bioptools).

* Install BioPerl
See the [BioPerl website](http://bioperl.org/INSTALL.html) for instructions

* Install Module::Builder
This module is required to check for and install perl dependecies. Run `cpan Module::Builder` on the command-line to install it.

# Setup

* Install external packages
In order to obtain the non-perl programs run by TCNutils, simply run `getexternalpackages`. Note that this will download the linux versions of these programs (other operating systems are currently unsupported).

* Install perl dependencies
Run `getperldeps.pl` to check for perl dependencies - any missing dependencies will be installed.

* Run tests
Run `runtests.pl`.

* Create `TCNPerlVars.pm`
Defaults for variables used by all TCNUtils modules are found in `TCNPerlVars.defaults`. This has to be copied to `TCNPerlVars.pm` for TCNUtils to work. At the command-line, run `cp lib/TCNPerlVars.defaults lib/TCNPerlVars.pm`.

* Get swissprot for blast searches (not necessary, but recommended!)
Running blast searches remotely is implemented in this library, but the NCBI cgi server can be unreliable and is always slow. To run blast searches against a local version of the swissprot database, run setup-blastdb (found in /bin). This will grab the latest version of swissprot from the NCBI FTP server and create the required database files so that blastall can be run. The files will be created at the location specified by the $TCNPerlVars::blastdbDir variable, found in TCNPerlVars.pm. Alternatively, if you already have a local version of the database, point the $TCNPerlVars::swissProtDB towards it.
