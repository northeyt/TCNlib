# About

This github is home to a bunch of perl modules and scripts that were written during the course of my PhD. Most modules are written using Moose, the perl library for object-orientation. The main bulk of the modules are used to manipulate PDB files and in particular to send or process data to and from [bioptools](https://github.com/ACRMGroup/bioptools) programs.

# Setup

These prerequisites are required for TCNlib

1. Install bioptools, available via [github](https://github.com/ACRMGroup/bioptools).

2. Install BioPerl. See the [BioPerl website](http://bioperl.org/INSTALL.html) for instructions. I recommend opening the cpan client on the command line and typing `d /bioperl/`. This will return a list of available BioPerl packages. Pick the most recent version and use force to install (it is likely that some tests will fail) by running `force install CJFIELDS/BioPerl-1.6.924.tar.gz`.

3. Install Module::Build, which is required to check for and install perl dependecies. Run `cpan Module::Build` on the command-line to install it.

4. Run setup.pl. This installs git submodules and external packages, checks for perl dependencies (installing any that are missing) and makes a working copy of lib/TCNPerlVars.pm from lib/TCNPerlVars.defaults.

5. Run tests using `runtests.pl`.

6. Add the following lines to your .bashrc:

```
export TCNlib=/path/to/TCNlib
export PERL5LIB=$TCNlib/lib:$PERL5LIB
```

Where `/path/to/TCNlib/` is replaced with the real path!

7. Get swissprot for blast searches (not necessary, but recommended!). Running blast searches remotely is implemented in this library, but the NCBI cgi server can be unreliable and is always slow. To run blast searches against a local version of the swissprot database, run setup-blastdb. This will grab the latest version of swissprot from the NCBI FTP server and create the required database files so that blastall can be run. The files will be created at the location specified by the $TCNPerlVars::blastdbDir variable, found in TCNPerlVars.pm. Alternatively, if you already have a local version of the database, point the $TCNPerlVars::swissProtDB towards it.