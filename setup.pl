#!/usr/bin/env perl
use strict;
use warnings;

print "Initialising sub-modules ...\n";
system("git submodule update --init --recursive");

print "Getting external packages ...\n";
system("./getexternalpackages");

print "Getting perl dependencies ...\n";
system("./getperldeps.pl");

print "Copying TCNPerlVars.defaults to TCNPerlVars.pm ...\n";
system("cp lib/TCNPerlVars.defaults lib/TCNPerlVars.pm");
