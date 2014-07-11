#!/usr/bin/perl -Iblib/lib -Iblib/arch -I../blib/lib -I../blib/arch
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl pdbFunctions.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl pdbFunctions.t'
#   Tom Northey <zcbtfo4@acrm18>     2014/06/24 17:07:15

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;
use Test::More qw( no_plan );
BEGIN { use_ok( 'pdb::pdbFunctions' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.


# test compareResSeqs
my $rsA = "52";
my $rsB = "52A";

is(pdb::pdbFunctions::compare_resSeqs($rsA, $rsB), -1,
   "compareResSeqs works okay");

my $riA = "C121";
my $riB = "C27";

is(pdb::pdbFunctions::compare_resids($riA, $riB), 1,
   "compare_resids works okay");
