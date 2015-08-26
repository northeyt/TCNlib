#!/usr/bin/env perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl asurf64.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl asurf64.t'
#   Tom Northey <zcbtfo4@acrm18>     2014/07/10 14:24:11

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;
use Test::More qw( no_plan );
use Test::Deep;
use File::Basename;

BEGIN { use_ok( 'pdb::solv' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my $testObj = new_ok('pdb::solv', [input => "1djs.pdb"]);

ok(my $atomASAHref = $testObj->getOutput(), "getOutput works ok");

cmp_deeply($atomASAHref->{1}, 44.61,
           "atom ASA value is correct");

is($testObj->resid2RelASAHref->{"A.147"}->{allAtoms}, "118.003",
   "getResid2RelASAHash works ok")
