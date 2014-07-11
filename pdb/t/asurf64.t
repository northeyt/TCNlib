#!/usr/bin/perl -Iblib/lib -Iblib/arch -I../blib/lib -I../blib/arch
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
BEGIN { use_ok( 'pdb::asurf64' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my $testObj = new_ok('pdb::asurf64', [input => "1djs.pdb"]);

is($testObj->runExec(), "/tmp/1djs.asa", "runExec works ok");

my $line
    = "ATOM      1  N   THR A 147      18.358  47.738  23.896  43.469  1.60\n";
my @exp = qw(1 43.469 1.60);

cmp_deeply([$testObj->parseLine($line)], \@exp, "parseLine works okay");

ok($testObj->getOutput(), "getOutput works ok);
