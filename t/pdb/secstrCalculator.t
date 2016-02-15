#!/usr/bin/env perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl secstr.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl secstr.t'
#   Tom Northey <zcbtfo4@acrm18>     2015/07/22 16:13:12

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;
use Test::More qw( no_plan );
BEGIN { use_ok( 'pdb::secstrCalculator' ); }
use FindBin qw($RealBin);
chdir($RealBin); # So test data files can be found

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

new_ok('pdb::secstrCalculator');

subtest "getResID2secStructHref" => sub {
    my $testFile = "1qok.pdb";
    my %gotHash  = pdb::secstrCalculator->new(input => $testFile)->getResID2secStructHref();

    ok(exists $gotHash{"A.27"},
       "first residue present with correctly formatted resID");
    is($gotHash{"A.27"}, "C", "first residue key has correct secstr value");

    ok(exists $gotHash{"A.267"},
       "last residue present with correctly formatted resID");
    is($gotHash{"A.267"}, "e", "last residue key has correct secstr value");

    is(scalar keys %gotHash, 227, "and hash has correct numbers of keys");
};
