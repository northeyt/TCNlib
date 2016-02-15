#!/usr/bin/env perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl kabatnum.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl kabatnum.t'
#   Tom Northey <zcbtfo4@acrm18>     2014/06/07 13:29:37

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use lib ("../..");
use strict;
use warnings;
use Test::More qw( no_plan );
use pdb;
BEGIN { use_ok( 'pdb::kabatnum' ); }
use FindBin qw($RealBin);
chdir($RealBin); # So test data files can be found

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my $testObj = pdb::kabatnum->new(input => "1afvH.pdb");         
$testObj->getOutput();

$testObj->input(getTestChain());
$testObj->sequenceChain();

is($testObj->input->resid_index->{'H.52'}->{'C'}->is_CDR, 1);

sub getTestChain {
    my $chain = chain->new(pdb_file => "1afv.pdb",
                           pdb_code => "1afv",
                           chain_id => "H");

    return $chain;
}

