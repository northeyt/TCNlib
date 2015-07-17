#!/usr/bin/env perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl blastSwissProt.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl blastSwissProt.t'
#   Tom Northey <zcbtfo4@acrm18>     2015/04/16 15:55:51

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More qw( no_plan );
use pdb::pdb;
BEGIN { use_ok( 'pdb::blastSwissProt' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my $testChain = chain->new(pdb_code => "1afv", pdb_file => "1afv.pdb",
                           chain_id => "A");
my $testObj = pdb::blastSwissProt->new(query => $testChain);

my @hits = $testObj->getHits(reliable => 1);

my $testAC   = "P12345";
my $testName = "sp|$testAC|GAG_HV1N5";
$hits[0]->name($testName);

is($testObj->parseACFromHit($hits[0]), $testAC,
   "parseACFromHitName works ok");

like($testObj->swissProtSeqFromHit($hits[0]), qr/[A-Z]+/,
     "swissProtSeqFromHit ok");
