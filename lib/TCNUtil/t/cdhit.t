#!/usr/bin/env perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl cdhit.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl cdhit.t'
#   Tom Northey <zcbtfo4@acrm18>     2014/06/30 14:54:15

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;
use Test::More qw( no_plan );
BEGIN { use_ok( 'TCNUtil::cdhit' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my $testObj = new_ok(cdhit => [input => 'test.fasta']);
ok($testObj->getClusters(), "getClusters() works okay");

ok(testBuildWordLength($testObj), "_buildWordLength works ok");

sub testBuildWordLength {
    my $testObj = shift;
    my %seqID2expWordLengths = ("0.5" => 2, "0.6" => 3, "0.7" => 4, "0.9" => 5);

    while (my($seqID, $wLen) = each %seqID2expWordLengths) {
        $testObj->seqIDThreshold($seqID);
        return 0 if $testObj->_buildWordLength() ne $wLen;
    }
    return 1;
}

