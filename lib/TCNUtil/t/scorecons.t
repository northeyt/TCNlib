#!/usr/bin/env perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl scorecons.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl scorecons.t'
#   Tom Northey <zcbtfo4@acrm18>     2015/04/22 11:06:42

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More qw( no_plan );
BEGIN { use_ok( 'TCNUtil::scorecons' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

sub test_startup {
    my $test     = shift;
    my $testFile = "test.aln.fasta";
    my $testObj  = new_ok("scorecons", ["inputAlignedSeqsStringFile", $testFile]);

    my @gotScores = $testObj->calcConservationScores();
    is(scalar @gotScores, 209, "all scores parsed");
    is($gotScores[0], 0.149, "score is correct");
    
    $testObj->targetSeqIndex(2);
    $testObj->calcConservationScores();
}
