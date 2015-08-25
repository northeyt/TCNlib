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

subtest "test scorecons" => sub {
    my $test     = shift;
    my $testFile = "test.aln.pir";
    my $testObj  = new_ok("scorecons", [inputAlignedSeqsStringFile => $testFile]);

    my @gotScores = $testObj->calcConservationScores();
    is(scalar @gotScores, 209, "all scores parsed");
    is(sprintf("%.2f", $gotScores[0]), 0.15, "score is correct");
};

subtest "test target seq index" => sub {
    my $test     = shift;
    my $testFile = "test.aln.pir";
    my $testObj  = new_ok("scorecons", [inputAlignedSeqsStringFile => $testFile,
                                        targetSeqIndex             => 0]);

    my @gotScores = $testObj->calcConservationScores();
    is(scalar @gotScores, 189, "only scores that map to first sequence are parsed");
    is(sprintf("%.2f", $gotScores[-1]), 0.220, "final score is correct");
};
