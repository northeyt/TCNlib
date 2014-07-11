#!/acrm/usr/local/bin/perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl MSA.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl MSA.t'
#   Tom Northey <zcbtfo4@acrm18>     2014/07/01 14:41:24

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;
use Test::More qw( no_plan );
BEGIN { use_ok( 'MSA' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my $testInput = "test.fasta";
my $testObj = new_ok('MSA', [input => $testInput]);

my $expStr = ">4hpyP\nELRDKKQKVHALFYKLDIV\n\n";
is($testObj->processedInput()->[0], $expStr, "_processInput works okay");

my $expEle
    = "---------------------------------------------------------------------"
    . "---------------------------------------------------------------------"
    . "---------------------------------------------------------------------"
    . "---------------------------------------------------------------------"
    . "-----------------ELRDKKQKV--------------------HALFYKLDIV-------------"
    . "---------------------------------------------------------------------"
    . "-----------------------------------------------------------------";

is([$testObj->align()]->[0], $expEle, "align works okay");

my $expSeqStr = "ELRDKKQKVHALFYKLDIV";
is(MSA::seqStrFromFASTAStr($expStr), $expSeqStr, "seqStrFromFASTAStr ok");
