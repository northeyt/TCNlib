#!/usr/bin/env perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl FOSTA.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl FOSTA.t'
#   Tom Northey <zcbtfo4@acrm18>     2015/07/28 11:38:17

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;
use Test::More qw( no_plan );
use Test::Deep;
use Test::Exception;
use Getopt::Long;
use Carp;

BEGIN { use_ok( 'TCNUtil::FOSTA' ); }
use FindBin qw($RealBin);
chdir($RealBin); # So test data files can be found

my $runLocalTests = 0;
GetOptions("l" => \$runLocalTests);
my $skipMessage = "Requires local db (supply -l opt to run)";

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

subtest "test getFOSTAFamIDAndReliability" => sub {
    plan skip_all => $skipMessage unless $runLocalTests;
    my $testObj = FOSTA::Factory->new(remote => 0)->getFOSTA();
    my $testReliableID = "CNTD1_HUMAN";
    cmp_deeply([$testObj->getFOSTAFamIDAndReliability($testReliableID)], [1, 13],
               "getFOSTAFamIDAndReliability returns reliable family id");
    
    my $testUnRelID = "CF211_XENTR";
    cmp_deeply([$testObj->getFOSTAFamIDAndReliability($testUnRelID)], [0, 18591],
               "getFOSTAFamIDAndReliability returns unreliable family id");
};

subtest "test getFEPIDsFromFamID" => sub {
    plan skip_all => $skipMessage unless $runLocalTests;
    my $testObj = FOSTA::Factory->new(remote => 0)->getFOSTA();
    my $famID = "CNTD1_HUMAN";
    cmp_deeply([$testObj->getFEPIDsFromFamID($famID, 13)], ["CNTD1_MOUSE"],
               "getFEPIDsFromFamID returns correct FEPIDs");
};

subtest "test getSequenceFromID" => sub {
    plan skip_all => $skipMessage unless $runLocalTests;
    my $testObj = FOSTA::Factory->new(remote => 0)->getFOSTA();
    my $famID = "CNTD1_HUMAN";
    my $gotSeq  = $testObj->getSequenceFromID($famID);
    is($gotSeq->string(), expSeqStr(), "getSequenceFromID retuns sequence ok");
};

subtest "test getReliableFEPSequencesFromSwissProtID" => sub {
    plan skip_all => $skipMessage unless $runLocalTests;
    my $testObj = FOSTA::Factory->new(remote => 0)->getFOSTA();
    my $famID = "CNTD1_HUMAN";    
    cmp_deeply([$testObj->getReliableFEPSequencesFromSwissProtID($famID)],
               array_each(isa("sequence")),
               "getReliableFEPSequencesFromSwissProtID returns sequences");
};

subtest "test FOSTA::Remote" => sub {
    my $testObj = FOSTA::Factory->new(remote => 1)->getFOSTA();
    my $famID   = "CNTD1_HUMAN";
    cmp_deeply([$testObj->getReliableFEPSequencesFromSwissProtID($famID)],
               array_each(isa("sequence")),
               "run remotely, getReliableFEPSequencesFromSwissProtID returns sequences");
};

sub expSeqStr {
    my $seq = <<EOF;
MDGPMRPRSASLVDFQFGVVATETIEDALLHLAQQNEQAVREASGRLGRFREPQIVEFVFLLSEQWCLEKSVSYQAVEIL
ERFMVKQAENICRQATIQPRDNKRESQNWRALKQQLVNKFTLRLVSCVQLASKLSFRNKIISNITVLNFLQALGYLHTKE
ELLESELDVLKSLNFRINLPTPLAYVETLLEVLGYNGCLVPAMRLHATCLTLLDLVYLLHEPIYESLLRASIENSTPSQL
QGEKFTSVKEDFMLLAVGIIAASAFIQNHECWSQVVGHLQSITGIALASIAEFSYAILTHGVGANTPGRQQSIPPHLAAR
ALKTVASSNT
EOF

    $seq =~ s/\s//gxms;
    return $seq;
}

1;
