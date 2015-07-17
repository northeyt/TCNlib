#!/usr/bin/env perl
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

ok(testClustalw(), "MSA::ClustalW ok");
ok(testMuscle(), "MSA:Muscle ok");
ok(testClustalO(), "MSA::ClustalO ok");
ok(testCalcConservationScores(), "testCalcConservationScores ok");
ok(test_getSequenceStringsFromInputFile(), "getSequenceStringsFromInputFile ok");

sub testClustalw {
    my $testInput = "test.fasta";
    my $testObj = new_ok('MSA::Clustalw', [inputSeqsFile => $testInput]);
    
    my $expEle
        = "---------------------------------------------------------------------"
        . "---------------------------------------------------------------------"
        . "---------------------------------------------------------------------"
        . "---------------------------------------------------------------------"
        . "-----------------ELRDKKQKV--------------------HALFYKLDIV-------------"
        . "---------------------------------------------------------------------"
        . "-----------------------------------------------------------------";
    
    is(scalar @{$testObj->getAlignedSeqStringAref}, 32,
       "all seqs present in align");
    is($testObj->getAlignedSeqStringAref->[0], $expEle,
       "align works okay");
}

sub testClustalO {
    my $testInput = "test.fasta";
    my $testObj = new_ok('MSA::ClustalO', [inputSeqsFile => $testInput]);
    
    is(scalar @{$testObj->getAlignedSeqStringAref}, 32,
       "all seqs present in align");
}

sub testMuscle {
    my $testInput = "test.fasta";
    my $muscleTestObj = new_ok('MSA::Muscle', [inputSeqsFile => $testInput]);

    my $expEle
        = "------------------------------------------------------------"
        . "------------------------------------------------------------"
        . "------------------------------------------------------------"
        . "------------------ELRDKKQK----VH----------------------------"
        . "------------------------------------------------------------"
        . "------------------------------------------------------------"
        . "------------------ALF---YKLDIV------------------------------"
        . "-------------------------------------------------";

    is(scalar @{$muscleTestObj->getAlignedSeqStringAref}, 32,
     "all seqs present in align");
    is($muscleTestObj->getAlignedSeqStringAref->[0], $expEle,
       "align works okay");
}

sub testCalcConservationScores {
    my $testInput = "test.fasta";
    use scorecons;
    my $testObj = new_ok('MSA::Clustalw', [inputSeqsFile => $testInput]);

    $testObj->consScoreCalculator(scorecons->new());
    ok($testObj->calculateConsScores(), "calculateConsScores ok");    
}

sub test_getSequenceStringsFromInputFile {
    my $testInput = "test.fasta";
    my $testObj = new_ok('MSA::Clustalw', [inputSeqsFile => $testInput]);

    my @seqStrings = $testObj->getSequenceStringsFromInputFile();
    is(scalar @seqStrings, 32, "all sequences from input file returned");

    is($seqStrings[0], 'ELRDKKQKVHALFYKLDIV', "first seq is correct");

    my $expLastSeq = 'MKHHHHHHHHHHSSDYKDDDDKGENLYFQGSKIEEGKLVIWINGDKGYNGLAEVGKK'
        . 'FEKDTGIKVTVEHPDKLEEKFPQVAATGDGPDIIFWAHDRFGGYAQSGLLAEITPDKAFQDKLYPFTW'
        . 'DAVRYNGKLIAYPIAVEALSLIYNKDLLPNPPKTWEEIPALDKELKAKGKSALMFNLQEPYFTWPLIA'
        . 'ADGGYAFKYENGKYDIKDVGVDNAGAKAGLTFLVDLIKNKHMNADTDYSIAEAAFNKGETAMTINGPW'
        . 'AWSNIDTSKVNYGVTVLPTFKGQPSKPFVGVLSAGINAASPNKELAKEFLENYLLTDEGLEAVNKDKP'
        . 'LGAVALKSYEEELAKDPRIAATMENAQKGEIMPNIPQMSAFWYAVRTAVINAASGRQTVDEALKDAQTN';
    
    is($seqStrings[31], $expLastSeq, "last seq is correct");
}
