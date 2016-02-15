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
BEGIN { use_ok( 'TCNUtil::MSA' ); }
use FindBin qw($RealBin);
chdir($RealBin); # So test data files can be found

use TCNUtil::sequence;
use TCNUtil::scorecons;

use File::Basename;
my ($testFile, $testDir, $suffix) = fileparse($0);
chdir($testDir);

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

subtest "testClustalw" => sub  {
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
};

subtest "testClustalO" => sub {
    my $testInput = "test.fasta";
    my $testObj = new_ok('MSA::ClustalO', [inputSeqsFile => $testInput]);
    
    is(scalar @{$testObj->getAlignedSeqStringAref}, 32,
       "all seqs present in align");
};

subtest "run muscle locally" => sub {
    my $testInput = "test.fasta";
    my $muscleTestObj = new_ok('MSA::Muscle::Local', [inputSeqsFile => $testInput]);
    _testMuscleAlignedSeqStringAref($muscleTestObj->getAlignedSeqStringAref());
};

subtest "run muscle remotely" => sub {
    my $muscleTestObj = new_ok('MSA::Muscle::Remote', [seqs => [_getSeqs()]]);
    _testMuscleAlignedSeqStringAref($muscleTestObj->getAlignedSeqStringAref());
};

subtest "testCalcConservationScores" => sub {
    my $testInput = "test.fasta";
    my $testObj = new_ok('MSA::Clustalw', [inputSeqsFile => $testInput]);

    $testObj->consScoreCalculator(scorecons->new());
    ok($testObj->calculateConsScores(), "calculateConsScores ok");
};

subtest "test_getSequenceStringsFromInputFile" => sub {
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
};

sub _testMuscleAlignedSeqStringAref {
    my $seqStrAref = shift;

    is(scalar @{$seqStrAref}, 32, "all seqs present in alignment");
    
    my $expSeq        = "ELRDKKQKVHALFYKLDIV";
    my $gotAlignedSeq = $seqStrAref->[0]; 
    $gotAlignedSeq    =~ s/-//g;
    is($gotAlignedSeq, $expSeq, "first sequence is fully aligned");
}

sub _getSeqs {
    my $inFile = "test.fasta";
    open(my $IN, "<", $inFile) or die "Cannot open file $inFile, $!";
    my $data;
    {
        local $/;
        $data = <$IN>;
    }
    my @FASTASeqs = split(/\n\n/, $data);
    return map {sequence->new($_)} @FASTASeqs;
}
