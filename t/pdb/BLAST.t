#!/usr/bin/env perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl runBlast.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl runBlast.t'
#   Tom Northey <zcbtfo4@acrm18>     2015/07/23 14:05:53

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use Test::More qw( no_plan );
use Test::Deep;
use Test::Exception;
use Bio::SearchIO::blast;
use IO::CaptureOutput qw(capture);
use Test::MockObject::Extends;
use pdb;

BEGIN { use_ok( 'pdb::BLAST' ); }
use FindBin qw($RealBin);
chdir($RealBin); # So test data files can be found

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

subtest "BLAST local" => sub {
    $ENV{'BLASTDB'} = "testDB";
    my $testObj = pdb::BLAST::Local->new(query         => _getTestChain(),
                                         db            => "1za7",
                                         reportHandler => pdb::BLAST::Report::PDBseq->new());
    can_ok($testObj, "runBlast");
    my $gotBlastReport = $testObj->runBlast();
    cmp_deeply($gotBlastReport, isa("Bio::SearchIO::blast"),
               "runBlast returns a blast report");
};

subtest "BLAST remote" => sub {
    my $testObj = pdb::BLAST::Remote->new(query         => _getTestChain(),
                                          db            => "swissprot",
                                          reportHandler => pdb::BLAST::Report::SwissProt->new());
    can_ok($testObj, "runBlast");
    my $gotBlastReport = $testObj->runBlast();
    cmp_deeply($gotBlastReport, isa("Bio::SearchIO::blast"),
               "runBlast returns a blast report");
};

subtest "BLAST Report" => sub {
    my $testObj = pdb::BLAST::Report::PDBseq->new(query => _getTestChain(),
                                                  report => _getMockBlastReport());
    my @hits = $testObj->getHits();
    
    is(scalar @hits, 3, "getHits - correct number of hits returned");
    cmp_deeply(\@hits, array_each(isa("Bio::Search::Hit::BlastHit")),
               "getHits: hits are Bio::Search::Hit::BlastHit objects");
};

subtest "BLAST PDBseq" => sub {
    my $chain = chain->new(pdb_code => "1za7",
                           pdb_file => "1za7.pdb",
                           chain_id => "A");
    
    my $testObj = pdb::BLAST::Factory->new(dbType => 'pdb')->getBlaster();
    $testObj->setQuery($chain);
    
    can_ok($testObj, "runBlast");
    $testObj->runBlast();
    
    my @hits = $testObj->reportHandler->getHits();
    
    my @expHitNames
        = qw(pdb1CWPA refNP_041199.1 pdb1YC6A);
    my @gotHitNames
        = map {join("",$testObj->reportHandler->_parseHitName($_))} @hits;

    cmp_deeply(\@gotHitNames, \@expHitNames, "_parseHitName works ok");

    my $hitChain = $testObj->reportHandler->getHitStructure($hits[0]);
    is($hitChain->pdbID(), "1CWPA", "getHitStructure works ok");

    my $badHit = $hits[1];
    dies_ok { $testObj->reportHandler->getHitStructure($badHit) }
        'getHitStructure dies when hit is not a pdb chain sequence';

    my %gotAlignMap = $testObj->reportHandler->getAlignMap($hits[0], $hitChain);

    cmp_deeply(\%gotAlignMap, {expAlignMap()}, "getAlignMap ok");

};

subtest "BLAST SwissProt" => sub {
    
    my $testChain = chain->new(pdb_code => "1afv", pdb_file => "1afv.pdb",
                               chain_id => "A");
    my $testObj
        = pdb::BLAST::Factory->new(dbType => 'swsprot')->getBlaster();
    $testObj->setQuery($testChain);

    $testObj->runBlast();
    my @hits = $testObj->reportHandler->getHits(reliable => 1);

    my $testAC   = "P12345";
    my $testName = "sp|$testAC|GAG_HV1N5";
    $hits[0]->name($testName);

    is($testObj->reportHandler->parseACFromHit($hits[0]), $testAC,
       "parseACFromHitName works ok");

    like($testObj->reportHandler->swissProtSeqFromHit($hits[0]), qr/[A-Z]+/,
         "swissProtSeqFromHit ok");
    
};


### SUBROUTINES ################################################################
################################################################################

sub _getTestChain {
    return chain->new(pdb_code => "1za7",
                      pdb_file => "1za7.pdb",
                      chain_id => "A");
}

sub _getMockBlastReport {
    my $report = Bio::SearchIO::blast->new();
    $report    = Test::MockObject::Extends->new($report);
    $report->mock("next_result",
                  sub { return _getOneMockBlastResultThenUndef() });
    return $report;
}

# Closure to keep track of number of blast result objects returned
{
    my $count = 0;
    sub _getOneMockBlastResultThenUndef {
        ++$count;
        if ($count > 1) {
            return undef;
        }
        else {
            # capture supresses warnings sent to STDERR
            my $result = capture sub{ Bio::Search::Result::BlastResult->new() };
            $result    = Test::MockObject::Extends->new($result);
            $result->mock("next_hit",
                          \&_getThreeBlastHitsThenUndef);
            return $result;
        }
    }
}

# Closure to keep track of number of blast hit objects returned
{
    my $count = 0;
    sub _getThreeBlastHitsThenUndef {
        ++$count;
        return $count > 3 ? undef : _getBlastHit();
    }
}

sub _getBlastHit {
    # capture supresses warnings sent to STDERR
    my $hit
        = capture sub { Bio::Search::Hit::BlastHit->new(name => "mockID",
                                                        frac_idendtical => 0.5) };
    return $hit;
}

sub expAlignMap {
    my @query = (1  .. 140);
    
    my @hit   = (26 .. 165);
    
    my %hash = ();
    @hash{@query} = @hit;
   
    return %hash;
}
