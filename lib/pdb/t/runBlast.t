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
use Bio::SearchIO::blast;
use IO::CaptureOutput qw(capture);
use Test::MockObject::Extends;

BEGIN { use_ok( 'pdb::runBlast' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

subtest "runBlast local" => sub {
    my $testObj = pdb::runBlastLocal->new(query => _getTestChain(),
                                          db => "testDB/pdbaa");
    can_ok($testObj, "runBlast");
    my $gotBlastReport = $testObj->runBlast();
    cmp_deeply($gotBlastReport, isa("Bio::SearchIO::blast"),
               "runBlast returns a blast report");
};


subtest "runBlast remote" => sub {
    my $testObj = pdb::runBlastRemote->new(query => _getTestChain(),
                                           db => "swissprot");
    can_ok($testObj, "runBlast");
    my $gotBlastReport = $testObj->runBlast();
    cmp_deeply($gotBlastReport, isa("Bio::SearchIO::blast"),
               "runBlast returns a blast report");
};

subtest "getHits" => sub {
    my $testObj = pdb::runBlastLocal->new(query => _getTestChain(),
                                          db => "testDB/pdbaa");    
    my @hits = $testObj->getHits(report => _getMockBlastReport());
    
    is(scalar @hits, 3, "getHits - correct number of hits returned");
    cmp_deeply(\@hits, array_each(isa("Bio::Search::Hit::BlastHit")),
               "getHits: hits are Bio::Search::Hit::BlastHit objects");
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
