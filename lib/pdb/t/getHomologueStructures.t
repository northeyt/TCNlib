#!/usr/bin/env perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl getHomologueStructures.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl getHomologueStructures.t'
#   Tom Northey <zcbtfo4@acrm18>     2015/02/11 15:04:06

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More qw( no_plan );
use Test::Deep;
use Test::Exception;

use pdb::pdb;

BEGIN { use_ok( pdb::getHomologueStructures ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

test_getHomologueStructures();

sub test_getHomologueStructures { 
    
    my $chain = chain->new(pdb_code => "1za7",
                           pdb_file => "1za7.pdb",
                           chain_id => "A");

    my $testObj
        = pdb::getHomologueStructures->new(query => $chain,
                                           db => "testDB/pdbaa");

    can_ok($testObj, "runBlastall");
           
    my @hits = $testObj->getHits();
  
    is(scalar @hits, 3, "getHits - correct number of hits returned");
    
    my $tests = all(
        isa("Bio::Search::Hit::BlastHit"),
    );
    
    cmp_deeply(\@hits, array_each($tests),
               "getHits: hits are Bio::Search::Hit::BlastHit objects");
    
    my @expHitNames = qw(pdb1CWPA refNP_041199.1 pdb1YC6A);
    my @gotHitNames = map {join("",$testObj->_parseHitName($_))} @hits;

    cmp_deeply(\@gotHitNames, \@expHitNames, "_parseHitName works ok");

    my $hitChain = $testObj->getHitStructure($hits[0]);
    is($hitChain->pdbID(), "1CWPA", "getHitStructure works ok");

    my $badHit = $hits[1];
    dies_ok { $testObj->getHitStructure($badHit) }
        'getHitStructure dies when hit is not a pdb chain sequence';

    my %gotAlignMap = $testObj->getAlignMap($hits[0], $hitChain);

    cmp_deeply(\%gotAlignMap, {expAlignMap()}, "getAlignMap ok");

}

sub expAlignMap {
    my @query = (1  .. 140);
    
    my @hit   = (26 .. 165);
    
    my %hash = ();
    @hash{@query} = @hit;
   
    return %hash;
}
