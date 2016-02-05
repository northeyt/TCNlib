#!/usr/bin/env perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl pdbsws.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl pdbsws.t'
#   Tom Northey <zcbtfo4@acrm18>     2013/11/11 10:32:01

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More qw( no_plan );
use Test::Deep;
BEGIN { use_ok( 'pdb::pdbsws' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.


subtest "test pdbsws locally" => sub {
    my $pdbsws = new_ok('pdb::pdbsws::Local');
    testpdbsws($pdbsws);    
};

subtest "test pdbsws remote" => sub {
    my $pdbsws = new_ok('pdb::pdbsws::Remote');
    testpdbsws($pdbsws);
};

sub testpdbsws {
    my $pdbsws = shift;
    
    cmp_deeply([$pdbsws->getACsFromPDBCodeAndChainID('4hou', 'A') ],
               ['P09914'], " getACsFromPDBCodeAndChainID returns correct ACs");
    
    cmp_deeply([$pdbsws->getIDsFromPDBCodeAndChainID('4hou', 'A') ],
               ['IFIT1_HUMAN'], " getIDsFromPDBCodeAndChainID returns correct IDs");
    
    testMapResSeq2SwissProtNum($pdbsws);
}

sub testMapResSeq2SwissProtNum {
    my $pdbsws = shift;
    my %gotMap = $pdbsws->mapResSeq2SwissProtNum("4hou", 'A', 'P09914');
    my $expMapSize = 255;
    
    is (scalar keys %gotMap, $expMapSize, "expected number of residues mapped");
    is($gotMap{10},  10,  "first chain residue is mapped correctly");
    is($gotMap{278}, 278, "last chain residue is mapped correctly");    
}
