#!/usr/bin/perl -Iblib/lib -Iblib/arch -I../blib/lib -I../blib/arch
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl chaincontacts.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl chaincontacts.t'
#   Tom Northey <zcbtfo4@acrm18>     2014/06/10 17:20:16

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use Test::More qw( no_plan );
use Test::Deep;
BEGIN { use_ok( 'pdb::chaincontacts' ); }
use pdb::pdb;

### Tests ######################################################################

my $testPDB = getTestPDB();

my $testObj = pdb::chaincontacts->new(input => $testPDB, threshold => 3);

cmp_deeply([pdb::chaincontacts::_parseLine(getTestLine())], [getExpLineArray()],
   "_parseLine works okay");

my $testResObj = $testObj->getOutput(); 
is(ref $testResObj, "pdb::chaincontacts::result", "getOutput ok");

cmp_bag($testResObj->chain2chainContacts(['A'], ['B']),
           getExpContacts(), "result::chain2chainContacts works okay");


### Subroutines ################################################################

sub getExpContacts {
    return [qw(B.17 B.58)];
}


sub getExpLineArray {
    return(qw(K 103 B 76 1));
}

sub getTestLine {
    return "Chain: K Res: 103  - Chain: B Res:  76  Contacts: 1";
}

# Returns array with H, L and antigen chain
sub getTestPDB {
    my $pdb_code = "1afv";
    my $pdb_file = $pdb_code . ".pdb";

    my $pdb = pdb->new(pdb_code => $pdb_code,
                       pdb_file => $pdb_file);

    return $pdb;
}

sub getTestChains {
    my $pdb = shift;

    my @chains = $pdb->create_chains('H', 'L', 'A');

    return(@chains);
}
