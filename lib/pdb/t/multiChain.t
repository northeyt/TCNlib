#!/acrm/usr/local/bin/perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl multiChain.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl multiChain.t'
#   Tom Northey <zcbtfo4@acrm18>     2014/07/10 11:52:21

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;
use Test::More qw( no_plan );
use pdb::pdb;

BEGIN { use_ok( 'pdb::multiChain' ); }

#########################

my @complexChains = getComplexChains();

pdb::multiChain::readASAb(\@complexChains);

use Data::Dumper;

is($complexChains[0]->resid2RelASAHref->{"L.130"}->{allAtoms},
   84.2,
   "resid2RelASAHref assigned correctly to chain");

my $testAtom = $complexChains[0]->atom_array()->[0];
is($testAtom->ASAb(), 42.473, "readASAb works ok");
sub getComplexChains {
    my $pdb = pdb->new(pdb_file => "1afv.pdb",
                            pdb_code => "1afv");

    my($L, $H) = $pdb->create_chains(qw(L H));
    return($L, $H);
}

my $matchingComplexes = [ [$complexChains[0], $complexChains[1]],
                          [$complexChains[1], $complexChains[0]] ];
is(pdb::multiChain::areComplexesIdentical($matchingComplexes), 1,
   "areComplexesIdentical true for identical complexes");
                         
my $misMatchingComplexes = [ [$complexChains[0], $complexChains[1]],
                             [$complexChains[0]] ];

is(pdb::multiChain::areComplexesIdentical($misMatchingComplexes), 0,
   "areComplexesIdentical false for non-identical complexes");
