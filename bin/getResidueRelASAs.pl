#!/usr/bin/env perl
use strict;
use warnings;
use pdb;
use pdb::pdbFunctions;

while (<>) {
    my $chain = chain->new(pdb_code => substr($_, 0, 4),
                           chain_id => substr($_, 4, 1),
                           solvent_cleanup => 1,);
    $chain->read_ASA();
    print join(":", ($chain->pdb_code, $chain->chain_id, substr($_, 2),
                     $chain->resid2RelASAHref->{$_}->{allAtoms})) . "\n"
                         foreach sort {pdb::pdbFunctions::compare_resids($a, $b)} $chain->getResIDs(); 
}
