#!/usr/bin/env perl
use strict;
use warnings;
use pdb;
use pdb::pdbFunctions;
use pdb::get_files;
use Getopt::Long;

my $inputIsPqsCode = 0;
GetOptions("p" => \$inputIsPqsCode);

my $fileType = $inputIsPqsCode ? "pqs_file" : "pdb_file";

while (<>) {
    chomp $_;
    my $chainID = chop $_;
    my $pdbCode = $_;
    my $pdbFile = pdb::get_files->new(pdb_code => $pdbCode)->$fileType();
    my $chain
        = chain->new(pdb_code => $pdbCode, chain_id => $chainID,
                     solvent_cleanup => 1, het_atom_cleanup => 1,
                     pdb_file => $pdbFile
                 );
    $chain->read_ASA();
    print join(":", ($chain->pdb_code, $chain->chain_id, substr($_, 2),
                     $chain->resid2RelASAHref->{$_}->{allAtoms})) . "\n"
                         foreach sort {pdb::pdbFunctions::compare_resids($a, $b)} $chain->getResIDs(); 
}
