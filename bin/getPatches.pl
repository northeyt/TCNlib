#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Carp;
use pdb;
use pdb::get_files;
use pdb::automatic_patches;

my $fileType = "pdb";
GetOptions("f=s", \$fileType);

@ARGV or Usage() and exit(1);
my $pdbCode = shift @ARGV;
my $chainID = shift @ARGV;

my $pdbGetFile = pdb::get_files->new();
$pdbGetFile->pdb_code($pdbCode);
my $file =  $fileType eq 'pdb' ? $pdbGetFile->pdb_file()
    :  $fileType eq 'pqs' ? $pdbGetFile->pqs_file()
    :  croak "Input from file type " . $fileType . " not implemented!";

my $chain = chain->new(pdb_code => $pdbCode, chain_id => $chainID,
                       pdb_file => $file);
my $ap = automatic_patches->new(pdb_object => $chain, patch_type => 'normal',
                                radius => 14, ASA_type => 'ASAb',
                                build_patches_from_parent => 1);

print map {$_->summary()} $ap->get_patches();

sub Usage {
    print <<EOF;
$0 -f [pdb|pqs] PDB_CODE CHAIN_ID
EOF
}
