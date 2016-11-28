#!/usr/bin/env perl
use strict;
use warnings;
use pdb;

@ARGV > 1 or Usage();
my $pdbCode = shift @ARGV;
my @resIDs = @ARGV;

my $pdb = pdb->new(pdb_code => $pdbCode);

foreach my $resID (@resIDs) {
    my @atoms = values %{$pdb->resid_index->{$resID}};
    my $resName  = $atoms[0]->resName();
    print "residue $resID: resName = $resName\n";
}

sub Usage {
    print "Please supply a pdb code and one or more resIDs\n";
    exit(1);
}
