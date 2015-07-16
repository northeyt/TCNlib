################################ CLASS #########################################
###################### pdb::getHomologueStructures #############################
################################################################################
package pdb::getHomologueStructures;

=head1 NAME

pdb::getHomologueStructures - do a BLAST search with a query chain against the
PDB.

=cut

=head1 SYNOPSIS

TODO

=cut


=head1 DESCRIPTION

pdb::getHomologueStructures allows you to use a chain object as a query for a
BLAST search against a database of sequences in the PDB. This package can return
the chain objects that correspond to each hit in the database. You can also get
an alignment mapping from query to hit, or hit to query.

=cut

use strict; 
use warnings;

use Moose;
use Carp;

use types;
use TCNPerlVars;
use pdb::pdb;

use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use File::Basename;

extends 'pdb::runBlast';

### Attributes #################################################################
################################################################################

=head1 Attributes

=over 12

=cut

=item C<db>

Name of formatted database of sequences in pdb, for blast to be run against.
File can either be present in blast path or a path can be given to the datbase,
e.g. /path/to/my/db (no file extension is neccessary)

=cut

### Attribute Builders #########################################################
################################################################################

# Note that this build method overwrites the pdb::runBastall method with a
# a database of sequences in the pdb
sub _build_db {
    return $TCNPerlVars::pdb_db
}

### Methods ####################################################################
################################################################################

=head1 Methods

=over 12

=cut

=item C<getHitStructure>

Retuns a chain object corresponding to the passed hit.

e.g.
 $hitChain = $getHomol->getHitStructure($hit);
 @hitChains = map {$getHomol->getHitStructure($_)} $getHomol->getHits();

=cut

sub getHitStructure {
    my $self = shift;
    my $hit  = shift;

    my $class = "Bio::Search::Hit::BlastHit";
    croak "getHitStructure must be passed a $class object!"
        if ref $hit ne $class;

    my ($type, $pdbCode, $chainID) = $self->_parseHitName($hit);

    croak "hit " . $hit->name() . " is not the sequence of a pdb chain!"
        if $type ne 'pdb';
    
    my $hitChain = chain->new(pdb_code => $pdbCode,
                              chain_id => $chainID);

    return $hitChain;
}

# This method parses a the pdb code and chain id from a hit.
sub _parseHitName {
    my $self = shift;
    my $hit = shift;

    my $hitName = $hit->name();
    
    # e.g. hitName = pdb|3GV2|A
    my ($type, $pdbCode, $chainID) = split(/\|/, $hitName);

    return($type, $pdbCode, $chainID);
}

=item C<getHitStructure>

This method returns a hash that maps query chain seq positions to hit query
chain seq positions. The hit object and corresponding hit chain must be passed.
The h2q arg can also be passed - if so, a map from hit to query is returned
(rather than query to hit)
e.g.
  $getHomolStructs->getAlignMap($hit, $hitChain);
  $getHomolStructs->getAlignMap($hit, $hitChain, h2q => 1)
=cut

sub getAlignMap {
    my $self     = shift;
    my $hit      = shift;
    my $hitChain = shift;

    my %arg = @_;
    
    my $queryID  = $self->query->pdbID();
    my $queryLen = scalar $self->query->get_sequence(include_missing => 1,
                                                     return_type => 1);
    my %q2h = ();
    
    while (my $hsp = $hit->next_hsp()) {
        my $aln = $hsp->get_aln();

        # Map from aln => query
        my %aln2query = map {_colFromResNum($aln, $queryID, $_) => $_}
            (1 .. $queryLen);

        # Map from aln => hit
        my %aln2hit = map {_colFromResNum($aln, $hit->name, $_) => $_}
            (1 .. $hsp->length('hit'));

        # Map from query => aln => hit
        for my $i (1 .. $aln->length()) {
            if (exists $aln2query{$i} && exists $aln2hit{$i}) {
                $q2h{$aln2query{$i}} = $aln2hit{$i};
            }
        }
    }

    if (exists $arg{h2q} && $arg{h2q}) {
        # Reverse mapping
        my %h2q = map {$q2h{$_} => $_} keys %q2h;
        return %h2q;
    }
    else {
        return %q2h;
    }
}

# This evals the wrapped column_from_residue_number call so that 0 is returned
# is residue position is not aligned, rather than an exception being thrown
sub _colFromResNum {
    my($aln, $seqName, $resNum) = @_;
    
    my $col = eval {$aln->column_from_residue_number($seqName, $resNum)};

    return $col ? $col : 0;
}

1;
__END__
 
