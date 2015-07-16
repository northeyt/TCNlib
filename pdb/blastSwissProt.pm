package pdb::blastSwissProt;

use strict; 
use warnings;
use Moose;
use sequence;

use Carp;

extends 'pdb::runBlast';

has 'pdbsws' => (
    isa  => 'pdb::pdbsws',
    is   => 'rw',
    lazy => 1,
    default => sub {pdb::pdbsws->new()},
);

sub _build_db {
    return $TCNPerlVars::swissProtDB;
}

# If reliable option has been passed, filter hits using isHitReliable.
around 'getHits' => sub {
    my $orig = shift;
    my $self = shift;

    my %arg = @_;
    
    my @hits = $self->$orig(%arg);
    
    if (exists $arg{reliable} && $arg{reliable}) {
        return grep {isHitReliable($_)} @hits;
    }
    else {
        return @hits;
    }
};

# Returns TRUE unless hit description contains any of the following words:
# putative, predicted or hypothetical
sub isHitReliable {
    my $hit = shift;

    my $desc = $hit->description;
    
    return $desc !~ /putative/i && $desc !~ /predicted/i &&
        $desc !~ /hypothetical/i;
}

sub swissProtSeqFromHit {
    my $self = shift;
    my $hit  = shift;
    
    my $ac  = $self->parseACFromHit($hit);
    my $seq = $self->pdbsws->seqFromAC($ac);
    return sequence->new(id => $ac, string => $seq);
}

sub parseACFromHit {
    my $self    = shift;
    my $hit     = shift;

    my ($ac) = $hit->name =~ /sp\|(\S{6})/;
    return $ac;
}

1;

=head1 NAME

pdb::blastSwissProt - blast pdb chains against a swissProt db

=head1 SYNOPSIS

   use pdb::blastSwissProt;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for pdb::blastSwissProt, 

Blah blah blah.

