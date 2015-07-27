package pdb::BLAST::ReportHandler;
use Moose::Role;

has 'query' => (
    is => 'rw',
);

has 'report'  => (
    is        => 'rw',
    isa       => 'Bio::SearchIO::blast',
    predicate => 'has_report',
);

has 'seqID' => (
    is => 'rw',
    isa => 'Num',
    default => 0,
);

=item C<getHits>

Returns an array of Bio::Search::Hit::BlastHit objects. If a
Bio::SearchIO::blast object is passed using report arg then this will be used -
otherwise, runBlastall is called to obtain a Bio::SearchIO::blast object.

A sequence identity threshold can also be passed with seq_id

# The @hits will be the same as @hits2 ...
$getHStruct = getHomologueStructures->new(query => $inChain);
@hits  = $getHStruct->getHits(seq_id => 0.95);
@hits2 = $getHStruct->getHits(report => $self->report(), seq_id => 0.95);


=cut

sub getHits {
    my $self   = shift;
    my $blastReport = $self->report();
        
    my @hits = ();

    # Used to check hit names against
    my $queryID
        = join("|", ("pdb", uc($self->query->pdb_code),
                     $self->query->chain_id));
    
    while (my $result = $blastReport->next_result ) {
        while (my $hit = $result->next_hit ) {
            # Avoid inclusion of query in hits array
            next if $hit->name eq $queryID;

            # Skip if hit has seq_id less than threshold
            if ($self->seqID) {
                next if $hit->frac_identical() < $self->seqID;
            }
            push(@hits, $hit);
        }
    }
    return @hits;
}

package pdb::BLAST::Report::PDBseq;
use Moose;
use Carp;

with 'pdb::BLAST::ReportHandler';

sub getHitStructures {
    my $self = shift;
    return map {$self->getHitStructure($_)} $self->getHits();
}

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

=item C<getAlignMap>

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
    my $queryLen = scalar $self->query()->get_sequence(include_missing => 1,
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

package pdb::BLAST::Report::SwissProt;
use Moose;
use pdb::pdbsws;
use sequence;

with 'pdb::BLAST::ReportHandler';

has 'pdbsws' => (
    isa  => 'pdb::pdbsws',
    is   => 'rw',
    lazy => 1,
    default => sub {pdb::pdbsws->new()},
);

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
