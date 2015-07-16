############################### CLASS ##########################################
############################ pdb::runBlast #####################################
################################################################################

package pdb::runBlast;

=head1 NAME

pdb::runBlast - do a BLAST search with a query chain against a blast db.

=cut

=head1 SYNOPSIS

TODO

=cut


=head1 DESCRIPTION

pdb::getHomologueStructures allows you to use a chain object as a query for a
BLAST search against a database of sequences. This is the base class of
pdb::getHomologueStructres.

=cut


use Moose;
use Carp;

use types;
use TCNPerlVars;
use pdb::pdb;

use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use File::Basename;


=item C<query>

Query chain for blast search

=cut

has 'query' => (
    isa => 'chain',
    is  => 'rw',
);


=item <evalue>

Maximum evalue threshold

=cut

has 'evalue' => (
    is => 'Num',
    is => 'rw',
    default => 10,
);

has 'flags' => (
    is => 'rw',
    isa => 'ArrayRef',
    default => sub { [] },
    lazy => 1,
);

has 'opts' => (
    is => 'rw',
    isa => 'HashRef',
    default => sub { {} },
    lazy => 1,
);

=item C<blastallExec>

blastall executable to be run

=cut

has 'blastallExec' => (
    isa => 'FileExecutable',
    is => 'rw',
    required => 1,
    lazy => 1,
    default => $TCNPerlVars::blastall,
);


=item C<db>

Name of formatted databasefor blast to be run against.
File can either be present in blast path or a path can be given to the datbase,
e.g. /path/to/my/db (no file extension is neccessary)

=cut

# Note that no FileReadable type is declared, as blastall only expects a base
# path that it then appends file extensions to actually find the database
# component files.
has 'db' => (
    is => 'rw',
    required => 1,
    lazy => 1,
    builder => '_build_db',
);

=item C<report>

Bio::SearchIO::blast object that is report of lastest blastall run.

=cut

has 'report' => (
    isa => 'Bio::SearchIO::blast',
    is => 'rw',
    predicate => 'has_report',
);

sub _build_db {
    croak "Database must be specified!";
}

### Methods ####################################################################
################################################################################


=item C<runBlastall>

runs blastall on query chain sequence. Assigns the resulting
Bio::SearchIO::blast object to result attribute.

getHomologueStructures->new(query => $inChain)->runBlastall();

=cut

sub runBlastall {
    my $self = shift;
    
    my($blastallExecName, $blastallExecPath, $suffix)
        = fileparse($self->blastallExec());
    
    $ENV{'BLASTDIR'} = $blastallExecPath;
    
    my @params = (-database => $self->db, -program => 'blastp',
                  -e => $self->evalue, @{$self->flags}, %{$self->opts});
    my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);
    
    my $query = $self->query();
    my $querySeq = join("", $query->get_sequence(include_missing => 1,
                                                 return_type => 1));
    
    my $seq = Bio::Seq->new(-id  => $query->pdbID,
                            -seq => $querySeq);
    
    my $blastReport = $factory->blastall($seq);

    $self->report($blastReport);
    
    return $blastReport;
}

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
    my %arg = @_;
        
    my $blastReport;
    if (! exists $arg{report}) {
        $blastReport = $self->runBlastall();
    }
    else {
        $blastReport = $arg{report};
    }

    my $seq_id = $arg{seq_id} ? $arg{seq_id} : 0;
    
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
            if ($seq_id) {
                next if $hit->frac_identical() < $seq_id;
            }
            push(@hits, $hit);
        }
    }
    return @hits;
}

############################## END OF CLASS ####################################
################################################################################


1;
__END__
 
