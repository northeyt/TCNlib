package pdb::blaster;
use Moose::Role;
use Carp;
use pdb::BLAST::Report;

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

has 'reportHandler' => (
    is   => 'rw',
    does => 'pdb::BLAST::ReportHandler',
    required => 1,
);

=item C<db>

Name of database for blast to be run against.

If BLAST is being run locally, then this is the location of the formatted db.
File can either be present in blast path or a path can be given to the datbase,
e.g. /path/to/my/db
Note that only the base of the db is expected, so do not include any file
extensions.

If BLAST is being run remotely, then this should be the name of the database
being searched remotely

=cut

has 'db'     => (
    is       => 'rw',
    required => 1,
);

sub _getSeqFromQuery {
    my $self = shift;
    confess "No query has been set!" if ! $self->getQuery();
    my $query = $self->getQuery();
    my $querySeq = join("", $query->get_sequence(include_missing => 1,
                                                 return_type => 1));
    return Bio::Seq->new(-id  => $query->pdbID, -seq => $querySeq);
}

requires 'runBlast';
requires 'getQuery';
requires 'setQuery';

package pdb::BLAST::Local;

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
use TCNUtil::types;
use TCNPerlVars;
use pdb;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use File::Basename;

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

=item C<query>

Query chain for blast search

=cut

has 'query' => (
    isa    => 'chain',
    is     => 'rw',
    reader => 'getQuery',
    writer => 'setQuery',
);

with 'pdb::blaster';

### Methods ####################################################################
################################################################################


=item C<runBlast>

runs blastall on query chain sequence. Assigns the resulting
Bio::SearchIO::blast object to result attribute.

getHomologueStructures->new(query => $inChain)->runBlastall();

=cut

sub runBlast {
    my $self = shift;
    
    my($blastallExecName, $blastallExecPath, $suffix)
        = fileparse($self->blastallExec());
    
    $ENV{'BLASTDIR'} = $blastallExecPath;
    
    my $blastReport = $self->_getFactory()->blastall($self->_getSeqFromQuery());

    $self->reportHandler->query($self->getQuery());
    $self->reportHandler->report($blastReport);
    return $blastReport;
}

sub _getFactory {
    my $self = shift;
    my @params = (-database => $self->db, -program => 'blastp',
                  -e => $self->evalue, @{$self->flags}, %{$self->opts});
    my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);
}

package pdb::BLAST::Remote;
use Moose;

use Bio::SeqIO;
use Bio::Tools::Run::RemoteBlast;
use Carp;

=item C<query>

Query chain for blast search

=cut

has 'query' => (
    isa    => 'chain',
    is     => 'rw',
    reader => 'getQuery',
    writer => 'setQuery',
);

with 'pdb::blaster';

sub runBlast {
    my $self = shift;

    my $factory = $self->_getFactory();
    $factory->submit_blast($self->_getSeqFromQuery());

    my ($remoteJobID) = $factory->each_rid;
    print STDERR "remote blast job id $remoteJobID: waiting...";
    while (1) {
        my $statusCodeOrReport = $factory->retrieve_blast($remoteJobID);
        if(_isStatusCode($statusCodeOrReport)) {
            
            croak "remote blast job $remoteJobID failed, error code: $statusCodeOrReport"
                if $statusCodeOrReport > 0;
            
            print STDERR ".";
            sleep 5;
        }
        else {
            print " finished\n";
            $self->reportHandler->query($self->getQuery());
            $self->reportHandler->report($statusCodeOrReport);
            return $statusCodeOrReport;
        }
    }
}

sub _isStatusCode {
    return ! ref ($_[0]);
}
    
sub _getFactory {
    my $self = shift;
    my @params = (-data => $self->db, -prog => 'blastp',
                  -expect => $self->evalue, -readmethod => 'SearchIO',
                  @{$self->flags}, %{$self->opts});
    return Bio::Tools::Run::RemoteBlast->new(@params);
}

package pdb::BLAST::Factory;
use Moose;
use Moose::Util::TypeConstraints;
use TCNPerlVars;

has 'remote' => (
    is       => 'rw',
    isa      => 'Bool',
    required => 1,
    default  => 0,
);

has 'dbType' => (
    is  => 'rw',
    isa => enum([qw(pdb swsprot)]),
);

has 'db' => (
    is      => 'rw',
    isa     => 'Str',
    lazy    => 1,
    builder => '_buildDB'
);

has 'reportHandler' => (
    is  => 'rw',
    isa => 'pdb::BLAST::ReportHandler',
    lazy => 1,
    required => 1,
    builder => '_buildReportHandler',
);

sub _buildReportHandler {
    my $self = shift;
    return $self->dbType eq 'pdb' ? pdb::BLAST::Report::PDBseq->new()
        : pdb::BLAST::Report::SwissProt->new();
}

sub _buildDB {
    my $self = shift;
    if ($self->remote) {
        return $self->dbType eq 'pdb' ? 'pdb'
            : 'swissprot';
    }
    else {
        return $self->dbType eq 'pdb' ? $TCNPerlVars::pdb_db
            : $TCNPerlVars::swissProtDB;
    }
}

sub getBlaster {
    my $self = shift;
    my @args = @_;

    if ($self->remote()) {
        return pdb::BLAST::Remote->new(db => $self->db,
                                       reportHandler => $self->reportHandler,
                                       @args);
    }
    else {
        return pdb::BLAST::Local->new(db => $self->db,
                                      reportHandler => $self->reportHandler,
                                      @args);
    }
}

1;
__END__
 
