package FOSTA::FEPFinder;
use Moose::Role;
use TCNUtil::sequence;

requires 'getReliableFEPSequencesFromSwissProtID';

package FOSTA::Remote;
use Moose;
use Carp;
use TCNUtil::WEB::FOSTA;
use UNIPROT;
use TCNUtil::sequence;

with 'FOSTA::FEPFinder';

sub getReliableFEPSequencesFromSwissProtID {
    my $self        = shift;
    my $swissProtID = shift;
    my ($errorCode, @fepIDs) = WEB::FOSTA::GetFOSTA($swissProtID);
    croak "Error code $errorCode returned by WEB::FOSTA::GetFOSTA!"
        if $errorCode;
    my @fepFASTAs = map {UNIPROT::GetFASTA($_, -remote => 1)} @fepIDs;
    return map {sequence->new($_)} @fepFASTAs;
}

package FOSTA::Local;
use Moose;
use Carp;
use DBI;

with 'FOSTA::FEPFinder';

has 'FOSTADBH' => (
    is => 'rw',
    required => 1,
    lazy => 1,
    builder => '_buildFOSTADBH',
);

sub _buildFOSTADBH {
    my $dbname3 = "fosta";
    my $dbserver3 = 'acrm8.biochem.ucl.ac.uk';
    my $datasource3 = "dbi:Pg:dbname=$dbname3;host=$dbserver3";
    
    my $fostadbh = DBI->connect ($datasource3)
        || die "Cannot connect to $dbname3 database.\n";

    return $fostadbh;
}

sub getReliableFEPSequencesFromSwissProtID {
    my $self     = shift;
    my $sprot_id = shift;
    
    #query FOSTA for all FEPs of that sprot ID
    my @FEPIDs = eval {$self->getFEPIDs($sprot_id)};
    croak "No FEP IDs returned for SwissProt ID $sprot_id: $@"
        if ! @FEPIDs;
    
    return $self->getSequences($sprot_id, @FEPIDs);
}
            
# From a SwissProt id, finds FOSTA family id and family reliability.
# Reliabilility = 1 if reliable, 0 = unreliable.
sub getFOSTAFamIDAndReliability {
    my $self = shift;
    my $id   = shift;
    
    my $sql3 = "SELECT fd.unreliable, f.fosta_family
                FROM fosta_descriptions fd, feps f
                WHERE fd.id='$id'
                AND f.id='$id'";
    my ($unreliable, $family) = $self->FOSTADBH->selectrow_array($sql3);

    # Return reliability rather than unrelability
    my $reliable = $unreliable ? 0 : 1;
    
    return $reliable, $family;
}

#checks whether entry is unreliable for orthofind and finds functionally equivalent proteins (FEPs) for reliable query proteins
#returns array, first element is query id, others are its FEPs
sub getFEPIDs {
    my $self = shift;
    my $id   = shift;
    
    my @familyFEPIDs = ();  

    chomp $id;

    my ($reliable, $family) = $self->getFOSTAFamIDAndReliability($id);

    croak "No family id returned for id $id" if ! $family;
    croak "ID $id is not a reliable query"   if ! $reliable; 
    
    @familyFEPIDs = $self->getFEPIDsFromFamID($id, $family);
    
    if (! @familyFEPIDs){
        croak "No FEP IDs returned with family id $family, query id $id";
    }    
    
    return @familyFEPIDs;
}

sub getFEPIDsFromFamID {
    my $self   = shift;
    my $id     = shift;
    my $family = shift;
    
    #ask Lisa is this query OK
    my $sql4 = "SELECT id FROM feps
                WHERE fosta_family = '$family'
                AND id != '$id'
                AND NOT unreliable;";

    my @FEPIDs = ();
    
    my $fostasth = $self->FOSTADBH->prepare($sql4);
    if($fostasth->execute){
        while(my ($fep_id) = $fostasth->fetchrow_array){
            push (@FEPIDs, $fep_id);                        
        }
    }
    return @FEPIDs;
}    

sub getSequenceFromID {
    my $self   = shift;
    my $seq_id = shift;
    my $seqStr = $self->getSequenceStrFromID($seq_id);
    croak "No sequence string found for id $seq_id" if ! $seqStr;
    return sequence->new(string => $seqStr, id => $seq_id)
}

sub getSequenceStrFromID {
    my $self   = shift;
    my $seq_id = shift;
    
    my $seq_sql = "SELECT sequence
                   FROM fosta_sequences
                   WHERE id='$seq_id'";
    my $seq = $self->FOSTADBH->selectrow_array($seq_sql);
    return $seq;
}

sub getSequences {
    my $self = shift;
    
    my @seq_ids = @_;
    my @seqs    = ();
    foreach my $seq_id (@seq_ids){
        my $seq = eval {$self->getSequenceFromID($seq_id)};
        next if ! $seq;
        push(@seqs, $seq);
    }
    return @seqs;
}

package FOSTA::Factory;
use Moose;

has 'remote' => (
    is       => 'rw',
    isa      => 'Bool',
    required => 1,
    default  => 0,
);

sub getFOSTA {
    my $self = shift;
    my @args = @_;

    if ($self->remote) {
        return FOSTA::Remote->new(@args);
    }
    else {
        return FOSTA::Local->new(@args);
    }
}

1;
