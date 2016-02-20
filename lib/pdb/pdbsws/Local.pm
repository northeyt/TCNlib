package pdb::pdbsws::Local;

use Moose;
use TCNUtil::types;
use DBI;
use TCNPerlVars;
use Carp;

with 'pdb::pdbsws::querier';

has 'dbhost' => (
    isa => 'Str',
    is => 'rw',
    default => $TCNPerlVars::pghost,
);

has 'dbname' => (
    isa => 'Str',
    is => 'rw',
    default => 'pdbsws',
);

has 'dbh' => (
    isa => 'Ref',
    is => 'ro',
    lazy => 1,
    builder => '_get_dbh', 
);

sub _get_dbh {
    my $self = shift;

    my $dbname = $self->dbname();
    my $dbhost = $self->dbhost();
    
    my $dbh = DBI->connect("dbi:Pg:dbname=$dbname;host=$dbhost");
    croak "Could not connect to pdbsws" if ! $dbh;

    return $dbh;
}

sub getIDsFromPDBCodeAndChainID {
    my $self     = shift;
    my $pdbCode  = shift;
    my $chainID  = shift;
    my @acs      = $self->getACsFromPDBCodeAndChainID($pdbCode, $chainID);
    return map { $self->_getSwissProtIDFromAC($_) } @acs;
}

sub getACsFromPDBCodeAndChainID {
    my $self     = shift;
    my $pdbCode  = shift;
    my $chain    = shift;
    
    my $sql = "SELECT ac
               FROM pdbsws
               WHERE pdb = '$pdbCode'
               AND chain = '$chain'
               AND valid = 't'
               AND aligned = 't'
               AND ac != 'SHORT'
               AND ac != 'DNA'
               AND ac != 'ERROR';";

    my $sth = $self->dbh->prepare($sql);

    my @ac = ();
    
    if ($sth->execute) {
        while ( my ($pdb_ac) = $sth->fetchrow_array ){
            push(@ac, $pdb_ac);
        }
    }
    return @ac;
}

sub get_ac {
    my $self = shift;
    croak "get_ac is no longer supported! please replace with getSwissProtACFromPDBCodeAndChainID"
        . " and note that pdb code and chain id are now passed separately"
            . " \ne.g. getSwissProtACFromPDBCodeAndChainID('4hou', 'A')";
}

sub _getSwissProtIDFromAC {
    my $self     = shift;
    my $sprot_ac = shift;
    
    my $sql = "SELECT i.id FROM idac i, acac a WHERE a.altac = '$sprot_ac' AND i.ac = a.ac;";
    my $sprot_id = $self->dbh->selectrow_array($sql);
    return $sprot_id;
}

sub mapResSeq2SwissProtNum {
    my $self     = shift; 
    my $pdbCode  = shift;
    my $chainID  = shift;
    my $targetAC = shift;
    
    my $sql = "SELECT resid, pdbaa, ac, swsaa, swscount
               FROM alignment
               WHERE pdb = '$pdbCode'
               AND chain = '$chainID';";
    
    my $pdbswssth = $self->dbh->prepare($sql);

    my %resSeq2SprotResNum = ();
    if($pdbswssth->execute){
        while (my ($pdbResSeq, $pdbRes, $sprotAC, $sprotRes, $sprotResNum)
                   = $pdbswssth->fetchrow_array){
            if ($sprotAC eq $targetAC) {
                $resSeq2SprotResNum{$pdbResSeq} = $sprotResNum;
            }
        }
    }
    return %resSeq2SprotResNum;
}

1;
