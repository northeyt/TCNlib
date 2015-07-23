package pdb::hbondFinder;
use Moose;
use Moose::Util::TypeConstraints;
use TCNPerlVars;
use pdb::pdbFunctions;
use Carp;

with 'roles::fileExecutor';

has 'input' => (
    is => 'rw',
    predicate => 'has_input',
);

# Currently can only handle protein-protein hbonds
has 'bondType' => (
    is => 'rw',
    isa => enum([qw(pp)]),
    default => 'pp'
);

sub _buildExecPath {
    $TCNPerlVars::pdbhbond;
}

sub cmdStringFromInputs {
    my $self = shift;

    croak "no input supplied!" if ! $self->has_input();
    my $inputFile = pdb::pdbFunctions::getPDBFile($self->input());
    my $exec      = $self->execFilePath();

    return "$exec $inputFile";
}

sub getHbonds {
    my $self = shift;

    croak "Error while running exec" . $self->execFilePath().  ": " . $self->stderr()
        if ! $self->runExec();

    return $self->_parseHydrogenBondsFromOutput();
}

sub _parseHydrogenBondsFromOutput {
    my $self = shift;
    my @lines = split("\n", $self->stdout);

    my @hbs = ();
    
    my $inType = 0;
    foreach my $line (@lines) {
        if ($line =~ /^TYPE: (\w+)hbonds/) {
            $inType = $1 eq $self->bondType() ? 1 : 0;
            next;
        }

        if ($inType && $line !~ /^#/) {
            push(@hbs, $self->_parseHbFromLine($line));
        }
    }
    return @hbs;
}

sub _parseHbFromLine {
    my $self = shift;
    my $line = shift;

    my ($donorSerial, $acceptorSerial) = $line =~ /(\d+) \s+ (\d+)/xms;
    return pdb::Hb->new(donorSerial => $donorSerial,
                        acceptorSerial => $acceptorSerial);
}

package pdb::Hb;
use Moose;

has 'donorSerial'   => (isa => 'Int', is => 'rw', required => 1);
has 'acceptorSerial' => (isa => 'Int', is => 'rw', required => 1);

1;
