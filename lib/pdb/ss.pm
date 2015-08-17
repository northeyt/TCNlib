package pdb::ssFinder;
use Moose;
use TCNUtil::types;
use TCNPerlVars;
use pdb::pdbFunctions;
use Carp;

with 'TCNUtil::roles::fileExecutor';

has 'input' => (
    is => 'rw',
    predicate => 'has_input',
);

sub _buildExecPath {
    $TCNPerlVars::pdbsslist
}

sub cmdStringFromInputs {
    my $self = shift;

    croak "no input supplied!" if ! $self->has_input();
    my $inputFile = pdb::pdbFunctions::getPDBFile($self->input());
    my $exec      = $self->execFilePath();

    return "$exec $inputFile";
}

sub getssArray {
    my $self = shift;

    croak "Error while running exec" . $self->execFilePath().  ": " . $self->stderr()
        if ! $self->runExec();

    return map {$self->_getssFromLine($_)} split("\n", $self->stdout());
}

sub _getssFromLine {
    my $self = shift;
    my $line = shift;

    # example line:
    #    A23  Atom   164 :   A88  Atom   682 : 2.076
    
    my($resIDA, $atomSerialA,
       $resIDB, $atomSerialB,
       $distance) = $line =~ /(\S+) \s* Atom \s* ([0-9]+) \s* : \s*
                              (\S+) \s* Atom \s* ([0-9]+) \s* : \s*
                              ([0-9.]+)/xms;
    
    return pdb::ss->new(atomSerialPairAref => [$atomSerialA, $atomSerialB],
                        resIDPairAref      => [$resIDA,      $resIDB],
                        distance => $distance);
}

package pdb::ss;
use Moose;

has 'atomSerialPairAref' => (
    is => 'rw',
    isa => 'ArrayRef[Int]',
    required => 1,
);

has 'resIDPairAref' => (
    is => 'rw',
    isa => 'ArrayRef[Str]',
    required => 1,
);

has 'distance' => (
    is => 'rw',
    isa => 'Num',
    required => 1,
);

1;
