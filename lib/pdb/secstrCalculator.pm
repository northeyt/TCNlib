package pdb::secstrCalculator;

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
    $TCNPerlVars::pdbsecstr;
}

sub cmdStringFromInputs {
    my $self = shift;

    croak "no input supplied!" if ! $self->has_input();
    my $inputFile = pdb::pdbFunctions::getPDBFile($self->input());
    my $exec      = $self->execFilePath();

    return "$exec $inputFile";
}

sub getResID2secStructHref {
    my $self = shift;

    croak "Error while running exec" . $self->execFilePath().  ": " . $self->stderr()
        if ! $self->runExec();

    my @resIDSecStrArefs
        = map {$self->_parseResIDAndSecStrFromLine($_)} split("\n", $self->stdout());

    return map {$_->[0] => $_->[1]} @resIDSecStrArefs
}

sub _parseResIDAndSecStrFromLine {
    my $self = shift;
    my $line = shift;
    my ($resID, $resName, $secStr) = $self->_parseLine($line);

    return [$self->_formatResID($resID), $self->_formatSecStr($secStr)];
}

sub _formatSecStr {
    my $self   = shift;
    my $secStr = shift;

    # Replace "-" notation with "C" for coil
    $secStr =~ s/-/C/;
    return $secStr;
}

sub _formatResID {
    my $self  = shift;
    my $resID = shift;
    unless ($resID =~ /\./) {
        $resID = $self->_addDelimiterToResID($resID);
    }
    return $resID;
}

sub _addDelimiterToResID {
    my $self  = shift;
    my $resID = shift;
    
    croak "Delimiter is already present in resID!" if $resID =~ /\./;
    my $chainID = substr($resID, 0, 1);
    my $resSeq  = substr($resID, 1);
    return "$chainID.$resSeq";
}

sub _parseLine {
    my $self = shift;
    my $line = shift;
    return split(/\s+/, $line);
}

1;
