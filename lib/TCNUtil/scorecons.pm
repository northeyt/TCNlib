package scorecons;
use Moose;
use TCNPerlVars;
use TCNUtil::write2tmp;
use Carp;
use TCNUtil::roles::consScoreCalculating;

with 'TCNUtil::roles::fileExecutor';
with 'roles::consScoreCalculating';

has 'inputAlignedSeqsString' => (
    isa => 'Str',
    is  => 'rw',
    lazy => 1,
    builder => '_buildStrFromSeqs',
);

has 'inputAlignedSeqsStringFile' => (
    isa => 'FileReadable',
    is  => 'rw',
    lazy => 1,
    builder => '_writeSeqsStr2File',
);

sub _buildStrFromSeqs {
    my $self = shift;
    
    return join("\n", map {$_->getFASTAStr()} @{$self->seqs});
}

sub _writeSeqsStr2File {
    my $self = shift;

    # Write string to temporary file
    return write2tmp->new(data => [$self->inputAlignedSeqsString])->file_name();
}
    
sub _buildExecPath {
    return $TCNPerlVars::scorecons;
}

sub cmdStringFromInputs {
    my $self = shift;

    my $inputFile = $self->inputAlignedSeqsStringFile();
    my $exec      = $self->execFilePath();
    my $flags     = $self->getFlags();

    if ($self->hasTargetSeqIndex) {
        $self->opts->{"--focus"} = $self->targetSeqIndex();
    }
    
    my $opts      = $self->getOpts();
    
    return "$exec $flags $opts $inputFile";
}

sub calcConservationScores {
    my $self = shift;

    croak "scorecons was not successful! Command run: "
        . $self->cmdStringFromInputs . " STDERR: " . $self->stderr
            if ! $self->runExec();

    return $self->parseConservationScores();
}

sub parseConservationScores {
    my $self = shift;

    my @lines = split("\n", $self->stdout);

    return map {/^([0-9.-]+)/} @lines
}

1;
