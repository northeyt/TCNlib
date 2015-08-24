package scorecons;
use Moose;
use TCNPerlVars;
use TCNUtil::write2tmp;
use Carp;

with 'TCNUtil::roles::fileExecutor';
with 'TCNUtil::roles::consScoreCalculating';

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

    print "DEBUG: seq string\n" . $self->inputAlignedSeqsString();
    
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
    my $self            = shift;
    my @lines           = split("\n", $self->stdout);
    my @results         = map {PositionResult->new($_)}  @lines;
    my @filteredResults = grep {$self->_resultPassesFilter($_)} @results; 
    return map {$_->score} @filteredResults; 
}

sub _resultPassesFilter {
    my $self   = shift;
    my $result = shift;
    return ! $self->hasTargetSeqIndex ? 1
        :  $result->doesMapToTargetWithIndex($self->hasTargetSeqIndex)
}

package PositionResult;
use Moose;
use TCNUtil::GLOBAL qw(rm_trail);

has 'score' => (
    isa => 'Num',
    is => 'rw',
    required => 1
);

has '_sequencePositionStr' => (
    isa => 'Str',
    is => 'rw',
    required => 1
);

sub doesMapToTargetWithIndex {
    my $self  = shift;
    my $index = shift;
    return substr($self->_sequencePositionStr(), $index, 1) eq '-' ? 0 : 1;
}

# Allows lazy object creation: PositionResult->new($myString)
around BUILDARGS => sub {
    my $orig  = shift;
    my $class = shift;
    
    if (@_ == 1) {
        my ($score, $seqPosStr) = _parseString($_[0]);
        return $class->$orig(score => $score,
                             _sequencePositionStr => $seqPosStr);
    }
    else {
        return $class->$orig(@_);
    }
};

sub _parseString {
    my $string = rm_trail(shift);
    my ($position, $score, $seqPositionString) = split(/\s+/, $string);
    return ($score, $seqPositionString);
}

1;
