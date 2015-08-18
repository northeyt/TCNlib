package TCNUtil::roles::consScoreCalculating;
use Moose::Role;
use Carp;

requires 'calcConservationScores';

has 'seqs' => (
    is      => 'rw',
    isa     => 'ArrayRef[canGetFASTAStr]',
    lazy    => 1,
    builder => '_buildSeqsFromOtherInput',
);

has 'targetSeqIndex' => (
    is  => 'rw',
    isa => 'Int',
    predicate => 'hasTargetSeqIndex',
);

sub _buildSeqsFromOtherInput {
    croak "_buildSeqsFromOtherInput has not been implemented for this class!";
}

1;
