package TestsFor::scorecons;
use Test::Class::Moose;

use TCNUtil::scorecons;

has 'testObj' => (
    isa => 'scorecons',
    is  => 'rw',
);
    
sub test_startup {
    my $test     = shift;
    my $testFile = "test.aln.fasta";
    my $testObj  = scorecons->new(inputAlignedSeqsStringFile =>  $testFile);
    print $test->test_report->current_class()->name() . "\n";
    $test->testObj($testObj);
}

sub test_calcConservationScores {
    my $test = shift;

    my @gotScores = $test->testObj->calcConservationScores();
    is(scalar @gotScores, 209,
       "A score is parsed for each MSA position");
    is($gotScores[0], 0.149,
       "and first score is correct");
}

sub test_targetSeqIndex {
    my $test = shift;
    
    $test->testObj->targetSeqIndex(2);
    my @gotScores = $test->testObj->calcConservationScores();

    is(scalar @gotScores, 189,
       "A score for each position of the target sequence only is returned");
    is($gotScores[0], 0.247,
       "and first score is correct");
}

1;
