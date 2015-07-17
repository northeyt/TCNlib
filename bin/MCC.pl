#!/usr/bin/env perl

@ARGV == 4 or die "Supply four numbers: TP, TN, FP, FN\n";
my ($TP, $TN, $FP, $FN) = @ARGV;

my $MCC
    = ( ($TP * $TN) - ($FP * $FN) )
    / sqrt ( ($TP+$FP)*($TP+$FN)*($TN+$FP)*($TN+$FN) );

print $MCC . "\n";

