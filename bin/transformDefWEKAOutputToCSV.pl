#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib";
use TCNUtil::WEKA;

die "Please supply a weka output file in default format."
    if ! @ARGV;

my @additionalFields = ();
while (@ARGV > 1) {
    push(@additionalFields, shift @ARGV);
}

my $inFile = shift @ARGV;
open(my $IN, "<", $inFile) or die "Cannot open in file $inFile, $!";

my $baseNumOfFields = 4;
my $weka            = WEKA->new();
my @CSVLines        = ();
my $finalNumOfAdditionalFields = 0;

while (my $line = <$IN>) {
    chomp $line;
    next if ! $line;
    
    my @fields   = $weka->parseDefaultOutputLine($line);
    my $numOfAdditionalFields = scalar @fields - $baseNumOfFields;
    $finalNumOfAdditionalFields = $numOfAdditionalFields
        if $numOfAdditionalFields > $finalNumOfAdditionalFields;
    
    my ($instNum, $valueLabel, $predLabel, $predValue, @remainingFields)
        = @fields;
    
    my $error   = $valueLabel eq $predLabel ? '' : '+';
    my $csvLine = join(",", ($instNum, $valueLabel, $predLabel, $error,
                             $predValue, @remainingFields));
    push(@CSVLines, "$csvLine\n");
}

my @fields = ('inst#', 'actual', 'predicted', 'error', 'prediction');

if (@additionalFields) {
    push(@fields, @additionalFields);
}
else {
    push (@additionalFields, '?') for (1 .. $finalNumOfAdditionalFields);
}   

print join(",", @fields) . "\n";
print @CSVLines;

