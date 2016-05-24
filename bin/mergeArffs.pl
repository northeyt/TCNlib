#!/usr/bin/env perl
# mergeArffs.pl --- merge multiple arff files into a single arff file
# Author: Tom Northey <zcbtfo4@acrm18>
# Created: 27 Mar 2015
# Version: 0.01

use warnings;
use strict;
use TCNUtil::ARFF;
use Getopt::Long;

my $outputFile = "";
GetOptions("o=s", \$outputFile);

my $OUT = $outputFile ? getFH($outputFile) : *STDOUT;

my @arffs = @ARGV;

my $firstArff = "";

foreach my $arffFile (@arffs) {
    my $arff = ARFF->new(file => $arffFile);
    if (! $firstArff) {
        $firstArff = $arff;
        next;
    }

    $firstArff->mergeArff($arff);
}

print {$OUT} $firstArff->arff2String;

sub getFH {
    my $fName = shift;
    open(my $OUT, ">", $fName) or die "Cannot open file $fName, $!";
    return $OUT;
}
