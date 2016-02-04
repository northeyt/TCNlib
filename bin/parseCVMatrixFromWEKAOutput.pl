#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my $file;
GetOptions("f=s", \$file);
my $string = "";

if ($file) {
    # Read from file
    open(my $IN, "<", $file) or die "Cannot open file $file, $!\n";
    $string = join("", <$IN>);
}
else {
    # Read from STDIN
    $string = join("", <>);
}

my ($subStr) = $string =~ /(?:=== Stratified cross-validation)* .* === Confusion Matrix === (.*)/g;
my ($tp, $fn, $fp, $tn) = $subStr =~ /[^=](?= (\d+))/g;
print `MCC.pl $tp $tn $fp $fn`;
