#!/acrm/usr/local/bin/perl -s

use strict;
use warnings;

my $string = "";

if ($::f) {
    # Read from file
    die "File name must be supplied with -f option!\n" if ! @ARGV;
    open(my $IN, "<", $ARGV[0]) or die "Cannot open file $ARGV[0], $!\n";

    $string = join("", <$IN>);
}
else {
    # Read from STDIN
    $string = join("", <>);
}

my ($subStr) = $string =~ /(?:=== Stratified cross-validation)* .* === Confusion Matrix === (.*)/g;

my ($tp, $fn, $fp, $tn) = $subStr =~ /[^=](?= (\d+))/g;

print `~/scripts/lib/MCC.pl $tp $tn $fp $fn`;
