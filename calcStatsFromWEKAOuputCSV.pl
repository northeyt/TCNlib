#!/acrm/usr/local/bin/perl
use strict;
use warnings;
use WEKA;
use Carp;
use Getopt::Long;

@ARGV == 3 or Usage();

my ($negLabel, $posLabel, $csv) = @ARGV;

open(my $IN, "<", $csv) or die "Cannot open file $csv, $!";

my $WEKA = WEKA->new(posLabel => $posLabel, negLabel => $negLabel);

my $table = $WEKA->parseTableFromOutput([<$IN>]);
$table->print_all();

sub Usage {
    print <<EOF;
$0 usage : negativeClass positiveClass CSV

 negativeClass: negative class label (e.g. "no disease")
 positiveClass: positive class label (e.g. "disease")

This script calculates performance statistics from a .csv file output from WEKA
EOF
    exit(1)
}
