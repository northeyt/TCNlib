#!/usr/bin/env perl
use strict;
use warnings;
use TCNUtil::WEKA;
use Carp;
use Getopt::Long;

my $undefLabel = "U";
GetOptions("u=s", \$undefLabel);

@ARGV == 3 or Usage();

my ($negLabel, $posLabel, $csv) = @ARGV;

open(my $IN, "<", $csv) or die "Cannot open file $csv, $!";

my $WEKA = WEKA->new(posLabel => $posLabel, negLabel => $negLabel,
                     undefLabel => $undefLabel);

my $table = $WEKA->parseTableFromOutput([<$IN>]);
$table->print_all();

sub Usage {
    print <<EOF;
$0 usage : [-u undefLabel] negativeClass positiveClass CSV
 undefLabel   : string that indicates a undefined class label.
                DEFAULT = U
 negativeClass: negative class label (e.g. "no disease")
 positiveClass: positive class label (e.g. "disease")

This script calculates performance statistics from a .csv file output from WEKA
EOF
    exit(1)
}
