#!/usr/bin/env perl
use strict;
use warnings;
use TCNUtil::WEKA;
use Carp;
use Getopt::Long;

my $undefLabel      = "U";
my $defaultOutForm  = 0;
GetOptions("u=s", \$undefLabel,
           "f",   \$defaultOutForm);

@ARGV == 3 or Usage();

my ($negLabel, $posLabel, $csv) = @ARGV;

open(my $IN, "<", $csv) or die "Cannot open file $csv, $!";

my $WEKA = WEKA->new(posLabel => $posLabel, negLabel => $negLabel,
                     undefLabel => $undefLabel);

my $ouputFormat = $defaultOutForm ? 'DEF' : 'CSV';
my $table       = $WEKA->parseTableFromOutput([<$IN>], $ouputFormat);
$table->print_all();

sub Usage {
    print <<EOF;
$0 usage : [-u undefLabel] -f negativeClass positiveClass CSV
 -u : string that indicates a undefined class label. DEFAULT = U
 -f : set this flag if WEKA output file is in the WEKA default format,
      rather than CSV.

 negativeClass: negative class label (e.g. "no disease")
 positiveClass: positive class label (e.g. "disease")

This script calculates performance statistics from a .csv file output from WEKA
EOF
    exit(1)
}
