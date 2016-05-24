#!/usr/bin/env perl
use strict;
use warnings;
use TCNUtil::WEKA;
use Carp;
use Getopt::Long;

my $undefLabel      = "U";
my $defaultOutForm  = 0;
my $csvOutput       = 0;
my $noHeader = 0;

GetOptions("u=s", \$undefLabel,
           "f",   \$defaultOutForm,
           "c",   \$csvOutput,
           "n",   \$noHeader);

@ARGV or Usage();

my ($negLabel, $posLabel, @inCSVs) = @ARGV;

my $outCSVHeader;
my @fields;

foreach my $inputCSV (@inCSVs) {
    open(my $IN, "<", $inputCSV) or die "Cannot open file $inputCSV, $!";
    
    my $WEKA = WEKA->new(posLabel => $posLabel, negLabel => $negLabel,
                         undefLabel => $undefLabel);

    
    my $ouputFormat = $defaultOutForm ? 'DEF' : 'CSV';
    my $table       = $WEKA->parseTableFromOutput([<$IN>], $ouputFormat, $noHeader);
    if ($csvOutput) {
        if (! $outCSVHeader) {
            @fields = ("input_file", $table->metrics_array);
            $outCSVHeader = join(",", @fields);
            print $outCSVHeader . "\n";
        }
        my %valueForField = $table->hash_all(printable => 1);
        my @metrics = map {$valueForField{$_}} @fields[1..$#fields];
        print join(",", ($inputCSV, @metrics)) . "\n";
    }
    else {
        $table->print_all();
    }
}

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
