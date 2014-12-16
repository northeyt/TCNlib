#!/acrm/usr/local/bin/perl
use strict;
use warnings;
use confusion_table;
use Carp;
use Getopt::Long;

my $tp = 0;
my $fp = 0;
my $tn = 0;
my $fn = 0;

GetOptions("t=i", \$tp,
           "p=i", \$fp,
           "f=i", \$fn,
           "n=i", \$tn);

my %type2Freq
    = (TruePos => $tp, FalsePos => $fp, TrueNeg => $tn, FalseNeg => $fn); 

my %type2ValueAndPred = (TruePos  => {value => 1, prediction => 1},
                         FalsePos => {value => 0, prediction => 1},
                         FalseNeg => {value => 1, prediction => 0},
                         TrueNeg  => {value => 0, prediction => 0});

@ARGV == 3 or Usage();

my $table = confusion_table->new(item_class => 'instance');

my ($negLabel, $posLabel, $csv) = @ARGV;

open(my $IN, "<", $csv) or die "Cannot open file $csv, $!";

my $reachedHeader = 0;
my %labelMap = ($negLabel => 0, $posLabel => 1);

while (my $line = <$IN>) {
    if (! $reachedHeader) {
        $reachedHeader = 1 if $line =~ /^inst#/;
        next;
    }

    next if $line =~ /^\n$/;
    
    my($instNum, $valueLabel, $predLabel, $err, $predValue)
        = split(",", $line);
    
    ($valueLabel, $predLabel) = map {[split(":", $_)]->[1]}
        ($valueLabel, $predLabel);

    foreach my $label ($valueLabel, $predLabel) {
        croak "labels do not match user-input labels!"
            if ! exists $labelMap{$label};
    }

    my $obj = bless {}, 'instance';
    
    my $instance = datum->new(object => $obj,
                              value => $labelMap{$valueLabel},
                              prediction => $labelMap{$predLabel});

    $table->add_datum($instance);
   
}

foreach my $type (keys %type2ValueAndPred) {
    my $range = $type2Freq{$type} - 1; # This ensures that no values are added if freq = 0
    
    my %valAndPred = %{$type2ValueAndPred{$type}};

    for (0 .. $range) {
        my $obj = bless {}, 'instance';
        my $instance = datum->new(object => $obj, %valAndPred);
        
        $table->add_datum($instance);
    }
}

$table->print_all();

sub Usage {
    print <<EOF;
$0 usage : [-t INT] [-p INT] [-f INT] [-n INT] negativeClass positiveClass CSV

 -t, -p, -f, -n : specify a number of instances to add to confusion table,
                  where -t = true positive,  -p = false positive,
                        -f = false negative, -n = true negative

 negativeClass: negative class label (e.g. "no disease")
 positiveClass: positive class label (e.g. "disease")

This script calculates performance statistics from a .csv file output from WEKA
EOF
    exit(1)
}
