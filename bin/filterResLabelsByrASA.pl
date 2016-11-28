#!/usr/bin/env perl
use strict;
use warnings;
use Carp;
use Getopt::Long;
my $relASAFile;
my $relASAMin;
my $outFile;

GetOptions("r=s", \$relASAFile,
           "o=s", \$outFile,
           "m=s", \$relASAMin);

$relASAMin = 10 if ! defined $relASAMin;
$outFile = "relASAFiltered-min$relASAMin.residue.labels" if ! defined $outFile;
open(my $OUT, ">", $outFile) or die "Cannot open file $outFile, $!";

@ARGV or Usage();
my $resLabelsFile = shift @ARGV;

if ($relASAFile) {
    open(my $IN, "<", $relASAFile) or die "Cannot open file $relASAFile, $!";
    my %resID2RelASA;
    while (my $line = <$IN>) {
        chomp $line;
        my ($resID, $relASA) =  $line =~ /(.*):(.*)/;
        $resID = uc(substr($resID, 0, 4)) . substr($resID, 4);
        $resID2RelASA{$resID} = $relASA;
    }
    my @residueLabels = getResidueLabelsSortedByResID($resLabelsFile);
    foreach my $resLabel (@residueLabels) {
        my ($resID, $resLabel) =  $resLabel =~ /(.*):(.*)/;
        if (! exists $resID2RelASA{$resID}) {
            print {*STDERR} "WARNING: $resID has no relASA value!\n";
            next;
        }
        my $relASA = $resID2RelASA{$resID};
        print {$OUT} join(":", $resID, $resLabel) . "\n" if $relASA > $relASAMin;
    }
}
else {
    print "Getting residue relASAs: this can take a while if your data set is large ...\n";
    my @residueRelASAs = getResidueRelASAsSortedByResID($resLabelsFile);
    my @residueLabels  = getResidueLabelsSortedByResID($resLabelsFile);
    
    croak "Number of residue relASA values does not match number of residue labels"
        if @residueRelASAs != @residueLabels;
    
    print "Filtering residue labels by relASA\n";
    for (my $i = 0 ; $i < @residueRelASAs ; ++$i) {
        my $resIDAndRelASA = $residueRelASAs[$i];
        my $resIDAndLabel  = $residueLabels[$i];
        
        my ($firstResID, $relASA) = $resIDAndRelASA =~ /(.*):(.*)/;
        my ($secndResID, $resLab) = $resIDAndLabel  =~ /(.*):(.*)/;
        
        croak "Residue labels $firstResID and $secndResID do not match!"
            if $firstResID ne $secndResID;
        
        print {$OUT} "$resIDAndLabel\n" if $relASA > $relASAMin;
    }
}

close $OUT;
print "Finished\n";

sub getResidueRelASAsSortedByResID {
    my $resLabelsFile = shift;
    my $resRelASAs = `cut -d: -f1,2 $resLabelsFile | uniq | tr -d : | getResidueRelASAs.pl | sortResidues.sh`;
    return split(/\n/, $resRelASAs);
}

sub getResidueLabelsSortedByResID {
    my $resLabelsFile = shift;
    my $sortedResLabels = `cat $resLabelsFile | sortResidues.sh`;
    return split(/\n/, $sortedResLabels);
}

sub Usage {
    print <<EOF;
$0 -r relASA-file -o out-file residue-labels-file

Opts:
    -r file containing residue relASA values. If not supplied, relASA values
       will be calculated
    -o output file. If not supplied, default of
       relASAFiltered-min{MIN}.residue.labels will be used, where {MIN} is
       relASA minimum.
    -m relASA minimum used to filter residues. Default = 10
EOF
    exit(1);
}
