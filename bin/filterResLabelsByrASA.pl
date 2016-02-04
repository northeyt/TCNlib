#!/usr/bin/env perl
use strict;
use warnings;
use Carp;

my $relASAMin = 10;

@ARGV or die "Please supply a residue labels file";
my $resLabelsFile = shift @ARGV;
print "Getting residue relASAs: this can take a while if your data set is large ...\n";
my @residueRelASAs = getResidueRelASAsSortedByResID($resLabelsFile);
my @residueLabels  = getResidueLabelsSortedByResID($resLabelsFile);


croak "Number of residue relASA values does not match number of residue labels"
    if @residueRelASAs != @residueLabels;

my $outFile = "relASAFiltered-min$relASAMin.residue.labels";

print "Filtering residue labels by relASA\n";
open(my $OUT, ">", $outFile) or die "Cannot open file $outFile, $!";
for (my $i = 0 ; $i < @residueRelASAs ; ++$i) {
    my $resIDAndRelASA = $residueRelASAs[$i];
    my $resIDAndLabel  = $residueLabels[$i];

    my ($firstResID, $relASA) = $resIDAndRelASA =~ /(.*):(.*)/;
    my ($secndResID, $resLab) = $resIDAndLabel  =~ /(.*):(.*)/;

    croak "Residue labels $firstResID and $secndResID do not match!"
        if $firstResID ne $secndResID;

    print {$OUT} "$resIDAndLabel\n" if $relASA > $relASAMin;
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
