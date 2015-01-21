#!/usr/bin/perl -w
# transformPatch2ResiduePred.pl --- This script will transform per-patch
# predictions into per residue predictions
# Author: Tom Northey <zcbtfo4@acrm18>
# Created: 08 Dec 2014
# Version: 0.01

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use Carp;
use Intf::Pred::lib::editOutputCSV;

my $minimal;
my $vote;
my $scoreAverage;
my $scoreThresh;

GetOptions("m"   => \$minimal,
           "v"   => \$vote,
           "s"   => \$scoreAverage,
           "t=f" => \$scoreThresh);

croak "Score threshold must be between 0 and 1"
    if $scoreThresh < 0 || $scoreThresh > 1;

$scoreThresh = 0.5 if ! $scoreThresh;

my $option
    = $minimal ? 0 :
    $vote ? 1 :
    $scoreAverage ? 2 : 2; # Default = 2

my $inputCSV       = shift @ARGV;
my $patchesDir     = shift @ARGV;
my $resIDLabelFile = shift @ARGV;

### MAIN #######################################################################
################################################################################

my %resID2LabelMap = mapResID2Label($resIDLabelFile);

my @csvLines = Intf::Pred::lib::editOutputCSV::getLinesFromCSVFile($inputCSV);
my %patchID2PredMap = mapPatchID2Pred(@csvLines);
my %patchID2ResSeqMap = mapPatchID2ResSeq($patchesDir);

my %resID2PredAndValueMap
    = resID2PredAndValue(\%patchID2PredMap, \%patchID2ResSeqMap, $option,
                         $scoreThresh);

print scalar (keys %resID2PredAndValueMap) . " residues found in patches\n";
print scalar (keys %resID2LabelMap) . " residues found in resID label file\n";

ammendValueLabels(\%resID2PredAndValueMap, \%resID2LabelMap);

my $CSVheader = Intf::Pred::lib::editOutputCSV::getCSVHeader($inputCSV);
# In header, replace patchID with residueID
$CSVheader =~ s/patchID/residueID/;

printCSV($CSVheader, \%resID2PredAndValueMap);

### SUBROUTINES ################################################################
################################################################################

sub printCSV {
    my $CSVheader = shift;
    my $resID2PredAndValueHref = shift;

    print $CSVheader;

    my $i = 1;
    foreach my $resID (keys %{$resID2PredAndValueHref}){
        my $infoHref = $resID2PredAndValueHref->{$resID};
        my $error = $infoHref->{value} eq $infoHref->{prediction} ? "" : '+';
        
        my @ordered = ($i, $infoHref->{value}, $infoHref->{prediction},
                       $error, $infoHref->{score}, $resID);

        print join(",", @ordered) . "\n";
    }
}

sub ammendValueLabels {
    my $resID2PredAndValueHref = shift;
    my $resID2LabelHref = shift;

    foreach my $resID (keys %{$resID2LabelHref}) {
        if (! exists $resID2PredAndValueHref->{$resID}) {
            # Add resID to map
            $resID2PredAndValueHref->{$resID}
                = {prediction => "2:S",
                   value => $resID2LabelHref->{$resID},
                   score => 0.00};
        }
        else {
            $resID2PredAndValueHref->{$resID}->{value}
                = $resID2LabelHref->{$resID};
        }
    }
}

sub mapResID2Label {

    my $resIDLabelFile = shift;
    
    open(my $IN, "<", $resIDLabelFile)
        or die "Cannot open file $resIDLabelFile, $!";

    my @lines = <$IN>;
    close $IN;
    
    my %labelConversion = (0 => '2:S', 1 => '1:I');
    
    my %map = ();

    foreach my $line (@lines) {
        chomp $line;
        my ($resID, $label) = split(":", $line);
        $map{$resID} = $labelConversion{$label};
    }
    return %map;
}

sub resID2PredAndValue {
    my $patchID2PredMap    = shift;
    my $patchID2ResSeqMap  = shift;
    my $option = shift;
    my $scoreThresh = shift;
    
    my %resID2PredAndValue = ();

    foreach my $patchID (keys %{$patchID2PredMap}) {
        
        my @patchIDSplit = split(":", $patchID);
        my $chainID = $patchIDSplit[0] . $patchIDSplit[1];
        
        my @resSeqs = @{$patchID2ResSeqMap->{$patchID}};
        
        foreach my $resSeq (@resSeqs) {
            my $resID = $chainID . $resSeq;

            if (! exists $resID2PredAndValue{$resID}) {
                $resID2PredAndValue{$resID}
                    = {value => $patchID2PredMap->{$patchID}->{value},
                       prediction => [$patchID2PredMap->{$patchID}->{prediction}],
                       score => [$patchID2PredMap->{$patchID}->{score}]};
            }
            else {
                push(@{$resID2PredAndValue{$resID}->{prediction}},
                         $patchID2PredMap->{$patchID}->{prediction});
                
                push(@{$resID2PredAndValue{$resID}->{score}},
                     $patchID2PredMap->{$patchID}->{score});
            }
        }
    }

    # Process predictions and scores
    averageScores(\%resID2PredAndValue);
    decidePredictionLabel(\%resID2PredAndValue, $option, $scoreThresh);
        
    return %resID2PredAndValue;
}

sub decidePredictionLabel {
    my $resID2PredAndValueHref = shift;
    my $option = shift;
    my $scoreThresh = shift;

    $scoreThresh = 0.5 if ! $scoreThresh;
    
    foreach my $resID (keys %{$resID2PredAndValueHref}) {
        my %label2Freq = ();
        map {++$label2Freq{$_}}
            @{$resID2PredAndValueHref->{$resID}->{prediction}};
        
        if ($option == 0) {
            $resID2PredAndValueHref->{$resID}->{prediction}
                = exists $label2Freq{"1:I"} ? "1:I" : "2:S";
        }
        elsif ($option == 1) {
            $resID2PredAndValueHref->{$resID}->{prediction}
                = exists $label2Freq{"1:I"} && $label2Freq{"1:I"} >= $label2Freq{"2:S"} ? "1:I"
                    : "2:S";
        }
        elsif ($option == 2) {
            $resID2PredAndValueHref->{$resID}->{prediction}
                = $resID2PredAndValueHref->{$resID}->{score} >= $scoreThresh ?
                    "1:I" : "2:S";
        }
    }
}

sub averageScores {
    my $resID2PredAndValueHref = shift;

    foreach my $resID (keys %{$resID2PredAndValueHref}) {
        my @scores = @{$resID2PredAndValueHref->{$resID}->{score}};
        my $sum;
        map {$sum += $_} @scores;
        my $avg = $sum / scalar @scores;
        $resID2PredAndValueHref->{$resID}->{score} = $avg;
    }
}

sub mapPatchID2Pred {
    my @csvLines = @_;

    my %map = ();

    my @infoArefs
        = map {Intf::Pred::lib::editOutputCSV::parseCSVLine($_)} @csvLines;
    
    foreach my $infoAref (@infoArefs) {
        my ($patchID, $value, $prediction, $score) = @{$infoAref}; 
        $map{$patchID} = {value => $value, prediction => $prediction,
                          score => $score};

        # Transform prediction score so that it ranges from 0 - 1,
        # where >= 0.5 = Interface
        if ($map{$patchID}->{prediction} eq "2:S") {
            $map{$patchID}->{score}
                = 1 - $map{$patchID}->{score};
        }
            
    }
    return %map;
}

sub mapPatchID2ResSeq {
    my $patchesDir = shift;

    opendir(my $DH, $patchesDir) or die "Cannot open dir $patchesDir, $!";

    my %map = ();
    
    while (my $fileName = readdir($DH)) {
        next if $fileName =~ /^\./;

        my $pdbCode = substr($fileName, 0, 4);
        my $chainID = substr($fileName, 4, 1);

        my $filePath = "$patchesDir/$fileName";
        open(my $IN, "<", $filePath) or die "Cannot open file $filePath, $!";

        while (my $line = <$IN>) {
            my ($patchCentre, @resSeqs) = parsePatchLine($line);
            my $patchID = join(":", (lc($pdbCode), $chainID, $patchCentre));
            $map{$patchID} = \@resSeqs;
        }
    }
    return %map;
}

sub parsePatchLine {
    my $line = shift;

    # example line: <patch A.24> A:23 A:24 A:25 A:26 A:27 A:28
    my @resids = $line =~ /[:.](\w+)/g;

    return @resids;
}

__END__

