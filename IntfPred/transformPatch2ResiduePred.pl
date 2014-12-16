#!/usr/bin/perl -w
# transformPatch2ResiduePred.pl --- This script will transform per-patch
# predictions into per residue predictions
# Author: Tom Northey <zcbtfo4@acrm18>
# Created: 08 Dec 2014
# Version: 0.01

use warnings;
use strict;
use Data::Dumper;

my $inputCSV       = shift @ARGV;
my $patchesDir     = shift @ARGV;
my $resIDLabelFile = shift @ARGV;

### MAIN #######################################################################
################################################################################

my %resID2LabelMap = mapResID2Label($resIDLabelFile);

my @csvLines = getLinesFromCSVFile($inputCSV);
my %patchID2PredMap = mapPatchID2Pred(@csvLines);
my %patchID2ResSeqMap = mapPatchID2ResSeq($patchesDir);

my %resID2PredAndValueMap
    = resID2PredAndValue(\%patchID2PredMap, \%patchID2ResSeqMap);

print scalar (keys %resID2PredAndValueMap) . " residues found in patches\n";
print scalar (keys %resID2LabelMap) . " residues found in resID label file\n";
exit;

ammendValueLabels(\%resID2PredAndValueMap, \%resID2LabelMap);

my $CSVheader = getCSVHeader($inputCSV);

printCSV($CSVheader, \%resID2PredAndValueMap);

### SUBROUTINES ################################################################
################################################################################

sub printCSV {
    my $CSVheader = shift;
    my $resID2PredAndValueHref = shift;

    my $outputFile = "transformed.csv";
    open(my $OUT, ">", $outputFile) or die "Cannot open file $outputFile, $!";

    print {$OUT} $CSVheader;

    my $i = 1;
    foreach my $resID (keys %{$resID2PredAndValueHref}){
        my $infoHref = $resID2PredAndValueHref->{$resID};
        my $error = $infoHref->{value} eq $infoHref->{prediction} ? "" : '+';
        
        my @ordered = ($i, $infoHref->{value}, $infoHref->{prediction},
                       $error, $infoHref->{score}, $resID);

        print {$OUT} join(",", @ordered) . "\n";
    }
}

sub ammendValueLabels {
    my $resID2PredAndValueHref = shift;
    my $resID2LabelHref = shift;

    foreach my $resID (keys %{$resID2LabelHref}) {
        if (! exists $resID2PredAndValueHref->{$resID}) {
            print "TESTING: adding residue $resID\n";
            # Add resID to map
            $resID2PredAndValueHref->{$resID}
                = {prediction => "2:S",
                   value => $resID2LabelHref->{$resID},
                   score => 0.00};
        }
        else {
            print "TESTING: residue $resID has a record\n";
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
                       prediction => $patchID2PredMap->{$patchID}->{prediction},
                       score => $patchID2PredMap->{$patchID}->{score}};
            }
            elsif ($resID2PredAndValue{$resID}->{prediction} eq "2:S") {
                # Only ammend prediction if resID has not yet been labelled
                # positive.
                $resID2PredAndValue{$resID}->{prediction}
                    = $patchID2PredMap->{$patchID}->{prediction};
            }
        }
    }
    return %resID2PredAndValue;
}

sub getCSVHeader {
    my $inputCSV = shift;

    open(my $CSV, "<", $inputCSV) or die "Cannot open file $inputCSV, $!";
    
    while (my $line = <$CSV>) {
        if ($line =~ /^inst#/){
            return $line;
        }
    }
    close $CSV;

    die "Did not parse header from CSV file!";
}

sub getLinesFromCSVFile {
    my $file = shift;

    open(my $IN, "<", $file) or die "Cannot open file $file, $!";

    my @lines = ();

    my $reachedHeader = 0;
    
    while (my $line = <$IN>) {
        if (! $reachedHeader) {
            $reachedHeader = 1 if $line =~ /^inst#/;
            next;
        }
        next if $line =~ /^\n$/;

        push(@lines, $line);
    }
    return @lines;
}

sub parseCSVLine {
    my $line = shift;
    chomp $line;
    my ($inst, $value, $prediction, $err, $score, $patchID) = split(",", $line);

    return [$patchID, $value, $prediction, $score];
}

sub mapPatchID2Pred {
    my @csvLines = @_;

    my %map = ();
    
    foreach my $infoAref (map {parseCSVLine($_)} @csvLines) {
        my ($patchID, $value, $prediction, $score) = @{$infoAref}; 
        $map{$patchID} = {value => $value, prediction => $prediction,
                          score => $score};
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

