#!/usr/bin/perl -w
# calcPatchClassStats.pl --- This script creates AA intf/surface propensity
# files from a patches directory and a .labels file. This script is required
# when a pre-existing patches directory is supplied to createDataSet.pl
# Author: Tom Northey <zcbtfo4@acrm18>
# Created: 12 Nov 2014
# Version: 0.01

use warnings;
use strict;
use pdb::pdb;
use pdb::get_files;
use Carp;
use File::Basename;

### Cmd line processing + Init #################################################
################################################################################

Usage() unless @ARGV == 4;

my $patchesDir = shift @ARGV;
die "$patchesDir is not a directory!" if ! -d $patchesDir;

my $labelsFile = shift @ARGV;

my %patchLabelsHash = processLabelsFile($labelsFile);

my $intfStatsOutputFile = shift @ARGV;
my $surfStatsOutputfile = shift @ARGV;

my ($intfStatHref, $surfStatHref) = getIntfAndSurfStatHrefs();

### Main #######################################################################
################################################################################

foreach my $patchFile (getPatchFiles($patchesDir)) {
    my ($pdb, $chainID) = getPDBAndChainID($patchFile);

    print "Processing patches for pdbid $pdb$chainID\n";
    
    my @patchSummaryLines = getPatchSummaryLines($patchFile);
    
    my $chain = chain->new(pdb_code => $pdb,
                           chain_id => $chainID,
                           hydrogen_cleanup => 1);
    
    $chain->read_ASA();
    
    foreach my $patchLine (@patchSummaryLines) {

        my $patch = patch->new(summary => $patchLine,
                               parent_pdb => $chain,);

        # Get patch label
        my $key = join(":", ($pdb, $chainID, $patch->central_atom->resSeq()));

        print "TEST: $key\n";
        
        croak "Could not find label for patch: " . $patch->id() . "\n"
            unless exists $patchLabelsHash{$key};
        
        my $patchLabel = $patchLabelsHash{$key};
        
        foreach my $resid (keys %{$patch->resid_index()}) {
            my @residAtoms = values %{$patch->resid_index->{$resid}};
            
            my $statHref = $patchLabel eq 'I' ? $intfStatHref : $surfStatHref;

            my $resid2ndForm
                = join(".", ($chain->chain_id(), $residAtoms[0]->resSeq()));

            if (! exists $chain->resid2RelASAHref->{$resid2ndForm}) {
                print "No ASA for $resid2ndForm!\n";
                next;
            }
            elsif ($chain->resid2RelASAHref->{$resid2ndForm} < 5){
                print "Skiping residue, low ASA ...\n";
                next;
            }
            
            
            my $resName = $residAtoms[0]->resName();

            # Skip non-natural amino acids like MSE
            next if ! exists $statHref->{$resName};
            
            # Add residue ASAb to sum for residue type
            map {$statHref->{$resName}->{sum} += $_->ASAm()} @residAtoms;
            #$statHref->{$resName}->{sum} += $patch->total_ASAm();

            # Increment residue type count
            ++$statHref->{$resName}->{count};
        }
    }
}


printStatsToOut($patchesDir, $labelsFile,
                $intfStatsOutputFile, $intfStatHref);
printStatsToOut($patchesDir, $labelsFile,
                $surfStatsOutputfile, $surfStatHref);

print "Finished!";

### Subroutines ################################################################
################################################################################

sub printStatsToOut {
    
    my $patchesDir = shift;
    my $labelsFile = shift;
    
    my $outputFile = shift;
    my $statHref = shift;
    
    open(my $OUT, ">", $outputFile) or die "Cannot open file $outputFile, $!";

    my $header
        = "input is patches dir $patchesDir with labels file $labelsFile\n";

    print {$OUT} $header;

    my $linePrefix = "[type:sum:count:mean]";
    
    foreach my $resType (keys %{$statHref}) {
        
        my $total = $statHref->{$resType}->{sum};
        my $count = $statHref->{$resType}->{count};

        my $mean = $total / $count;
        
        printf {$OUT} $linePrefix . "%s:%12.2f:%6d:%10.4f\n",
            (lc($resType), $total, $count, $mean);
    }

    my $ASAtotal = calcASATotalFromStatHref($statHref);
    printf {$OUT} "total ASA for the dataset is %.4f.\n", $ASAtotal;
}

sub calcASATotalFromStatHref {
    my $statHref = shift;

    my $total = 0;

    map {$total += $_->{sum}} values %{$statHref};

    return $total;
}


sub getResPropensityHref {
    my @AAs = qw(ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG
                 SER THR VAL TRP TYR);

    my %intfStatsHash = map {$_ => {sum => 0, count => 0, mean => 0}} @AAs;
    my %surfStatsHash = %intfStatsHash;

    return (\%intfStatsHash, \%surfStatsHash);
}


sub processLabelsFile {
    my $labelsFile = shift;

    my %labelsHash = ();
    
    open(my $IN, "<", $labelsFile) or die "Cannot open file $labelsFile, $!";

    while (my $line = <$IN>) {
        # Example line format: 2nz9:A:622:S
        chomp $line;

        my ($pdbCode, $chainID, $resSeq, $label) = split(":", $line);

        my $key = join(":", ($pdbCode, $chainID, $resSeq));

        $labelsHash{$key} = $label;
    }
    return %labelsHash;
}


sub getPatchSummaryLines {
    my $file = shift;

    open(my $IN, "<", $file) or die "Cannot open file $file, $!";

    my @patchSummaryLines = <$IN>;

    return @patchSummaryLines;
}


sub getPDBAndChainID {
    my $fileName = shift;

    my ($name, $path, $suffix) = fileparse($fileName, ".patches");
                                           
    my $pdb = lc substr($name, 0, 4);
    my $chainID = uc substr($name, 4, 1);

    return ($pdb, $chainID);
}


sub getPatchFiles {
    my $patchesDir = shift;

    $patchesDir .= '/' unless $patchesDir =~ /\/$/;
    
    opendir(my $DH, $patchesDir) or die "Cannot open directory $patchesDir, $!";

    my @patchFiles = ();
    
    while (my $file = readdir($DH)) {
        next unless $file =~ /\.patches$/;
        
        push(@patchFiles, $patchesDir . $file);
    }

    die "No patch files were parsed from $patchesDir!" if ! @patchFiles;
    
    return @patchFiles;
}


sub Usage {
    print <<EOF;
$0 USAGE: patchesDir labelsFile intfStatsOutputFile surfStatsOutputfile

This script creates AA intf/surface propensity files from a patches directory
and a .labels file. This script is required when a pre-existing patches
directory is supplied to createDataSet.pl
EOF
    exit(1);
}

