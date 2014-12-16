#!/usr/bin/perl -w
# create_arff.pl --- given a output csv from Anja's feature processing steps,
# creates an .arff compatible with RF models that emulate Anja's predictor, on training
# sets that do include patch ids.
# Author: Tom Northey <zcbtfo4@acrm18>
# Created: 05 Feb 2014
# Version: 0.01

use strict;
use warnings;
use Carp;
use Getopt::Long;
use IO::CaptureOutput qw( qxx );

## INITIALIZATION ######################################################################
########################################################################################

my $csv;
my $weka_jar;
my $training_arff;
my $rem_missing_blast = 0;
my $keep_intermeds = 0;
my $standardizeOnSelf;

GetOptions( 'c=s' => \$csv,
            'w=s' => \$weka_jar,
            't=s' => \$training_arff,
            'b' => \$rem_missing_blast,
            'd' => \$keep_intermeds,
            't' => \$standardizeOnSelf
        );

Usage() if ! ($csv && $weka_jar);

# Ensure that weka is in java class-path
$ENV{'CLASSPATH'} = $weka_jar;

## MAIN ################################################################################
########################################################################################

# Create raw .arff from csv
my $raw_arff_file = 'raw.arff';

convert_CSV($csv, $raw_arff_file);

# Reorder secondary_str and intf_class attribute strings to ensure they match Anja's
# training arff, as well as ensuring that scorecons are defined as numeric
my $matching_attr_strs_arff_file = 'matching_attr_strs.arff';

format_attr_strings($raw_arff_file, $matching_attr_strs_arff_file);

# Convert nominal secondary structure attribute to four binary attributes
my $binary_arff_file = 'binary.arff';

secstruct2binary($matching_attr_strs_arff_file, $binary_arff_file);

# If specified, remove instances with no blast conscore
my $to_be_standardized = '';

if ($rem_missing_blast) {
    # First convert patch id column to a nominal attribute (otherwise filter will not work)
    my $patchid_nominal_arff = 'patchid_nominal.arff';
    
    patchid2nominal($binary_arff_file, $patchid_nominal_arff);
    
    # Remove instances missing blast conscores
    my $no_blast_filtered_arff_file = 'no_blast_filtered.arff';
    remove_noblast_instances($patchid_nominal_arff, $no_blast_filtered_arff_file);

    # Convert patch id column back to a string attribute
    my $no_blast_patchid_string_arff = 'no_blast_patchid_string.arff';

    patchid2string($no_blast_filtered_arff_file, $no_blast_patchid_string_arff);

    # Remove intermediate files
    unless ($keep_intermeds) {
        `rm $patchid_nominal_arff $no_blast_filtered_arff_file`;
    }
    
    $to_be_standardized = $no_blast_patchid_string_arff;
}
else {
    $to_be_standardized = $binary_arff_file;
}

# Standardize our arff using given training arff
my $final_file = "";

my $std_file = 'standardised.arff';
my $training_std_arff = 'training_std.arff';

if ($standardizeOnSelf) {
    # Standardize against self, rather than using another set to standardize with
    standardize_arff("", "", $to_be_standardized, $std_file);

    $final_file = $std_file;
}
elsif($training_arff) {

    standardize_arff($training_arff, $training_std_arff, $to_be_standardized, $std_file);

    $final_file = $std_file;
    
    unless ($keep_intermeds) {
        `rm $training_std_arff`;
    }
}
else {
    $final_file = $to_be_standardized;
}

# Print final file to stdout
print `cat $final_file`;
 
# Unless specified, remove all intermediate .arff files
my @to_delete = ($raw_arff_file, $matching_attr_strs_arff_file, $binary_arff_file, $to_be_standardized,
                 $final_file);

unless ($keep_intermeds) {
    foreach my $file (@to_delete) {
        if (-e $file) {
            `rm $file`;
        }
    }
}

exit;

## SUBROUTINES #########################################################################
########################################################################################

sub Usage{
    print <<EOF;
$0: -b -d -c FILE -w FILE [-t FILE -n]
 -c : input csv to convert to arff compatable with Anja's RF Classifier
 -w : weka.jar file so that various weka commands can be run
 -t : .arff file used to standardize final .arff file. This file
      must not be standardized! If ommitted, output is not standardized.
 -n : Standardize file against itself - this is for the purposes of creating a
      training set, rather than a test set. This opt will overwrite -t opt.

 -b : remove instances missing blast conscores
 -d : keep all intermediate .arff files

Given a output csv from Anja's feature processing steps, this script creates an .arff
compatible with testing on RF models that emulate Anja's RF predictor, but have been
trained on sets that do include patch ids.
The option is also given to create a filtered .arff that does not contain instances
missing blast conscores.
EOF
    exit;
    
}

sub format_attr_strings {
    my $input_arff_file = shift
        or die "reorder_attr_strings must be passed an input arff file";
    my $output_arff_file = shift
        or die "reorder_attr_strings must be passed an output arff file";
    
    open(my $IN, '<', $input_arff_file)
        or die "Cannot open file $input_arff_file to read, $!";

    my $arff_str;
    {
        local $/;
        $arff_str = <$IN>;
    }
    
    close $IN;

    # Ensure that secondary structure attribute lists classes are correctly ordered
    # to match Anja's training sets
    $arff_str =~ s/secondary_str {.*?}/secondary_str {H,EH,E,C}/;
    $arff_str =~ s/intf_class {.*?}/intf_class {I,S}/;

    # Sometimes fosta/blast scorecons are defined as string, rather than numeric ...
    # so ensure that this is corrected
    $arff_str =~ s/_scorecons string/_scorecons numeric/;
    
    open(my $OUT, '>', $output_arff_file)
        or die "Cannot open file $output_arff_file to write, $!";

    print {$OUT} $arff_str;
    close $OUT;
    
}

sub convert_CSV {
    my $input_CSV = shift or die "convert_CSV must be passed an input CSV";
    my $output_arff = shift or die "convert_CSV must be passed an output .arff";

    my $cmd = "java weka.core.converters.CSVLoader -S first $csv > $raw_arff_file";

    my $success = eval { run_WEKA($cmd); 1; };

    croak "convert_CSV failed: $@" if ! $success;
}


sub secstruct2binary {
    my $input_arff = shift
        or die "secstruct2binary must be passed an input .arff file";
    my $output_arff = shift
        or die "secstruct2binary must be passed an output .arff file";
    
    my $cmd = "java weka.filters.unsupervised.attribute.NominalToBinary -R 5 -i $input_arff -o $output_arff";

    my $success = eval { run_WEKA($cmd); 1; };
    
    croak "secstruct2binary failed: $@" if ! $success;
}

sub remove_noblast_instances {
    my $input_arff = shift
        or die "remove_noblast_instances must be passed an input .arff file";
    my $output_arff = shift
        or die "remove_noblast_instances must be passed an output .arff file";

    my $cmd = "java weka.filters.unsupervised.instance.SubsetByExpression -E \"not ismissing(ATT11)\" -i $input_arff > $output_arff";

    my $success = eval { run_WEKA($cmd); 1; };
    
    croak "remove_noblast_instances failed: $@" if ! $success;
}

sub patchid2nominal {
    my $input_arff = shift or die "patchid2nominal must be passed an input arff";
    my $output_arff = shift or die "patchid2nominal must be passed an output arff";

    my $cmd = "java weka.filters.unsupervised.attribute.StringToNominal -R first -i $input_arff > $output_arff";
    
    my $success = eval { run_WEKA($cmd); 1; };
    
    croak "patchid2nominal failed: $@" if ! $success;
}

sub patchid2string {
    my $input_arff = shift or die "patchid2string must be passed an input arff";
    my $output_arff = shift or die "patchid2string must be passed an output arff";

    my $cmd = "java weka.filters.unsupervised.attribute.NominalToString -C first -i $input_arff > $output_arff";
    
    my $success = eval { run_WEKA($cmd); 1; };
    
    croak "patchid2string failed: $@" if ! $success;
}


sub standardize_arff{
    my ($training_arff, $training_std_arff, $to_be_standardized, $std_file) = @_;

    # If a training set has been supplied, use this to standardize against.
    # Otherwise, set is standardized against itself
    my $opts
        = $training_arff ? " -b -i $training_arff -o $training_std_arff -r $to_be_standardized -s $std_file"
            : "-i $to_be_standardized -o $std_file";
    
    my $cmd = "java weka.filters.unsupervised.attribute.Standardize $opts";
    
    my $success = eval { run_WEKA($cmd); 1; };
    
    croak "remove_noblast_instances failed: $@" if ! $success;
}


sub run_WEKA {
    my $cmd = shift or die "run_WEKA must be passed a command string!";
    
    my($stdout, $stderr, $success) = qxx($cmd);
    
    if (! $success) {
        croak "run_WEKA failed -  Command run: '$cmd'\nErrors; '$stderr'";
    }

    return $stdout;
}


__END__
