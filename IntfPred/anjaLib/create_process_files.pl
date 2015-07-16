#!/usr/bin/perl -w
# create_process_files.pl --- creates required files for running Anja's WEKA classifier models
# Author: Tom Northey <zcbtfo4@acrm17>
# Created: 03 Jan 2014
# Version: 0.01

use warnings;
use strict;
use Carp;
use File::chdir;
use File::Spec;
use Getopt::Long;

use pdb::get_files;
use TCNPerlVars;

use parse_how_file qw( parse_how_file );

use IO::CaptureOutput qw( capture_exec );

my $USAGE = <<EOF;
$0 USAGE :  [-how_file -chain_file FILE] -pdb_dir DIR -xmas_dir DIR

This script creates pdb and xmas files modified to contain only those
chains specified in the input .how or .chain file. These pdb and xmas
files are used to correctly calculate ASA values for the chains specified.
This script also creates a prots.PQS.mimic files that is required for the
rest of the process.

A .how file can be specified as input (see anja_predictor/datasets for example),
or a .chain file. Chain files are a whitespace-separated list of antigen
chains, e.g.
    4houA 1afvB 2dfgDC

Note that multi-chain inputs  are specified with their multiple chain IDs.
If a .chain file is specified, note that their must also have been a
class labels file supplied to the overall process, if real epitope labels are
required (see run_from_how.pl).

EOF

@ARGV or die $USAGE;

my $how_file;
my $chain_file;
my $pdb_dir;
my $xmas_dir;

GetOptions(
    'xmas_dir=s' => \$xmas_dir,
    'pdb_dir=s'  => \$pdb_dir,
    'how_file=s' => \$how_file,
    'chain_file=s' => \$chain_file,
);

# Files and directories for creating xmas files
my $pdb2xmas = $TCNPerlVars::pdb2xmas;
my $required_bin
    = "/home/bsm/martin/acrm/CONSULTANCY/inpharmatica/software/bin/";
my $solv = $required_bin . 'solv';
my $ss   = $required_bin . 'ss';
my $hb   = $required_bin . 'hb';

# Open prots.PQS.mimic for output
my $output_file = 'prots.PQS.mimic';

open( my $OUT, '>', $output_file )
    or die "Cannot open output file '$output_file', $!";

if ($how_file) {
    createProcessFilesFromHowFile($how_file, $OUT,
                                  $pdb_dir, $xmas_dir);
}
elsif ($chain_file) {
    createProcessFilesFromChainFile($chain_file, $OUT,
                                    $pdb_dir, $xmas_dir);
}
else {
    croak "You must specify either a .how or .chain file!";
}

close $OUT;
print "create_process_files.pl finished\n";
exit;

### SUBROUTINES ################################################################
################################################################################

sub createProcessFilesFromChainFile {
    my $chain_file = shift;
    my $OUT = shift;
    my $pdb_dir = shift;
    my $xmas_dir = shift;

    my @pdbAndChainIDs = pdbIDsFromChainFile($chain_file);

    foreach my $arrRef (@pdbAndChainIDs) {
        my $pdb_code = $arrRef->[0];
        my @chainIDs = @{$arrRef->[1]};

        # parse pdb lines for those starting with relevant chains
        my $get_files = pdb::get_files->new(pdb_code => $pdb_code);
        
        my $pdb_file = eval {$get_files->pdb_file()};
        
        if (! $pdb_file) {
            print $@;
            next;
        }

        open(my $fh, '<', $pdb_file ) or die "Cannot open '$pdb_file', $!";

        my $complex_name = $pdb_code . "_" . join("", @chainIDs);
        
        createPDBAndXMASFiles($pdb_file, \@chainIDs, $pdb_dir, $xmas_dir,
                              $complex_name);
        
        # add file to prots.PQS.mimic
        my $chains_string = join(':', @chainIDs);
        
        my $line = 'pqs/' . uc ( $pdb_code ) . ".mmol:$chains_string\n";
        
        print {$OUT} $line;
    }
}

sub createProcessFilesFromHowFile {
    my $how_file = shift;
    my $OUT = shift;
    my $pdb_dir = shift;
    my $xmas_dir = shift;

    my @all_info = parse_how_file($how_file);

    foreach my $complex_hr (@all_info) {

        my $pdb_code = $complex_hr->{pdb_code};
        
        my @ab_chain_ids = map { uc $_ } @{ $complex_hr->{antibody_chain_ids} }
            or croak "No antibody chain ids parsed for pdb $pdb_code";
        
        my $ag_chain_id = uc $complex_hr->{antigen_chain_id}
            or croak "No antigen chain id parsed for pdb $pdb_code";
        
        my $complex_name = "$pdb_code" . "_$ag_chain_id";
        
        print "pdb: $pdb_code antigen chain: $ag_chain_id, "
            . "antibody chains: @ab_chain_ids\n";
        
        # parse pdb lines for those starting with relevant chains
        my $get_files = pdb::get_files->new(pdb_code => $pdb_code);
        
        my $pdb_file = eval {$get_files->pdb_file()};
        
        if (! $pdb_file) {
            print $@;
            next;
        }
                    
        my @chains = (@ab_chain_ids, $ag_chain_id);

        createPDBAndXMASFiles($pdb_file, \@chains, $pdb_dir, $xmas_dir,
                              $complex_name);
            
        # add file to prots.PQS.mimic
        my $chains_string = join( ':', @chains );
        
        my $line = 'pqs/' . uc ( $pdb_code ) . ".mmol:$ag_chain_id\n";
        
        print {$OUT} $line;
    }
}

sub createPDBAndXMASFiles {
    my $pdb_file = shift;
    my $chainIDAref = shift;
    my $pdb_dir = shift;
    my $xmas_dir = shift;
    my $complex_name = shift;
    
    my @chains = @{$chainIDAref};

    # Test that pdb file contains ATOM lines for all chains - otherwise
    # complex is invalid
    
    my %test_chain_lines_exist = map { $_ => 0 } @chains;

    my @chain_lines = ();
    
    open(my $fh, '<', $pdb_file ) or die "Cannot open '$pdb_file', $!";

    while (my $line = <$fh>) {
        if ($line =~ /^ATOM/ || $line =~ /^HETATM/) {
            my $chain_id = substr ($line, 21, 1);
            
            if (grep /^$chain_id$/, @chains){
                
                push(@chain_lines, $line);
                
                $test_chain_lines_exist{$chain_id} = 1;
            }
        }
        else {
            # Keep all other data line types in case used downstream by
            # AM's scripts
            push(@chain_lines, $line);
        }
    }
    
    my $invalid_complex = 0;
    
    foreach my $chain_id (keys %test_chain_lines_exist) {
        if (! $test_chain_lines_exist{$chain_id}){
            print "\tcomplex contains invalid chain $chain_id: "
                . "no ATOM lines parsed containing chain id $chain_id\n";
            $invalid_complex = 1;
        }
    }
    
    if ($invalid_complex){
        print "\tInvalid complex, skipping ...\n";
        next;
        }
    
    # write parsed lines to pdb file
    my $mod_pdb_file = "$pdb_dir/$complex_name.pdb";
    
    open(my $PDB, '>', $mod_pdb_file)
            or die "Cannot open output file '$pdb_file, $!";
    
    print {$PDB} @chain_lines;
    close $PDB;
    
    my @output = ();
        
    # run pdb file through pdb2xmas
    # get abs pdb path before changing dir
    my $abs_pdb_path = File::Spec->rel2abs($mod_pdb_file);
    {
        # change working dir to bin to access required data files
        local $CWD = $required_bin;
        
        my $cmd = "$pdb2xmas $abs_pdb_path | $solv | $ss | $hb";
        
        my ($stdout, $stderr) = capture_exec($cmd);
        
        croak "pdb2xmas produced no output. Command run:\n$cmd\n"
            if ! $stdout;
        
        print $stderr;
        
        @output = split (/(?<=\n)/, $stdout);
    }
    
    # write output to file within xmas directory
    my $xmas_file = "$xmas_dir/$complex_name.xmas";
    
    open(my $XMAS, '>', $xmas_file)
        or die "Cannot open output file '$xmas_file', $!";
    
    print {$XMAS} @output;
    close $XMAS;
}

sub pdbIDsFromChainFile {
    my $antigen_file = shift;

    my $str;
    
    {
        local $/;
        open(my $IN, "<", $antigen_file)
            or die "Cannot open file $antigen_file, $!";
        $str = <$IN>;
    }

    my @pdbIDs = split(/\s+/, $str);
    
    my @pdbAndChainIDs
        = map { [uc(substr($_, 0, 4)), [split("", substr($_, 4))]] } @pdbIDs;

    return @pdbAndChainIDs;
}

# Run PQS.mimic.txt through modified new.extract.intf.surf.residues
    
__END__

=head1 NAME

create_process_files.pl - Describe the usage of script briefly

=head1 SYNOPSIS

create_process_files.pl [options] args

      -opt --long      Option description

=head1 DESCRIPTION

Stub documentation for create_process_files.pl, 

=head1 AUTHOR

Tom Northey, E<lt>zcbtfo4@acrm17E<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014 by Tom Northey

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
