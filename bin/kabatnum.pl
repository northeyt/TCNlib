#!/acrm/usr/local/bin/perl
#*************************************************************************
#
#   Program:    KabatNum
#   File:       kabatnum.pl
#   
#   Version:    V2.0
#   Date:       18.04.08
#   Function:   Apply standard Kabat numbering to a PDB file
#   
#   Copyright:  (c) Dr. Andrew C. R. Martin 1995-2008
#   Author:     Dr. Andrew C. R. Martin
#   Address:    Biomolecular Structure & Modelling Unit,
#               Department of Biochemistry & Molecular Biology,
#               University College,
#               Gower Street,
#               London.
#               WC1E 6BT.
#   EMail:      andrew@bioinf.org.uk
#               
#*************************************************************************
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If 
#   someone else breaks this code, I don't want to be blamed for code 
#   that does not work! 
#
#   The code may not be sold commercially or included as part of a 
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************
#
#   Description:
#   ============
#
#*************************************************************************
#
#   Usage:
#   ======
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0  09.08.95 Original    By: ACRM
#   V1.1  10.11.95 Allowed -c flag for Chothia numbering
#   V2.0  18.04.08 Modified to use Abhi's numbering program and generally
#                  cleaned up.
#
#*************************************************************************
use strict;
use warnings;

# Path for other programs we use
my $bin = "/home/bsm/martin/bin";
my $abhibin = "/home/bsm/martin/abnum/installed/numbering";

# System programs we use
my $cp  = "/bin/cp";
my $cat = "/bin/cat";

# And specify the names of these programs
my $chainpdb    = "$bin/chainpdb";
my $pdb2pir     = "$bin/pdb2pir";
my $patchpdbnum = "$bin/patchpdbnum";
my $kabatseq    = "$abhibin/kabnum_wrapper.pl";

# Turn on flushing
$| = 1;

# Initialise flags and variables
my $flags   = " ";
my $outfile = " ";
my $relabel = 0;

# Check for flags
while(substr($ARGV[0],0,1) eq "-")
{
    if($ARGV[0] eq "-c")
    {
        $flags = " -c ";
    }
    elsif($ARGV[0] eq "-k")
    {
        $flags = " -k ";
    }
    elsif($ARGV[0] eq "-a")
    {
        $flags = " -a ";
    }
    elsif($ARGV[0] eq "-l")
    {
        $relabel = 1;
    }
    elsif($ARGV[0] eq "-h")
    {
        Usage();
        exit 0;
    }
    else
    {
        Usage();
        exit 1;
    }
    shift @ARGV;
}

# Check the input PDB file has been specified 
if($#ARGV < 0)
{
    Usage();
    exit 0;
}

# Record the input file and output file if specified
my $infile = $ARGV[0];
if($#ARGV == 1)
{
    $outfile = $ARGV[1];
}

# Create names for temporary files
my $patchfile = "/tmp/seq.$$";
my $pdbfile   = "/tmp/pdb.$$";
my $resfile   = "/tmp/res.$$";
my $pirfile   = "/tmp/pir.$$";

# Create a temp file containing the Kabat numbered sequence
`$pdb2pir $infile >$pirfile`;
`$kabatseq $pirfile $flags >$patchfile`;

# Run the Patch program on the input file using the sequence file
system("$patchpdbnum $patchfile $infile $pdbfile");

# If the -l flag was specified, run this through chainpdb
if($relabel)
{
    if($outfile eq " ")
    {
        system("$chainpdb $pdbfile");
    }
    else
    {
        system("$chainpdb $pdbfile $outfile");
    }
}
else    # Not -l
{
    if($outfile eq " ")
    {
        system("$cat $pdbfile");
    }
    else
    {
        system("$cp $pdbfile $outfile");
    }
}

# Clean up by removing temp files
unlink($patchfile);
unlink($pdbfile);
unlink($resfile);
unlink($pirfile);
        
exit 0;




#*************************************************************************
#  sub Usage()
#  -----------
#  Print a usage message
#
#  09.08.95 Original    By: ACRM
#  10.11.95 V1.1, Added -c
#  18.04.08 V2.0 Rewritten to use Abhi's numbering program
#
sub Usage
{
    print <<__EOF;

KabatNum V2.0 (c) 1995-2008 Dr. Andrew C.R. Martin, UCL

Usage: kabatnum [-k|-c|-a] [-l] in.pdb [out.pdb]
        -k Do Kabat numbering (default)
        -c Do Chothia numbering rather than Kabat
        -a Do improved Chothia numbering rather than Kabat
        -l Relabel chains as A and B

Applies standard numbering to a PDB file

Output to stdout if file not specified. (Sorry this can't be piped
into at the moment.)

If the PDB file contains non-antibody chains these will not appear in
the output.

If the file contains more than one antibody, each will be labelled
with L and H unless -l is specified in which case the chains will be
labelled A,B,C,D, etc.

__EOF
}
