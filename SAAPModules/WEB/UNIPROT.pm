package WEB::UNIPROT;
#*************************************************************************
#
#   Program:    
#   File:       
#   
#   Version:    
#   Date:       
#   Function:   
#   
#   Copyright:  (c) Dr. Andrew C. R. Martin, UCL, 2011
#   Author:     Dr. Andrew C. R. Martin
#   Address:    Institute of Structural and Molecular Biology
#               Division of Biosciences
#               University College
#               Gower Street
#               London
#               WC1E 6BT
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
#   Web services-based routines to access FOSTA and UniProt FASTA files
#
#*************************************************************************
#
#   Usage:
#   ======
#   $sequence = GetFASTAWeb($swissprotID);
#      Extracts the sequence for a given SwissProt ID - web services
#      access to Uniprot
#
#   CONFIG:
#   =======
#   Requires a config.pm module which specifies...
#     $LocalSwissProt - 1: use local files, 2: use web
#   If using local files, the following must be specified
#     $swissprot      - SwissProt file if using local files
#     $sprotCacheDir  - Directory for SwissProt index
#     $sprotIndex     - Index file (full path)
#     $getsprot       - getsprot indexing program
#     $indexsprot     - SwissProt indexing program
#
#*************************************************************************
#
#   Revision History:
#   =================
#
#*************************************************************************
use strict;
use WEB;

#-------------------------------------------------------------------------
my $uniprotURLfasta = "http://www.uniprot.org/uniprot/%s.fasta";
my $uniprotURLsprot = "http://www.uniprot.org/uniprot/%s.txt";

#-------------------------------------------------------------------------
sub GetFASTA
{
    my ($id) = @_;

    # Grab the FASTA sequence from the UniProt web service
    my $url = sprintf($uniprotURLfasta, $id);
    my $ua = WEB::CreateUserAgent("");
    my $req = WEB::CreateGetRequest($url);
    my $content = WEB::GetContent($ua, $req);

    return($content);
}

sub GetSwissProt
{
    my($id) = @_;

    # Grab the SwissProt entry from the UniProt web service
    my $url = sprintf($uniprotURLsprot, $id);
    my $ua = WEB::CreateUserAgent("");
    my $req = WEB::CreateGetRequest($url);
    my $content = WEB::GetContent($ua, $req);

    return($content);
}

1;
