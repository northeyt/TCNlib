#!/acrm/usr/local/bin/perl -w

# TCN THIS SCRIPT IS MODIFIED FROM VERSION AT ~/anya/inft_prediction/code_new. ALL MODIFICATIONS BY ME (TCN) ARE COMMENTED WITH A STARTING TCN
use strict;
use warnings;
use Carp;
#this script is old create.patches.pl (modified) from ~/intf_prediction/pqs2xmas/patch_creation
#01.07.09 changed xmas file dir from /acrm/data/xmas/pqs/ to /acrm/data/xmas/pqsatom

# TCN: added cmd-line input for specifying a text file matching 'multichain.prots.PQS.txt' format and a directory of xmas files

my $USAGE = <<EOF;
USAGE: PQS.txt xmasdir 
Please specify a text file matching 'multichain.prots.PQS.txt' format and
a directory containing Antibody-Antigen complex xmas files
EOF

croak $USAGE
    if ! @ARGV;

my $input_file = $ARGV[0] or croak $USAGE;
my $xmas_dir   = $ARGV[1] or croak $USAGE; 

####variables and handles####
# TCN: Altered file for IN file handle to $ARGV[0] from a static file name
open(IN, $ARGV[0]) || die "extract.residues: Cannot open $ARGV[0]!!!\n";
open(LOG, ">log.patchCreation.txt") || die "extract.residues: Cannot open LOG!!!\n";
open(CENT, ">patch.centres") || die "extract.residues: Cannot open CENT!!!\n";
open(SURF, ">patch.residues") || die "extract.residues: Cannot open SURF!!!\n";
open(INTF, ">intf.residues") || die "extract.residues: Cannot open INTF!!!\n";
open(CORE, ">core.residues") || die "extract.residues: Cannot open CORE!!!\n";

####   main   #### 

#1. list all relevant xmas files
#1.a) eliminate xmas files for obsolete PDB files
my %xmas_files = ();
my %empty_xmas_files = ();

print "[list_xmas_files]...";
# TCN: list_xmas_files is now sent the previously specified xmas dir
&list_xmas_files(\%xmas_files, \%empty_xmas_files, $xmas_dir);
print "done.\n";


#2. create 3 lists for all multichain files
while (my $line = <IN>)
{
    chomp $line;
    my ($file, $chainstring) = split(':', $line, 2);
    $file =~ /pqs\/(\S+)\.mmol/;
    # TCN: changed id to contain antigen chain id
    my $id = $1 . "_$chainstring";

    my $xmas_empty
        = &is_xmas_empty($id, \%xmas_files, \%empty_xmas_files );

    if ($xmas_empty == 1)
    {
        print "processing $id...";

        # Altered $xmas_filename var to include non-static $xmas_dir
        my $xmas_filename = "$xmas_dir/$id.xmas";

        my %central_residues = ();
        my %surface_residues = ();
        my %interface_residues = ();
        my %core_residues = ();

        #finds all residues with:
        #(i) rel.access. in monomer >25%
        #(ii) rel.access. in monomer >5%
        #(iii) rel.access. difference between monomer and PQS >10%
        my $chains = &FindPatchCentral($xmas_filename, $chainstring, \%central_residues, \%surface_residues, \%interface_residues, \%core_residues);
        
        my $cent_res_count = scalar keys %central_residues;
        my $surf_res_count = scalar keys %surface_residues;
        my $intf_res_count = scalar keys %interface_residues;
        my $core_res_count = scalar keys %core_residues;

        printf LOG "For %s found %d centres of patches, %d surface residues, %d core residues and %d interface residues.\n", $file, $cent_res_count, $surf_res_count, $core_res_count, $intf_res_count;
        
        print CENT "$id->(";
        foreach my $res (sort keys %central_residues)
            {
                print CENT "$res,";
            }
        print CENT ")\n";
            
        print SURF "$id->(";
        foreach my $res (sort keys %surface_residues)
            #while (my ($res, $value) = each %surface_residues)
            {
                print SURF "$res,";
            }
        print SURF ")\n";
        
        print CORE "$id->(";
        foreach my $res (sort keys %core_residues)
            {
                print CORE "$res,";
            }
        print CORE ")\n";
        
        print INTF "$id->(";
        foreach my $res (sort keys %interface_residues)
            #while (my ($res, $value) = each %interface_residues)
            {
                print INTF "$res,";
            }            
        print INTF ")\n";
    }
   
    print " done.\n";
}
#else case prints to LOG, done in &is_xmas_empty

close(IN);
close(LOG);
close(CENT);
close(SURF);
close(INTF);
close(CORE);
exit;

########################################################################
#SUBROUTINES
########################################################################
#separates xmas files into empty ones and the rest
sub list_xmas_files
{
    # TCN: list_xmas_files is now sent $xmas_dir
    my ($nonempty, $empty, $xmas_dir) = @_;

    # ls cmd now uses $xmas_dir variable
    my $ls = `ls -ltr $xmas_dir`;
    my @ls_array = split(/\n/, $ls);
    
    foreach my $xmas_file (@ls_array)
    {
        if ($xmas_file !~ /total/)
        {
            my @details = split(/\s+/, $xmas_file);
            my $size = $details[4];
            my $xmas_filename = $details[8];
            
            if ($size)
            {
                $$nonempty{$xmas_filename} = $size; 
            }
            else
            {
                $$empty{$xmas_filename} = $size; 
            }
        }
    }
} 

########################################################################
#returns 1 if xmas file is in $nonempty hash, 0 if it is in $empty hash, 999 else
sub is_xmas_empty
{
    my ($id, $nonempty, $empty, $xmas_dir) = @_;
    my $value = 999;
    # TCN: xmas_filename now contains ag-specific suffix
    my $xmas_filename = $id.'.xmas';

    if ($$nonempty{$xmas_filename})
    {
        $value = 1;
    }
    else
    {
        print LOG "$id not proccessed - empty xmas file\n";
        #print LOG "$id not proccessed - smth wrong with xmas file\n";
    }

    return $value;
}


########################################################################
#INPUT: xmas filename (here: from PQS file) 
#finds all residues with relative ASA>=25%
#OUTPUT: hash containing centres of patches, key is residue [chain]resnum[insert], value is rel. accessibility of residue _in monomer_
sub FindPatchCentral
{
    my ($xmas_file, $chain_string, $central, $surface, $interface, $core) = @_;
    my $inatoms_datatype = 0;
    my $inmolecules = 0;
    my $indata = 0;
    my $type = 0;
    my $chain = 'not defined';

    $chain_string = "\L$chain_string";
    my @array = split(':', $chain_string);
    my $chain_num = @array;	


    #open file, return a hash with all residues (key as chain, resnum, insert, rest as value)
    open(XMAS, "$xmas_file") || die "in FindPatchCentral: $xmas_file not opened, $!"; 
	
    while (my $line = <XMAS>)	
    {
        chomp $line;	        
        $line = "\L$line";
        
        # Clear flags on leaving atoms data section
        if ( ($inatoms_datatype) && ($line =~ /^<\/data/) )
        {
            $inatoms_datatype = 0;
            $inmolecules = 0;
        }
        
        if ($line =~ /^<data type=molecules>$/)
        {
            $inmolecules = 1;
        }
        
        #going through molecules tab
        if ( ($inmolecules) && (!$inatoms_datatype) )
        {
            if ($line =~ /(\S+)\s+protein\s+(\S+)\s+(\S+)\s+(\S+)/) 
            {
                my $chain_id = $3;
                my $mol_id = $1;
               
                #$chains{$chain_id} = $mol_id;
            }
        }
        
        #Check for entering data atoms section
        if ($line =~ /<data type=atoms>/)
        {
            $inatoms_datatype = 1;
        }		
            
        # Handle ATOMS data section
        if ( ($inatoms_datatype) && ($inmolecules) )
        {   
            if ($line =~ /<chain>/)
            { 	
                $line =~ /<chain>(.*)<\/chain>/;
                $chain = $1;
                $chain =~ s/\s//g;
                $chain = "" if($chain eq "|sp|");
            }
            elsif ($line =~ /<type>(.*)<\/type>/)
            {				
                $type = 0;
                $type = 1 if($line =~ /atom/);
            }
            
            if ($type)
            {
                #parsing accessibility data from residue entry
                # TCN: changed regex to deal with 8th SS column
                if ( ($chain_string =~ $chain) && 
                     ($line =~ /^<residue>\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+<\/residue>/) )                
                {
                    
                    my $resnam = $1;
                    my $resnum = $2;
                    my $molnum = $3;
                    my $a1 = $4;
                    my $a2 = $5;
                    my $a3 = $6;
                    my $a4 = $7;

                    $resnum =~ s/\.//g;
                    #$resnum =~ s/[a-z]//g;
                    my $resid = $chain.":".$resnum.":".$resnam;
                    
                    #detects centres of patches
                    if ($a4 > 25) 
                    {
                        $$central{$resid} = "$a4";                        
                    }

                    #detects residues that could constitute patches
                    if ($a4 > 5)
                    {
                        $$surface{$resid} = "$a4";
                    }
                    #non-surface residues are core...
                    else
                    {
                        $$core{$resid} = "$a4";
                    }

                    #detects residues that are in protein-protein interface
                    my $access_drop = $a4 - $a2;
                    if ($access_drop > 10)
                    {
                        $$interface{$resid} = "$access_drop";
                    }
                }
            }
        }              
    }        
    close(XMAS);

    return $chain_num;
}
########################################################################
