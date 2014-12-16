#!/acrm/usr/local/bin/perl -w

use strict;
use Getopt::Long;  
use lib '/home/bsm/anya/perllib'; #sets dir where Perl_paths.pm is
use Perl_paths; 
 
#### Define command line inputs ####
my $patch_dir;
my $xmas_dir;
my $out_dir;
my $log_dir;     
my $help = '';
my $v_flag = '';


#### Checking command line options ####
#in is dir with <pqs_id>.patches files, out is outfile 
GetOptions('patch_dir=s' => \$patch_dir,
           'xmas_dir=s' => \$xmas_dir,
           'log_dir=s' => \$log_dir,
           'out_dir=s' => \$out_dir, 
           'help' => \$help,
           'v_flag' => \$v_flag);
#i.e. /acrm/usr/local/bin/perl calc_patch_hydrophobicity.pl -patch_dir ~/intf_prediction/datasets_new/my_dataset/patches_11/ -log_dir log_files/ -out_dir hydrophobicity/ -v_flag

if ( (!$patch_dir) || (!$log_dir) || 
     (!$out_dir) || (!$xmas_dir) )
{
    $help = 1;
}

if ($help)
{
    &Usage();
    exit;
}

##################################
#this code will produce a file for every patch-containing pqs_id
#25.02.2010. files are in /home/bsm/anya/intf_prediction/datasets_new/my_dataset/patches_11/
# for every patch in that pqs_id it will output a line:
#    -<patch A.105> X
#     where A.105 is the central residue for that patch
#     X is an average hydrophobicity of that patch (based on Kyle & Doolittle scale), calculated as sum of hydrophobicity scores of residues H(res) in that patch, divided by the #residues in the patch


####   variables and handles   ####
###################################
my $Rfile = "patches_hydrophobicities_r11.Rdata";
open (R, ">$Rfile") || die "calc_patch_hydrophobicity.pl cannot open R!!!\n";


#### main ####
##############
my @filelist = `ls $patch_dir`;

#my $file = "1bdgA.patches";
foreach my $file (@filelist)
{
    chomp $file;
    $file =~ /(.*)\.patches/;
    my $pqs_id = $1;
    my $chain = chop($pqs_id);

    print "$file\n" if ($v_flag);

    my $logfile = "$log_dir/$pqs_id.$chain.hydrophobicity.log";
    open (LOG, ">$logfile") || die "intf_sorting_patches.pl cannot open LOG!!!\n";
    my $outfile = "$out_dir/$pqs_id.$chain.hydr.patches"; 
    open (OUT, ">$outfile") || die "intf_sorting_patches.pl cannot open OUT!!!\n";

    #key is $chain:$resnum, value is aa_type:absASA
    #[residues] B:443A->ARG:92.207
    my %residues = &res_types_from_xmas($pqs_id, $chain, $xmas_dir);

    my %patches = &read_patches($file, $patch_dir, \%residues);

    foreach my $k (sort keys %patches)
    {
        #$A is area of that patch
        printf OUT "<patch %s> %.4f\n", $k, $patches{$k};
               
        #the graph is in ~/R/patches/patches_hydrophobicities_r11.Rdata
        if ($patches{$k})
        {
            printf R "$pqs_id:$k\t%.2f\n", $patches{$k};
        }
    }    
    close(OUT);
    close(LOG);
    
    if(&file_empty($logfile))
    {
        print "rm -f $logfile\n";
        `rm -f $logfile`;
    }
    else
    {
        print "$logfile had errors!!!\n";
    }   
}

close(R);
exit;


##############################################################
# SUBROUTINES
##############################################################
#read in a pqsatoms XMAS file
#outputs hash of residues (key is uc(chain):resnum), value is res_type:absASA
sub res_types_from_xmas
{
    my ($pqs_id, $chain, $xmas_dir) = @_;
    my %aatypes = ();

    my $inatoms_datatype = 0;
    my $inmolecules = 0;
    my $type = 0;
    my $xmas_chain;

#    my $xmas_file = "/grid/gridstore/xmas/pqsatoms/pqs".$pqs_id.".xmas";
    my $xmas_file = "$xmas_dir/" . uc ("$pqs_id" . "_$chain") . ".xmas";
    open (XMAS, "$xmas_file") || die "sub res_types_from_xmas cannot open $xmas_file!!!\n";
        
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
                $xmas_chain = $1;
                $xmas_chain =~ s/\s//g;
                $xmas_chain = "" if($xmas_chain eq "|sp|");
            }
            elsif ($line =~ /<type>(.*)<\/type>/)
            {				
                $type = 0;
                $type = 1 if($line =~ /atom/);
            }
            
            if ($type)
            {
                #parsing accessibility data from residue entry
                if ($line =~ /^<residue>\s+(.*)\s+(.*)\s+(.*)\s+(.*)\s+(.*)\s+(.*)\s+(.*)\s+(.*)\s+<\/residue>/)                
                {
                    my $resnam = $1;
                    my $resnum = $2;
                    my $molnum = $3;
                    my $a1 = $4;
                    my $a2 = $5;
                    my $a3 = $6;
                    my $a4 = $7;
                    my $second_str = $8;

                    $resnum =~ s/\.//g;

                    #removes insert labels and causes lot of trouble
                    #$resnum =~ s/[a-z]//g;

                    #to make insert symbol upper case again
                    $resnum = uc($resnum);
                    $resnam = uc($resnam);
     
                    my $key = uc($xmas_chain).":".$resnum;
                    my $value = $resnam.":".$a3;
                    $aatypes{$key} = $value; #abs ASA in monomer
                }
            }
        }              
    }        

    close(XMAS);
    return %aatypes;
}


##############################################################
#
sub read_patches
{
    my ($in, $patch_dir, $r_residues, $r_ln, $r_surf_ASAmean) = @_;
    my $infile = "$patch_dir/$in";

    #Kyte & Doolittle   J.Mol.Biol. 157(1982)105-132    
    my %Hphob_values = ('ala' => '1.80',
                        'cys' => '2.50',
                        'asp' => '-3.50',
                        'glu' => '-3.50',
                        'phe' => '2.80',
                        'gly' => '-0.40',
                        'his' => '-3.20',
                        'ile' => '4.50',
                        'lys' => '-3.90',
                        'leu' => '3.80',
                        'met' => '1.90',
                        'asn' => '-3.50',
                        'pro' => '-1.60',
                        'gln' => '-3.50',
                        'arg' => '-4.50',
                        'ser' => '-0.80',
                        'thr' => '-0.70',
                        'val' => '4.20',
                        'trp' => '-0.90',
                        'tyr' => '-1.30');

    open (IN, "$infile") || die "sub read_patches cannot open $infile!!!\n";
    my %patches = ();

    my $linecount = 0;
    
    #runs for every patch, assigns average hydrophobicity value
    while (my $line = <IN>)
    {
        chomp $line;
        
        if ($line =~ /<patch\s(\S+)>\s(.*)$/)
        {
            #central in A.100 format, list is A:73  A:74  A:75  A:76...
            my $central = $1;
            my $list = $2; 
            my $patch_hydrophobicity;
            my $hydrophobicity_sum = 0;

            $central =~ s/\./:/g;            
            $list =~ s/\s+$//;
            my @patch_res = split(' ', $list);
            my $patch_size = scalar @patch_res;

            foreach my $r (@patch_res)
            {
                if (my $val = $$r_residues{$r})
                {
                    my ($aa, $asa) = split(':', $val);
                    $aa = "\L$aa";

                    #eliminates rare cases of UNK as AA type
                    if ($Hphob_values{$aa})
                    {
                        my $res_hydrophobicity = $Hphob_values{$aa};
                        $hydrophobicity_sum += $res_hydrophobicity;
                    }
                    else
                    {
                        print LOG "Residue $r is type $aa, not a standard AA!!!\n";
                        $patch_size--;
                    }
                }
                else
                {
                    print LOG "Residue $r (from patch res list) was not found in XMAS file.\n";
                }
            }
            
            if (!$patch_size)
            {
                print LOG "patch $central has no residues!!!\n";
            }
            else
            {
                $patch_hydrophobicity = $hydrophobicity_sum / $patch_size;
                $patches{$central} = $patch_hydrophobicity;
            }
        }
    }    
    close(IN);
       
    return %patches;
}


##############################################################
sub Usage
{
    print"\nUSAGE:\n";
    print"\t/acrm/usr/local/bin/perl calc_patch_hydrophobicity.pl -patch_dir <dir> -xmas_dir <dir> -log_dir <dir> -out_dir <outfile>\n\n"; 
    print"i.e.   /acrm/usr/local/bin/perl calc_patch_hydrophobicity.pl -patch_dir /home/bsm/anya/intf_prediction/datasets_new/my_dataset/patches_11/ -log_dir log_files/ -out_dir hydrophobicity/\n\n";
    print "This code will produce a file for every patch-containing pqs_idch\n"; 
    print "for every patch in that pqs_id it will output a line:\n";
    print "\t-<patch A.105> X\n";
    print "\twhere A.105 is the central residue for that patch\n";
    print "\tX is average hydrophobicity value, based on Kyle&Doolittle(1982)\n";
    print "\n";
}


##############################################################
#Returns 1 if file is empty, 0 if it is not empty
sub file_empty
{
    my $file = $_[0];
    open(F, "$file") || die "outfile_not_empty cannot open $file!!!\n";

    my @array = <F>;
    my $size = scalar @array;
    close(F);

    if ($size==0)
    {
        return 1;
    }
    return 0;
}
