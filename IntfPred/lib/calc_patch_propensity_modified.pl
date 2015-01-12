#!/acrm/usr/local/bin/perl -w

use strict;
use Getopt::Long;  
use lib '/home/bsm/anya/perllib'; #sets dir where Perl_paths.pm is
use Perl_paths; 
use Carp;

#### Define command line inputs ####
my $patch_dir;
my $xmas_dir;
my $out_dir;
my $log_dir;                        
my $intf_summary;
my $surf_summary;
my $help = '';
my $v_flag = '';


#### Checking command line options ####
#in is dir with <pqs_id>.patches files, out is outfile 
GetOptions('patch_dir=s' => \$patch_dir,
           'xmas_dir=s', => \$xmas_dir,
           'intf_summary=s' => \$intf_summary,
           'surf_summary=s' => \$surf_summary,
           'log_dir=s' => \$log_dir,
           'out_dir=s' => \$out_dir, 
           'help' => \$help,
           'v_flag' => \$v_flag);
#i.e. /acrm/usr/local/bin/perl calc_patch_propensity.pl -intf_summary intf.prepare.out -surf_summary surf.prepare.out -patch_dir ~/intf_prediction/datasets_new/my_dataset/patches_11/ -log_dir log_files/ -out_dir propensities/ -v_flag

if ( (!$patch_dir) || (!$log_dir) || 
     (!$out_dir) || (!$intf_summary) || 
     (!$surf_summary) || (!$xmas_dir) )
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
#23.02.2010. files are in /home/bsm/anya/intf_prediction/datasets_new/my_dataset/patches_11/
# for every patch in that pqs_id it will output a line:
#    -<patch A.105> X
#     where A.105 is the central residue for that patch
#     X is an average propensity of that patch, calculated as sum of propensities of residues P(res) in that patch, divided by the #residues in the patch
#     [P(res) = (ln{P_res_intf / P_res_surf}) * (S_res / S_res_ave)]
#     - P_res is % of rASA in intf/surf dataset belonging to that res type
#     - S_res is actual rASA extracted from xmas file for that residue
#     - S_res_aver is average rASA for res type in SURFACE dataset

####   variables and handles   ####
###################################
my $Rfile = "patches_propensities_r11.Rdata";
open (R, ">$Rfile") || die "calc_patch_propensity.pl cannot open R!!!\n";

my %AA_surf_ASAmean = ('ala' => '0',
                       'cys' => '0',
                       'asp' => '0',
                       'glu' => '0',
                       'phe' => '0',
                       'gly' => '0',
                       'his' => '0',
                       'ile' => '0',
                       'lys' => '0',
                       'leu' => '0',
                       'met' => '0',
                       'asn' => '0',
                       'pro' => '0',
                       'gln' => '0',
                       'arg' => '0',
                       'ser' => '0',
                       'thr' => '0',
                       'val' => '0',
                       'trp' => '0',
                       'tyr' => '0');
#### main ####
##############
my @filelist = `ls $patch_dir`;

 my %intf_averages = &extract_percentages($intf_summary);
my %surf_averages = &extract_percentages_and_averages($surf_summary, \%AA_surf_ASAmean);
 
my %ln = &calc_ln_part(\%intf_averages, \%surf_averages);
                                
foreach my $k (keys %ln)
{
    print "[ln]$k\t$ln{$k}\n";
}

#my $file = "1bdgA.patches";
foreach my $file (@filelist)
{
    chomp $file;
    $file =~ /(.*)\.patches/;
    my $pqs_id = $1;
    my $chain = chop($pqs_id);

    print "$file\n" if ($v_flag);

    my $logfile = "$log_dir/".$pqs_id.$chain.".propensity.log";
    open (LOG, ">$logfile") || die "intf_sorting_patches.pl cannot open LOG!!!\n";
    my $outfile = "$out_dir/".$pqs_id.$chain.".propensity.patches"; 
    open (OUT, ">$outfile") || die "intf_sorting_patches.pl cannot open OUT!!!\n";

    #key is $chain:$resnum, value is aa_type:absASA
    #[residues] B:443A->ARG:92.207
    my %residues = eval {&res_types_from_xmas($pqs_id, $chain, $xmas_dir)};
    if (! %residues) {
        if ($@ =~ /cannot open/) {
            print "No XMAS file found for $pqs_id, skipping ...\n";
            next;
        }
        else {
            croak $@;
        }
    }
    
    my %patches = &read_patches($file, $patch_dir, \%residues, \%ln, \%AA_surf_ASAmean);

    foreach my $k (sort keys %patches)
    {
        #$A is area of that patch
        printf OUT "<patch %s> %.4f\n", $k, $patches{$k};
               
        #the graph is in ~/R/patches/patches_propensities_r11.Rdata
        if ($patches{$k})
        {
            printf R "$pqs_id:$k\t%.2f\n", $patches{$k};
        }
    }    
    close(OUT);
    close(LOG);
    
    if(&file_empty($logfile))
    {
        print "Removing $logfile ...\n";
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
#reads patches in a hash, 
#adds res.type and ASA value in monomer to each hash residue
sub read_patches
{
    my ($in, $patch_dir, $r_residues, $r_ln, $r_surf_ASAmean) = @_;
    my $infile = "$patch_dir/$in";

    open (IN, "$infile") || die "sub read_patches cannot open $infile!!!\n";
    my %patches = ();

    my $linecount = 0;
    
    #runs for every patch, assigns three values:
    #is central res of a patch an intf res
    #% of patch area in intf
    #%of intf area in that patch
    while (my $line = <IN>)
    {
        chomp $line;
        
        if ($line =~ /<patch\s(\S+)>\s(.*)$/)
        {
            #central in A.100 format, list is A:73  A:74  A:75  A:76...
            my $central = $1;
            my $list = $2; 
            my $patch_propensity;
            my $propensity_sum = 0;

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
                    if ($$r_surf_ASAmean{$aa})
                    {
                        my $res_propensity = $$r_ln{$aa} * ($asa / $$r_surf_ASAmean{$aa});
                        $propensity_sum += $res_propensity;
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
            
            if ( (!$propensity_sum) || (!$patch_size) )
            {
                print LOG "patch $central has no propensity value or no residues!!!\n";
            }
            else
            {
                $patch_propensity = $propensity_sum / $patch_size;
                $patches{$central} = $patch_propensity;
            }
        }
    }    
    close(IN);
       
    return %patches;
}

##############################################################
#read in a pqsatoms XMAS file
#outputs hash of residues (key is uc(chain):resnum), value is res_type:absASA
sub res_types_from_xmas
{
    my ( $pqs_id, $chain, $xmas_dir ) = @_;
    my %aatypes = ();

    my $inatoms_datatype = 0;
    my $inmolecules = 0;
    my $type = 0;
    my $xmas_chain;

#    my $xmas_file = "/grid/gridstore/xmas/pqsatoms/pqs".$pqs_id.".xmas";;
    my $xmas_file
        = "$xmas_dir/" . uc ($pqs_id) . uc ("_$chain") . '.xmas';
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
                    my $sec_str = $8;
                    
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


##############################################################
sub Usage
{
    print"\nUSAGE:\n";
    print"\t/acrm/usr/local/bin/perl calc_patch_propensity.pl -patch_dir DIR -xmas_dir DIR -intf_summary FILE -surf_summary FILE -log_dir DIR -out_dir DIR\n\n"; 

    print"i.e.   /acrm/usr/local/bin/perl calc_patch_propensity.pl -intf_summary intf.prepare.out -surf_summary surf.prepare.out -patch_dir /acrm/home/anya/interface_prediction/patches/patches_11/ -log_dir log_files/ -out_dir /acrm/home/anya/interface_prediction/propensities/\n\n";

    print "This code will produce a file for every patch-containing pqs_idch\n"; 
    print "for every patch in that pqs_id it will output a line:\n";
    print "\t-<patch A.105> X\n";
    print "\twhere A.105 is the central residue for that patch\n";
    print "\tX is the averaged patch propensity\n\n";
}


##############################################################
sub extract_percentages
{
    my $file = $_[0];
    my %averages = ();
    my $total = 0;

    open(IN, "$file") || die "sub extract_averages: Cannot open IN!!!\n";

    while (my $line = <IN>)
    {
        chomp $line;
        if ($line =~ /\[type:sum:count:mean\](.*)/)
        {
            my ($aa_type, $ASAsum, $aa_count, $ASAmean) = split(/:/, $1);
            
            $ASAsum =~ s/\s+//g;
            $averages{$aa_type} = $ASAsum;
        }
        
        if ($line =~ /total ASA for the dataset is\s(\S+)\./)
        {
            $total = $1;
        }
    }
    close(IN);

    foreach my $k (keys %averages)
    {
        $averages{$k} /= $total;
    }

    return %averages;
}


##############################################################
sub extract_percentages_and_averages
{
    my ($file, $hash_ref) = @_;
    my %averages = ();
    my $total = 0;

    open(IN, "$file") || die "sub extract_averages: Cannot open IN!!!\n";

    while (my $line = <IN>)
    {
        chomp $line;
        if ($line =~ /\[type:sum:count:mean\](.*)/)
        {
            my ($aa_type, $ASAsum, $aa_count, $ASAmean) = split(/:/, $1);
            
            $ASAsum =~ s/\s+//g;
            $averages{$aa_type} = $ASAsum;

            $ASAmean =~ s/\s+//g;
            $$hash_ref{$aa_type} = $ASAmean;
        }
        
        if ($line =~ /total ASA for the dataset is\s(\S+)\./)
        {
            $total = $1;
        }
    }
    close(IN);

    foreach my $k (keys %averages)
    {
        $averages{$k} /= $total;
    }

    return %averages;
}

##########################################################################
#calculates log_base2_($n), default log function uses base e 
sub log2
{
    my $n = $_[0];
    my $log2 = log($n) / log(2);
    return $log2;
}

##########################################################################
sub calc_ln_part
{
    my ($r_intf, $r_surf) = @_;
    my %ln = ();

    foreach my $k (sort keys %$r_intf)
    {
        my $a = $$r_intf{$k} / $$r_surf{$k};
        my $ln_prop = log($a);
        $ln{$k} = $ln_prop;
    }

    return %ln;
}
