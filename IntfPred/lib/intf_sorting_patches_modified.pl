#!/acrm/usr/local/bin/perl -w

# TCN THIS SCRIPT IS MODIFIED FROM VERSION AT ~/anya/inft_prediction/code_new. ALL MODIFICATIONS BY ME (TCN) ARE COMMENTED WITH A STARTING TCN

use strict;
use warnings;

use Getopt::Long;  
use lib '/home/bsm/anya/perllib'; #sets dir where Perl_paths.pm is
use Perl_paths; 

#### Define command line inputs ####
my $xmas_dir; # added by TCN
my $patch_dir;
my $sort_dir;
my $rasa_dir; #added on 20.10.2011.
my $absasa_dir; #added on 20.10.2011.
my $log_dir;                        
my $intf;
my $help = '';
my $v_flag = '';


#### Checking command line options ####
#in is dir with <pqs_id>.patches files, out is outfile 
GetOptions(
    'xmas_dir=s'  => \$xmas_dir,
    'patch_dir=s' => \$patch_dir, 
    'intf=s' => \$intf,
    'log_dir=s'   => \$log_dir,
    'sort_dir=s'  => \$sort_dir,
    'rasa_dir=s'  => \$rasa_dir,
    'absasa_dir=s'=> \$absasa_dir,
    'help'        => \$help,
    'v_flag'      => \$v_flag);

if ( (!$patch_dir) || (!$log_dir) || (!$intf) || 
         (!$sort_dir) || (!$rasa_dir) || (!$absasa_dir)
             || (!$xmas_dir) )
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
#14.12.2009. files are in /home/bsm/anya/intf_prediction/datasets/my_dataset/patches_with_inserts/
# for every patch in that pqs_id it will output a line:
#    -<patch A.105> X Y Z
#     where A.105 is the central residue for that patch
#     X is 0/1, if central is/is not one of the intf residues in that pqs_id
#     Y (in %) shows how much of that patch ASA(abs) is in intf
#     Z (in %) shows how much of the intf ASA(abs) is in that patch

#18.12.2009.added INTF output file, where all the intf sizes will be stored, added print of average intf size at the end

#02.02.2010. patch files are in /home/bsm/anya/intf_prediction/datasets_new/my_dataset/patches_11/

#29.06.2010. patch files are in /acrm/home/anya/interface_prediction/patches/patches_11/

#20.10.2011. added print of patch-averaged ASA values in separate files
###################################

####   variables and handles   ####
###################################
#my $intf_residues_file = "/home/bsm/anya/intf_prediction/datasets_new/my_dataset/nonred.intf.residues";

my $Rfile = "patches_intf_overlap_1ratios_singleres.Rdata";
open (R, ">$Rfile") || die "intf_sorting_patches.pl cannot open R!!!\n";

open (INT1, ">intf_singleres.stats") || die "cannot open INT1!\n";

my @intf_totals;

#### main ####
##############
my @filelist = `ls $patch_dir`;

#my $file = "2vqeP.patches";
foreach my $file (@filelist)
{    
    chomp $file;
    $file =~ /(.*)\.patches/;
    my $pqs_id = $1;
    my $chain = chop($pqs_id);
    
    print "$file\n" if ($v_flag);

    # TCN: edited file names to avoid merging dir name with file name
    my $logfile = "$log_dir/".$pqs_id.$chain.".log";
    open (LOG, ">$logfile") || die "intf_sorting_patches.pl cannot open LOG!!!\n";
    my $outfile = "$sort_dir/".$pqs_id.$chain.".sorted.patches"; 
    my $absASAfile = "$absasa_dir/".$pqs_id.$chain.".absASA";
    my $rASAfile = "$rasa_dir/".$pqs_id.$chain.".rASA";
  
    open (OUT, ">$outfile") || die "intf_sorting_patches.pl cannot open OUT!!!\n";
    open (absASA, ">$absASAfile") || die "intf_sorting_patches.pl cannot open absASA!!!\n";
    open (rASA, ">$rASAfile") || die "intf_sorting_patches.pl cannot open rASA!!!\n";

    #key is $chain:$resnum, value is aa_type:absASA
    #[residues] B:443A->ARG:92.207
    my %residues = &res_types_from_xmas( $pqs_id . $chain, $xmas_dir );

    printf "found %d residues in xmas file, ", scalar keys %residues if ($v_flag);

    #key is $chain:$resnum, value is absASA
    #[intf_residue] A:354->79.175
    my %intf_residues = &find_intf_residues($pqs_id, $chain, $intf, \%residues);
    my $intf_with_asa = scalar keys %intf_residues;

    printf "intf has %d residues, ", $intf_with_asa if ($v_flag);

    #here we assume all intf residues belong to only one intf surface 
    #TEST CASES WHEN THIS IS NOT TRUE!!!!
    my $intf_area = &total_intf_area(\%intf_residues);
    
    if(!$intf_area)
    {
        print LOG "total intf area zero!!!!\n";
    }
    else
    {
        push(@intf_totals, $intf_area);

        printf "with intf area of %.2f.\nPatches are:\n", $intf_area if ($v_flag);

        printf INT1 "%d\t%.2f\n", $intf_with_asa, $intf_area;

        #key is $chain:$resnum, value is is_central, %patch in intf, %intf in patch, total absASA patch area, total rASA patch area, number of residues in patch
        my %patches = &read_patches($file, $patch_dir, \%residues, \%intf_residues, $intf_area);

        foreach my $k (sort keys %patches)
        {
            #$R_ASA and $A_ASA are solvent-accessible areas of that patch
            my ($X, $Y, $Z, $A_ASA, $R_ASA, $res_in_patch) = split (':', $patches{$k});
            printf OUT "<patch %s> %d %.2f %.2f\n", $k, $X, $Y, $Z;

            #added on 20.10.2011.
            my $p_absASA = $A_ASA / $res_in_patch;
            printf absASA "<patch %s> %.2f\n", $k, $p_absASA;

            my $p_rASA = $R_ASA / $res_in_patch;
            printf rASA "<patch %s> %.2f\n", $k, $p_rASA;

            printf "<patch %s> %d %.2f %.2f %.2f %.2f %.2f\n", $k, $X, $Y, $Z, $A_ASA, $R_ASA, $res_in_patch if ($v_flag);

            #if the ratios are not (0,0), will save them to be plotted in R
            #the graph is in ~/R/patches/patch_intf_overlap_ratios.ps
            if ( ($Y) || ($Z) ) 
            {
                printf R "$pqs_id:$k\t%.2f\t%.2f\n", $Y, $Z;
            }
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

my ($mean_intf, $stdev_intf, $count) = &calculate_average_intf(@intf_totals);

printf "An average interface (based on %d intfces), has a mean area of %.2f with st.dev. of %.2f.\n", $count, $mean_intf, $stdev_intf;
close(R);
close(INT1);
exit;


##############################################################
# SUBROUTINES
##############################################################
#reads patches in a hash, 
#adds res.type and ASA value in monomer to each hash residue
sub read_patches
{
    my ($in, $patch_dir, $r_residues, $r_intf_res, $intf_area) = @_;
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

            # criterion 1: true if patch central residue is intf residue
            my $is_central = 0;
            #crit.2: percentage of patch area belonging to the interface
            my $perc_patch_in_intf = 'NA';
            #crit.3: percentage of interface area belonging to that patch
            my $perc_intf_in_patch = 'NA';
    
            my $absolute_patch_area = 0;
            my $relative_patch_area = 0;
            my $processed_residues_in_patch = 0;
            my $overlap_area = 0; #interface and patch residues
            $central =~ s/\./:/g;

            if ($$r_intf_res{$central})
            {
                $is_central = 1;
            }

            $list =~ s/\s+$//;
            my @patch_res = split(' ', $list);

            foreach my $r (@patch_res)
            {
                if (my $val = $$r_residues{$r})
                {
                    my ($aa, $absASA, $rASA) = split(':', $val);

                    #added on 20.10.2011.
                    $absolute_patch_area += $absASA;
                    $relative_patch_area += $rASA;
                    $processed_residues_in_patch++;

                    if ($$r_intf_res{$r})
                    {
                        $overlap_area += $absASA;
                    }
                }
                else
                {
                    print LOG "Residue $r (from patch res list) was not found in XMAS file.\n";
                }
            }
            
            if (!$absolute_patch_area)
            {
                print LOG "patch $central has area zero!!!\n";
            }
            else
            {
                $perc_patch_in_intf = ($overlap_area / $absolute_patch_area)*100;
                $perc_intf_in_patch = ($overlap_area / $intf_area)*100;
            }

            #modified on 20.10.2011.
            $patches{$central} = "$is_central:$perc_patch_in_intf:$perc_intf_in_patch:$absolute_patch_area:$relative_patch_area:$processed_residues_in_patch";
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
    my ($pqs_id, $xmas_dir) = @_;
    my %aatypes = ();

    my $inatoms_datatype = 0;
    my $inmolecules = 0;
    my $type = 0;
    my $xmas_chain;

    print "TEST: $pqs_id\n";
    
    # TCN: prcoess 'pqs_id' to match format of xmas file names
    my $pdb_code = uc substr($pqs_id, 0, 4);
    my $chain_id = uc substr($pqs_id, 4, 1);
    
    # TCN: xmas file is now specified using xmas dir
    my $xmas_file
        = "$xmas_dir/" . $pdb_code . '_' . $chain_id . '.xmas';
    
    open (XMAS, "$xmas_file") || die "sub res_types_from_xmas cannot open $xmas_file, $!";
        
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
                    my $a3 = $6;  # abs ASA in monomer
                    my $a4 = $7;  # rel ASA in monomer
                    my $sec_str = $8;
                    
                    $resnum =~ s/\.//g;

                    #removes insert labels and causes lot of trouble
                    #$resnum =~ s/[a-z]//g;

                    #to make insert symbol upper case again
                    $resnum = uc($resnum);
                    $resnam = uc($resnam);
     
                    my $key = uc($xmas_chain).":".$resnum;
                    my $value = $resnam.":".$a3.":".$a4;
                    $aatypes{$key} = $value; #abs and rel ASA in monomer
                }
            }
        }              
    }        

    close(XMAS);
    return %aatypes;
}

##############################################################
#
sub find_intf_residues
{
    my ($pqs_id, $pqs_ch, $intffile, $ref_reshash) = @_;
    my %intf_res = ();

    open (INTF, "$intffile") || die "[find_intf_residues] cannot open $intffile!!!\n";
    while (my $line = <INTF>)
    {
        chomp $line;
        if ($line =~ /$pqs_id/i)
        {
            $line =~ /\->\((.*)\)/;
            my $a = $1;
            chop $a; #removing last ,
            my @array = split(',', $a);

            foreach my $res (@array)
            {
                $res = uc($res);
                my ($ch, $num, $aa) = split(':', $res);

                if ($ch eq $pqs_ch)
                {
                    if (my $value = $$ref_reshash{"$ch:$num"})
                    {
                        my ($xmas_aa_type, $absASA, $rASA) = split(':', $value);
                        if ($xmas_aa_type ne $aa) 
                        {
                            print LOG "for residue $ch:$num aa_type in intf file is $aa and in xmas file is $xmas_aa_type!!!\n";
                        }
                        else        
                        {
                            $intf_res{"$ch:$num"} = "$absASA:$rASA";
                            #print "intf: $res $absASA\n";
                        }
                    }
                    else
                    {
                        print LOG "No record for intf residue $res.\n";
                    }
                }
            }
        }
    }
    close(INTF);

    return %intf_res;
}

##############################################################
#calculates total intf area, based on absolute ASA of residues in a hash
sub total_intf_area
{
    my $hash_ref = $_[0];
    my $total_area = 0;

    while (my ($k, $v) = each %$hash_ref)
    {
        my ($absASA, $rASA) = split(':', $v);
        $total_area += $absASA;
    }

    return $total_area;
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
sub calculate_average_intf
{
    my @intf_totals_array = @_;
    my $total = 0;
    my $count = 0;
    my $stdev_sum = 0;

    foreach my $a (@intf_totals_array)
    {
        $total += $a;
        $count++;
    }
    my $mean = $total / $count;
    printf "added %d intf areas, sum is %.2f, mean is %.2f.\n", $count, $total, $mean;
    
    foreach my $a (@intf_totals_array)
    {
        my $square = ($a-$mean)*($a-$mean);
        $stdev_sum += $square;
        #print "square is $square, stdev sum is $stdev_sum.\n" if ($v_flag);
    }

    my $stdev = sqrt( $stdev_sum / $count);

    return ($mean, $stdev, $count);
}

##############################################################
sub Usage
{
    print"\nUSAGE:\n";
    print"\t/acrm/usr/local/bin/perl intf_sorting_patches.pl -xmas_dir <dir> -patch_dir <dir> -intf <intf_file> -log_dir <dir> -sort_dir <dir> -rasa_dir <dir> -absasa_dir <dir>\n\n"; 
    print"i.e.   /acrm/usr/local/bin/perl intf_sorting_patches.pl -xmas_dir /acrm/data/xmas/pdb/ -patch_dir /home/bsm/anya/intf_prediction/datasets_new/my_dataset/patches_11/ -intf /acrm/home/anya/interface_prediction/patches/nonred2.intf.residues -log_dir log_files/ -sort_dir intf_sorting_patches/ -absasa_dir absASA/ -rasa_dir rASA/\n\n";
    print "This code will produce three files for every patch-containing pqs_idch (class file, rASA and absASA)\n"; 
    print "for every patch in that pqs_id in the sorting file it will output a line:\n";
    print "\t-<patch A.105> X Y Z\n";
    print "\twhere A.105 is the central residue for that patch\n";
    print "\tX is 0/1, if central is/is not one of the intf residues in that pqs_id\n";
    print "Y (in %) shows how much of that patch ASA(abs) is in intf\n";
    print "Z (in %) shows how much of the intf ASA(abs) is in that patch\n";
    print "rASA file has total rASA patch area, divided by the number of res in that patch\n";
    print "absASA file has total absASA patch area, divided by the number of res in that patch\n";
    print "TCN: added -xmas_dir\n";
}
