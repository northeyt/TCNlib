#!/acrm/usr/local/bin/perl -w
use strict;
use Getopt::Long;  
use lib '/home/bsm/anya/perllib'; #sets dir where Perl_paths.pm is
use Perl_paths; 

##### Description #####
#This code extracts S-S bonds, H-bonda and secondary structure info from the XMAS files. S-S bonds are done distance-based only, whereas calc_patch_SSgeom.pl hasan angle restriction as well and produces XXXX ?????? less S-S bond records.
#######################

##### Variables, handles and command line options #####
my $patch_dir;
my $xmas_dir;
my $SS_dir;
my $Hb_dir;
my $secstr_dir;
my $verbose;
my $help;

my $SS_flag = 0;
my $Hb_flag = 0;
#my $secstr_dir = "/acrm/home/anya/interface_prediction/properties_singleres/second_str/";

GetOptions( 'patch_dir=s' => \$patch_dir,
            'xmas_dir=s'  => \$xmas_dir,
            'SS=s'        => \$SS_dir,
            'Hb=s'        => \$Hb_dir,
            'secstr=s'    => \$secstr_dir,
            'verbose'     => \$verbose,
            'help'        => \$help);

if ( (!$patch_dir) || (!$xmas_dir) || (!$secstr_dir) )
{
    $help = 1;
}

if ($help)
{
    &Usage();
    exit;
}

if($SS_dir)
{
    $SS_flag = 1;
}

if($Hb_dir)
{
    $Hb_flag = 1;
}

##### Main #####
my @p_filelist = `ls $patch_dir`;

#my $p_file = "1aocA.patches";
foreach my $p_file (@p_filelist)
{
    chomp $p_file;
    print "$p_file...\n"; 
    my $count = 0;                         

    $p_file =~ /(.*)\.patches/;
    my $pqs_id = $1;
    my $chain = chop($pqs_id);

    my %SS = ();
    my %Hb = ();

    my %sec_str = &extract_SS_Hb_2D_from_XMAS($pqs_id, $chain, $SS_flag, \%SS, $Hb_flag, \%Hb, $xmas_dir);

    printf "Have secondary structures for %d residues in $pqs_id$chain\n", scalar keys %sec_str if ($verbose);

&output_patch_3params($p_file, $patch_dir, $secstr_dir, \%sec_str, $SS_flag, \%SS, $SS_dir, $Hb_flag, \%Hb, $Hb_dir);

}

exit;

##### SUBROUTINES #####
##############################################################
sub output_patch_3params
{
    my ($file, $indir, $secstr_dir, $r_2D_hash, 
        $SSflag, $r_SS_hash, $SS_dir, $Hbflag, $r_Hb_hash, $Hb_dir) = @_;

    my $path = "$indir/$file";
    open(P, "$path") || die "calc_SS_Hb_2Dstr.pl (P handle): $path not opened\n";

    my $idch = $file;
    $idch =~ s/\.patches//;

    my $secstr_path = "$secstr_dir/$idch." . ".2D";
    open(STR, ">$secstr_path") || die "calc_SS_Hb_2Dstr.pl (STR handle): $secstr_path not opened\n";

    if($SSflag)
    {
        my $SS_path = "$SS_dir/$idch".".SSbonds";
        open(SS, ">$SS_path") || die "calc_SS_Hb_2Dstr.pl (SS handle): $SS_path not opened\n";
    }

    if($Hbflag)
    {
        my $Hb_path = "$Hb_dir/$idch." . ".Hbbonds";
        open(Hb, ">$Hb_path") || die "calc_SS_Hb_2Dstr.pl (Hb handle): $Hb_path not opened\n";
    }

    while (my $line = <P>)
    {
        chomp $line;

        if (length($line))
        {
            $line =~ /patch\s(\S+)>\s(.*)$/;
            my $central_res = $1;
            my $patch_r = $2;
            my @patch_residues = split(' ', $patch_r);

            #Secondary structure part
            my $patch_2Dstr = &calc_2D_for_patch(\@patch_residues, $r_2D_hash);
            print STR "<patch $central_res> $patch_2Dstr\n";

            if($SSflag)
            {
                #SSbonds part
                my $patch_SS = &calc_bondcount_for_patch(\@patch_residues, $r_SS_hash);
                printf SS "<patch %s> %.6f\n", $central_res, $patch_SS; 
            }

            if($Hbflag)
            {
                #H-bonds part
                my $patch_Hb = &calc_bondcount_for_patch(\@patch_residues, $r_Hb_hash);
                printf Hb "<patch %s> %.6f\n", $central_res, $patch_Hb;
            }
        }
    }

    close(P);
    close(STR);

    if($SSflag)
    {             
        close(SS);
    }

    if($Hbflag)
    {
        close(Hb);
    }
}

##############################################################
#INPUT: xmas filename (here: from PQS file, contains only protein chains) 
#extracts H-bonds, SS and 2D data for every residue in a given chain
#OUTPUT: hash containing centres of patches, key is residue [chain]resnum[insert], value is 2D:Hb:SS
sub extract_SS_Hb_2D_from_XMAS
{
    my ($pqsid, $pqschain, $SSflag, $r_SShash,
        $Hbflag, $r_Hbhash, $xmas_dir) = @_;
    
    my %sec_str = ();

    my $xmas_file = "$xmas_dir/" . uc($pqsid . "_$pqschain") . '.xmas';

    if (-e $xmas_file)
    {
        open(XMAS, "$xmas_file") || die "calc_SS_Hb_2Dstr.pl (XMAS handle): $xmas_file not opened\n";

        my $inatoms_datatype = 0;
        my $inmolecules = 0;
        my $indata = 0;
        my $type = 0;
        my $chain = 'not defined';
        my $inSS_datatype = 0;     
        my $inHb_datatype = 0; 
        
        my $lc_pqschain = "\L$pqschain";
        
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
                    if ( ($lc_pqschain =~ $chain) && 
                         ($line =~ /^<residue>\s+(.*)\s+(.*)\s+(.*)\s+(.*)\s+(.*)\s+(.*)\s+(.*)\s+(.*)\s+<\/residue>/) )                
                    {
                        my $resnam = $1;
                        my $resnum = $2;
                        my $molnum = $3;
                        my $a1 = $4;
                        my $a2 = $5;
                        my $a3 = $6;
                        my $a4 = $7;
                        my $second_str = $8;
                        
                        $resnum =~ s/\.//g;     #removing dots
                        #$resnum =~ s/[a-z]//g; #watch out for inserts!!!
                        #my $resid = $pqschain.":".$resnum.":".$resnam;
                        my $resid = $pqschain.":".$resnum;
                        
                        $sec_str{$resid} = $second_str;
                    }
                }
            }


            ##### Optional SSbond section #####
            ###################################
            # Clear flags on leaving ssbond data section
            if ( ($inSS_datatype) && ($line =~ /^<\/data/) )
            {	
                $inSS_datatype = 0;
            }
            
            # Handle SSBOND data section
            if($inSS_datatype)
            {
                my ($atomnum1, $ch1, $resnum1, $atomnum2, $ch2, $resnum2) = split(' ', $line);
                
                $resnum1 =~ s/\.//g;
                $resnum2 =~ s/\.//g;
   
                if( ($ch1 eq $lc_pqschain) && ($ch2 eq $lc_pqschain) )
                {
                    $$r_SShash{"$pqschain:$resnum1"} = 1;
                    $$r_SShash{"$pqschain:$resnum2"} = 1;
                }
            } 
            
            #Check for entering data ssbond section
            if ( ($line =~ /<data type=ssbond>/) && ($SSflag) )
            {
                $inSS_datatype = 1;
            }
            
            
            ##### Optional Hbonds section #####
            ###################################
            # Clear flags on leaving ssbond data section
            if ( ($inHb_datatype) && ($line =~ /^<\/data/) )
            {	
                $inHb_datatype = 0;
            }
            
            # Handle PPHBONDS data section
            if($inHb_datatype)
            {

                # Fixed split error where columns where whitespace was
                # being assigned to vars (repalced /\s/ with ' ')
                # (I'm guessing Anja's Hbond data was buggy because of this)
                my ($atomnum1, $atomnum2, $resnam1, $ch1, $resnum1, $atomnam1, $resnam2, $ch2, $resnum2, $atomnam2) = split(' ', $line);
                
                $resnum1 =~ s/\.//g;
                $resnum2 =~ s/\.//g;
                $atomnam1 =~ s/\.//g;
                $atomnam2 =~ s/\.//g;

                if( ($ch1 eq $lc_pqschain) && ($ch2 eq $lc_pqschain) )
                {
                    $$r_Hbhash{"$pqschain:$resnum1"} = 1;
                    $$r_Hbhash{"$pqschain:$resnum2"} = 1;
                }
            }

            #Check for entering data pphbonds section
            if ( ($line =~ /<data type=pphbonds>/) && ($Hbflag) )
            {
                $inHb_datatype = 1;
            }
        }
        close(XMAS);
    
        if($verbose)
        {
            my $SScount = scalar keys %$r_SShash;
            printf "Found %d residues involved in SS bonds.\n", $SScount if ($SSflag);
            
            my $Hbcount = scalar keys %$r_Hbhash;
            printf "Found %d residues involved in H-bonds.\n",  $Hbcount if ($Hbflag);
        }
    }
    else
    {
        print "$xmas_file:$pqsid$pqschain not found\n";
    }

    return %sec_str;
}


##############################################################
#threshold is set to 20%
sub calc_2D_for_patch
{
    my($r_p_res_array, $r_2D_hash) = @_;
    
    my $Ecount = 0;
    my $Hcount = 0;
    my $total = 0;

    my $patch_2D_type = "N/A";

    foreach my $res (@$r_p_res_array)
    {
        my $res_2D_type = $$r_2D_hash{$res};
        
        if($res_2D_type)
        {
            $total++;

            if($res_2D_type eq 'e')
            {
                $Ecount++;
            }
            
            if($res_2D_type eq 'h')
            {
                $Hcount++;
            }
        }
    }

    if($total)
    {
        my $E = $Ecount/$total;
        my $H = $Hcount/$total;

        if( ($E > 0.20) && ($H > 0.20) )
        {
            $patch_2D_type = "EH";
        }
        elsif ($E > 0.20)
        {
            $patch_2D_type = "E";
        }
        elsif ($H > 0.20)
        {
            $patch_2D_type = "H";
        }
        else
        {
            $patch_2D_type = "C";
        }
    }

    return $patch_2D_type;
}


##############################################################
#N/A if no residues in the patch have SS/Hbonds value calculated 
#else averaged SS/Hbonds count for all residues in that patch
sub calc_bondcount_for_patch
{
    my ($r_p_res_array, $r_bonds_hash) = @_;
    my $rescount = scalar @$r_p_res_array;

    my $patch_bonds_count = 0;
    my $patch_bonds;

    foreach my $res (@$r_p_res_array)
    {
        if(exists $$r_bonds_hash{$res})
        {
            $patch_bonds_count++;
        }
    }
    
    if($rescount)
    {
        $patch_bonds = $patch_bonds_count/$rescount;
    }
    else
    {
        $patch_bonds = "N/A";
    }

    return $patch_bonds;
}


##############################################################
sub Usage
{
    print "\nUSAGE:\n";
    print "/acrm/usr/local/bin/perl calc_SS_Hb_2Dstr.pl -patch_dir dir -xmas_dir dir -secstr dir [-help] [-SS] [-Hb] \n\n";
    print "\t-patch_dir /acrm/home/anya/interface_prediction/patches/patches_11/\n";
    print "\t-SS_dir /acrm/home/anya/interface_prediction/SSbonds/\n";
    print "\t-Hb_dir /acrm/home/anya/interface_prediction/Hbonds/\n";
    print "\t-secstr_dir /acrm/home/anya/interface_prediction/second_str/\n\n";
}
