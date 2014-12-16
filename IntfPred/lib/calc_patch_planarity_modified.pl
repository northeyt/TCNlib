#!/acrm/usr/local/bin/perl -w

use strict;
use DBI;
use Getopt::Long;  
use lib '/home/bsm/anya/perllib'; #sets dir where Perl_paths.pm is
use Perl_paths;
use Carp;
use IO::CaptureOutput qw( qxx );

#### Define command line inputs ####
my $patch_dir;
my $pdb_dir;
my $out_dir;
my $log_dir;         
my $help = '';
my $v_flag = '';


#### Checking command line options ####
#in is dir with <pqs_id>.patches files,
GetOptions('patch_dir=s' => \$patch_dir,
           'pdb_dir=s' => \$pdb_dir,
           'log_dir=s' => \$log_dir,
           'out_dir=s' => \$out_dir, 
           'help' => \$help,
           'v_flag' => \$v_flag);
#i.e. /acrm/usr/local/bin/perl calc_patch_planarity.pl -patch_dir /acrm/home/anya/interface_prediction/patches/patches_11/ -log_dir log_files/ -out_dir planarity/ -v_flag

if ( (!$patch_dir) || (!$log_dir) || (!$out_dir)
         || (!$pdb_dir) )
{
    $help = 1;
}

if ($help)
{
    &Usage();
    exit;
}

####   variables and handles   ####
###################################

#### main ####
##############
my @p_filelist = `ls $patch_dir`;

#my $p_file = "1is3A.patches";
foreach my $p_file (@p_filelist)
{
    chomp $p_file;
    print "$p_file..."; 
    my $count = 0;                         

    $p_file =~ /(.*)\.patches/;
    my $pqs_id = $1;
    my $chain = chop($pqs_id);

    my $logfile = "$log_dir/$pqs_id.planarity.log";
    open(LOG, ">$logfile") || die "calc_patch_planarity: $logfile not opened\n";
    my $patchfile_path = "$patch_dir/$p_file";
    open(PFILE, "$patchfile_path") || die "calc_patch_planarity (PFILE handle): $p_file not opened\n";

    my $outfile = "$out_dir/$pqs_id.$chain.patch.rms";
    open(OUT, ">$outfile") || die "calc_patch_planarity (OUT handle): $outfile not opened\n";

    while (my $line = <PFILE>)
    {
        chomp $line;
        $count++;
        print "line is $line.\n" if ($v_flag);

        if (length($line))
        {
            $line =~ /patch\s(\S+)>\s(.*)$/;
            my $central_res = $1;
            my $patch_r = $2;
            my @patch_residues = split(' ', $patch_r);

            my $rms = &calc_patch_rms(\@patch_residues, $pqs_id, $chain, $central_res);    
            
            #write output
            print OUT "<patch $central_res> $rms\n";

            if ($rms eq 'N/A')
            {
                print LOG "<patch $central_res> $rms\n";
            }
        }  
        else
        {
            print LOG "Empty line in patchfile!\n";
        }
    }    
    close(OUT);
    close(LOG);
    close(PFILE); 

print " $count patches. Done!\n";

    #cleanup for that patch file...
    if(&file_empty($logfile))
    {                           
        print "rm -f $logfile\n\n";
        `rm -f $logfile`;
    }
    else
    {
        print "$logfile had errors!!!\n\n";
    } 
}
  
exit;

##############################################################
# SUBROUTINES
##############################################################
sub calc_patch_rms
{
    my ($res_array_ref, $pqs_id, $ch, $central) = @_;
    my $res_count = scalar @$res_array_ref;
    my $rms = "N/A";

    if (!$res_count)
    {
        print LOG " $pqs_id$ch $central patch had no residues!\n";
    }
    else
    {
        my $reslist_file = &patch_res2reslist($res_array_ref, $pqs_id, $ch, $central);
        my $pdb_file = &reslist2pdb($reslist_file, $pqs_id, $ch, $central);
        $rms = &pdbfile2rms($pdb_file);
    }
    return $rms;
}           


##############################################################
#inputs a file of residues, 
#extracts ATOM entries for those residues from the PQS file
sub reslist2pdb
{
    my ($file, $id, $ch, $central) = @_;
    
    #create PDB file with only reslist
    my $getres_exe = $Perl_paths::martin_bin."getresidues";
    # TCN: need to change this file to pdb dir
    my $pqs_file = "$pdb_dir/" . uc("$id" . "_$ch") . ".pdb";

    my $pdbformat =  $Perl_paths::tmp_dir.$id.$ch."_".$central.".reslist.pqs";

#    print "$getres_exe $file $pqs_file $pdbformat\n";
    `$getres_exe $file $pqs_file $pdbformat`;

    print "rm -f $file\n" if ($v_flag);
    `rm -f $file`;

    return $pdbformat;
}


##############################################################
#prints all residues to a file
sub patch_res2reslist
{
    my ($resarray_ref, $pqs_id, $ch, $central) = @_;

    #write all patch residues in a reslist file
    my $res_list_file = $Perl_paths::tmp_dir.$pqs_id.$ch."_".$central.".reslist";
    open(L, ">$res_list_file") || die "calc_patch_planarity (R handle): $res_list_file not opened\n";

    foreach my $r (@$resarray_ref)
    {
        $r =~ s/:/\./g;
        print L "$r\n";
    }
    close(L);

    printf "for %s%s %s printed %d residues.\n", $pqs_id, $ch, $central, scalar @$resarray_ref if ($v_flag);

    return $res_list_file;
}


##############################################################
#inputs a PDB file with all residues in a patch
#runs PRINCIP (from SURFNET package by R. Laskowski) and extracts RMS
sub pdbfile2rms
{
    my $file = $_[0];
    my $rms = "N/A";
    
    #output PRINCIP results
    # TCN: 64-bit princip now used
    my $planarity_exe = '/acrm/usr/local/apps/surfnet/bin64/princip';

    my ($result, $error, $success) = qxx("echo $file | $planarity_exe");  

    if ($error) {
        # croak if this script is not running on a 32-bit machine
        if ($error =~ /error while loading shared libraries/) {
            croak $error;
        }
        else {
            print $error;
        }
    }
    #parse output and extract RMS from best-fit plane
    my @array = split(/\n/, $result);
    foreach my $l (@array)
    {
        chomp $l;

        my $string = "RMS difference from best\-fit plane:";
        if ($l =~ /$string\s+(.*)$/)
        {
            $rms = $1;   
        }
    }  

    print "rm -f $file\n" if ($v_flag);
    `rm -f $file`;

    return $rms;
}


##############################################################
sub Usage
{
    print"\nUSAGE:\n";
    print"\t/acrm/usr/local/bin/perl calc_patch_planarity.pl -patch_dir <dir> -log_dir <dir> -out_dir <dir> [-v_flag]\n\n"; 
    print"i.e.   /acrm/usr/local/bin/perl calc_patch_planarity.pl -patch_dir /acrm/home/anya/interface_prediction/patches/patches_11/ -log_dir log_files/ -out_dir planarity/ -v_flag \n\n";
    print "This code will produce a file for every patch-containing pqs_idch\n"; 
    print "for every patch in that pqs_id it will output a line:\n";
    print "\t-<patch A.105> X\n";
    print "\twhere A.105 is the central residue for that patch\n";
    print "\tX is the RMS from the best-fit plane (R. Laskowski's PRINCIP in SURFNET package).\n";
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

##############################################################
