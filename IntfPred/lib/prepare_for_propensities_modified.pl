#!/acrm/usr/local/bin/perl -w

use strict;
use Getopt::Long;   
use lib '/home/bsm/anya/intf_prediction/code/properties/';
use propensities_hardcodes;

########################
#### Input handling ####
########################
my $in;          #parse here: nonred2.intf.residues
my $xmas_dir;
my $out;
my $help = '';   #prints usage
my $v_flag = ''; #for testing

GetOptions('in=s' => \$in, 
           'out=s' => \$out,
           'xmas_dir=s' => \$xmas_dir,
           'help' => \$help,
           'v_flag' => \$v_flag);

if ( (!$in) || (!$out) || (!$xmas_dir) )
{
    $help = 1;
}

if ($help)
{
    &Usage;
    exit;
}

###############################
#### Variables and handles ####
###############################
my @aas = ('ala', 'cys', 'asp', 'glu', 'phe', 'gly', 'his', 'ile', 'lys', 'leu', 'met', 'asn', 'pro', 'gln', 'arg', 'ser', 'thr', 'val', 'trp', 'tyr');

my %AA_sums = ('ala' => '0',
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

my %AA_counts = ('ala' => '0',  
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

open(LOG, ">prepare.log") || die "Cannot open LOG!!!\n";

#################
### Main code ###
#################
#1. get a list of all PQS_idch-s to process
#(patch.residues have the same PQS chains as the cluster.representatives)
open(IDCH, "$in") || die "Cannot open IDCH!!!\n";
my @in = <IDCH>;
close(IDCH);

#2. foreach pqs_idch, extract rASA values for all residues
#key is $chain:$resnum, value is aa_type:absASA
#[residues] B:443A->ARG:92.207

#my $line = "1bcp_1D->";
foreach my $line (@in)
{
    chomp $line;
    $line =~ /^(\S+)\->/;
    my $pqs_id = $1;
    my $chain = chop($pqs_id);
    print "$pqs_id$chain...\n";

    my %xmas_res = &res_types_from_xmas($pqs_id, $chain, $xmas_dir);

printf "[res_types_from_xmas] found %d residues for %s%s in xmas file.\n", scalar keys %xmas_res, $pqs_id, $chain if ($v_flag);

    #3. extract intf residues for that pqs_idch
    #key is $chain:$resnum, value is absASA
    #[intf_residue] A:354->ARG:79.175
    my %residues = &find_intf_residues($pqs_id, $chain, $in, \%xmas_res);

printf "[find_intf_residues] found %d residues for %s%s.\n", scalar keys %residues, $pqs_id, $chain if ($v_flag);

    #4. sum the area up for each AA type 
    #each array element will have the sum of intf area for that AA type
    &ASA_per_AA_type(\%residues, \%AA_sums, \%AA_counts);

    if ($v_flag)
    {
        foreach my $k (sort keys %AA_sums) 
        {
            print "[AA_sums:$pqs_id$chain]$k->$AA_sums{$k}\n";
        }        
   
        foreach my $k (sort keys %AA_counts) 
        {
            print "[AA_counts:$pqs_id$chain]$k->$AA_counts{$k}\n";
        }        
    } 
}              
                 
#5. calculate average rASA for every AA type in surface dataset
my %mean_ASA = &calculate_mean_ASA(\%AA_sums, \%AA_counts);

#6. calculate total area of the dataset (sum all AA types)
my $totalASA = &total_area(\%AA_sums);

#7. print to out 
&print_to_out($in, $out, \@aas, \%AA_sums, \%AA_counts, \%mean_ASA, $totalASA);

close(LOG);
exit;


####################
### Subroutines  ###
####################
#read in a pqsatoms XMAS file
#outputs hash of residues (key is uc(chain):resnum), value is res_type:absASA
sub res_types_from_xmas
{
    my ($pqs_id, $pqs_chain, $xmas_dir) = @_;
    my %aatypes = ();

    my $inatoms_datatype = 0;
    my $inmolecules = 0;
    my $type = 0;
    my $xmas_chain;
    
    my $xmas_file = "$xmas_dir/$pqs_id$pqs_chain.xmas";
    open (XMAS, "$xmas_file") || die "sub res_types_from_xmas cannot open $xmas_file!!!\n";
        
    while (my $line = <XMAS>)	
    {
        chomp $line;

        my $line = lc $line;
        
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
                if ( ($line =~ /^<residue>\s+(.*)\s+(.*)\s+(.*)\s+(.*)\s+(.*)\s+(.*)\s+(.*)\s+(.*)\s+<\/residue>/) &&
                     ($xmas_chain eq lc $pqs_chain ) )
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
#extracts intf residues from $intffile, outputs in A:354->ARG:79.175 format
#A=pqs chain, 354=resnum, ARG=res type, 79.175=absASA in monomer
sub find_intf_residues
{
    my ($pqs_id, $pqs_ch, $intffile, $ref_reshash) = @_;
    my %intf_res = ();

    my $idch = $pqs_id.$pqs_ch;

    open (INTF, "$intffile") || die "[find_intf_residues] cannot open $intffile!!!\n";
    while (my $line = <INTF>)
    {
        chomp $line;
        if ($line =~ /$idch/)
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
                        my ($xmas_aa_type, $absASA) = split(':', $value);
                        if ($xmas_aa_type ne $aa) 
                        {
                            print LOG "for residue $pqs_id$ch:$num aa_type in intf file is $aa and in xmas file is $xmas_aa_type!!!\n";
                        }
                        else        
                        {
                            $intf_res{"$ch:$num"} = $value;
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
#prints usage
sub Usage
{
    print "USAGE:\nperl prepare_for_propensities.pl -in <infile> -xmas_dir <xmas_dir> -out <outfile> [-v_flag] [-help]\n\n";
    print "\t-infile contains interface/surface residues, listed by PQS chains\n";
    print "\t i.e. /home/bsm/anya/intf_prediction/datasets_new/my_dataset/nonred2.intf.residues\n";
    print "-outfile will have: infile, summary per protein type [type:sum:count:mean], and total ASA area in the dataset\n";
    print "\t-v_flag used for testing...\n";
    print "\t-help will display this message\n";
    print "\nThis code prepares numners for calculate.propensities.pl (see usage of that code)\n";
    print "For a given dataset (intf or surface) will calculate:\n";
    print "\t-sum of solvent accessible surfaces for every AA type [sum] in infile\n";
    print "\t-count of residues of a AA type [count] in infile [count]\n";
    print "\t-average absASA for every AA type [mean=sum/count]\n";
    print "\t-total ASA for all AAs in infile\n";
}

##############################################################
#input hash is in format A:354->ARG:79.175
#adds ASA values to $sums_hash_ref, res. type counts to $counts_hash_ref
sub ASA_per_AA_type
{
    my ($data_hash_ref, $sums_hash_ref, $counts_hash_ref) = @_;

    foreach my $k (keys %$data_hash_ref)
    {
        my $a = "\L$$data_hash_ref{$k}";
        my ($aa_type, $asa) = split(":", $a);
        
        $$sums_hash_ref{$aa_type} += $asa;
        $$counts_hash_ref{$aa_type} ++;

        $totalASA += $asa;
    }

    return $totalASA;
}


##############################################################
#inputs sums of ASA and residue count, per AA type
#outputs average ASA per AA type
sub calculate_mean_ASA
{
    my ($sums_hash, $counts_hash) = @_;
    my %averages = ();
    my %prop_indices = ('ala' => '0',
                        'cys' => '1',
                        'asp' => '2',
                        'glu' => '3',
                        'phe' => '4',
                        'gly' => '5',
                        'his' => '6',
                        'ile' => '7',
                        'lys' => '8',
                        'leu' => '9',
                        'met' => '10',
                        'asn' => '11',
                        'pro' => '12',
                        'gln' => '13',
                        'arg' => '14',
                        'ser' => '15',
                        'thr' => '16',
                        'val' => '17',
                        'trp' => '18',
                        'tyr' => '19');

    foreach my $k (keys %$sums_hash)
    {
        if (!$$counts_hash{$k})
        {
            $averages{$k} = "N/A";
        }
        else
        {
            my $mean = $$sums_hash{$k} / $$counts_hash{$k};
            $averages{$k} = $mean;
        }
    }

    return %averages;
}

##############################################################
#inputs hash containing sum of areas for every standard AA type
#outputs their sum -> total area of residues in the dataset
sub total_area
{
    my $h_ref = $_[0];
    my $total = 0;

    foreach my $k (keys %$h_ref)
    {
        $total += $$h_ref{$k};
    }

    return $total;
}

##############################################################
#prints the summary to out
sub print_to_out
{
    my ($infile, $outfile, $r_AAarray, $r_sums, $r_counts, $r_means, $total) = @_;

        open(OUT, ">$outfile") || die "Cannot open OUT!!!\n";

        print OUT "input is $infile.\n";
        
        foreach my $AAtype (@$r_AAarray)
        {
            printf OUT "[type:sum:count:mean]%s:%12.2f:%6d:%10.4f\n", $AAtype, $$r_sums{$AAtype}, $$r_counts{$AAtype}, $$r_means{$AAtype};
        }
        printf OUT "total ASA for the dataset is %.4f.\n", $total;

        close(OUT);
    }
