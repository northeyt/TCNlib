#!/acrm/usr/local/bin/perl -w

use strict;
use DBI;
use Getopt::Long;  
use lib '/home/bsm/anya/perllib'; #sets dir where Perl_paths.pm is
use Perl_paths; 

#### Define command line inputs ####
my $patch_dir;
my $out_dir;
my $log_dir;
my $aln_dir;         
my $help = '';
my $v_flag = '';


#### Checking command line options ####
#in is dir with <pqs_id>.patches files,
GetOptions('patch_dir=s' => \$patch_dir, 
           'aln_dir=s' => \$aln_dir,
           'log_dir=s' => \$log_dir,
           'out_dir=s' => \$out_dir, 
           'help' => \$help,
           'v_flag' => \$v_flag);
#i.e. /acrm/usr/local/bin/perl calc_patch_scorecons.pl -patch_dir /acrm/home/anya/interface_prediction/patches/patches_11/ -aln_dir /acrm/home/anya/interface_prediction/fosta_alignments/ -log_dir log_files/ -out_dir scorecons/fosta/ -v_flag

if ( (!$patch_dir) || (!$log_dir) || 
     (!$out_dir) || (!$aln_dir) )
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
my $dbname = "pdbsws";
my $dbserver = 'acrm8';
my $datasource = "dbi:Pg:dbname=$dbname;host=$dbserver";
my $pdbswsdbh;

$pdbswsdbh = DBI->connect($datasource) || die "calc_patch_scorecons.pl: Cannot connect to $dbname database!!!\n";

#### main ####
##############
my @p_filelist = `ls $patch_dir`;

#my $p_file = "3hfzB.patches";
foreach my $p_file (@p_filelist)
{
    chomp $p_file;
    $p_file =~ /(.*)\.patches/;
    my $pqs_id = $1;

    my $logfile = "$log_dir/$pqs_id.scorecons.f.log";
    open(LOG, ">$logfile") || die "calc_patch_scorecons: $logfile not opened\n";

    print "TEST: grepping for $pqs_id\n";
    
    my $msa_file = `ls $aln_dir | grep $pqs_id`;
    if ($msa_file)
    {
        chomp $msa_file;
        my $chain = chop($pqs_id);

        #need to map pdb numbering to sprot numbering -> in a hash
        #hash is in 3fbi_1:B:43->P38633:40 format
        my %pdb2Sprot = &pdb2Sprot_numbering($pqs_id, $chain, $pdbswsdbh);

        printf "[pdb2sprot]For %s%s found %d residues mapped to Sprot.\n", $pqs_id, $chain, scalar keys %pdb2Sprot if ($v_flag);

        #need to map sprot numbering to msa column numbering -> in a hash
        my %pdb2Msa = &pdb2Msa_numbering($pqs_id, $chain, $aln_dir, \%pdb2Sprot);

        printf "[pdb2msa]For %s%s found %d PDB residues mapped to MSA file.\n", $pqs_id, $chain, scalar keys %pdb2Msa if ($v_flag);

        #need to extract all scorecons values per msa columns -> in a hash
        my %msa2scorecons = &msa2scorecons_numbering($pqs_id, $chain, $aln_dir);
        
        #calculate average seq.cons. score per patch 
        #then write that down in an outfile in out_dir 
        my $patchfile_path = "$patch_dir/$p_file";
        open(PFILE, "$patchfile_path") || die "calc_patch_scorecons (PFILE handle): $p_file not opened\n";
        
        my $outfile = "$out_dir/$pqs_id.$chain.patch.scorecons";
        open(OUT, ">$outfile") || die "calc_patch_scorecons (OUT handle): $outfile not opened\n";

        print "TESTING: $outfile\n";
        
        while (my $line = <PFILE>)
        {
            chomp $line;
            if (length($line))
            {
                $line =~ /patch\s(\S+)>\s(.*)$/;
                my $central_res = $1;
                my $patch_r = $2;
                my @patch_residues = split(' ', $patch_r);
                
                my $avg_cons_score = &calc_patch_cons_score($pqs_id, $chain, $central_res, \@patch_residues, \%pdb2Msa, \%msa2scorecons);
                print OUT "<patch $central_res> $avg_cons_score\n";
            }
        } 
        close(OUT);
        close(PFILE);
    }
    else
    {
        print "TESTING: NO MSA files found for $pqs_id!!!\n";
    }

    close(LOG);

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

$pdbswsdbh->disconnect();
exit;

##############################################################
# SUBROUTINES
##############################################################
sub pdb2Sprot_numbering
{
    my ($pqs_id, $ch, $dbhandle) = @_;
    my %pdb2Sprot = ();

    my $pdb_id = substr($pqs_id, 0, 4);

    #resid is resnum, pdbaa is 1-letter, resnam is 3-letter (both upper case)
    my $sql = "SELECT resid, pdbaa, ac, swsaa, swscount
               FROM alignment
               WHERE pdb = '$pdb_id'
               AND chain = '$ch';";

    my $pdbswssth = $dbhandle->prepare($sql);
    if($pdbswssth->execute)
    {
        my ($pdb_resnum, $pdb_aa, $sprot_ac, $sprot_aa, $sprot_resnum);
        while(($pdb_resnum, $pdb_aa, $sprot_ac, $sprot_aa, $sprot_resnum) = $pdbswssth->fetchrow_array)
        {
            #we currently ignore if sprot and pdb res. types are not the same!
            $pdb2Sprot{"$pqs_id:$ch:$pdb_resnum:$pdb_aa"} = "$sprot_ac:$sprot_resnum:$sprot_aa";

            print "$pqs_id:$ch:$pdb_resnum:$pdb_aa->$sprot_ac:$sprot_resnum:$sprot_aa\n" if ($v_flag);
        }
    }

    return %pdb2Sprot;
}


##############################################################
sub pdb2Msa_numbering
{
    my ($pqs_id, $ch, $aln_dir, $r_pdb2sprot_hash) = @_;
    my %pdb2Msa = ();

    my $idch = $pqs_id.$ch;
    my $msa_file = `ls $aln_dir | grep $idch`;

    chomp $msa_file;
    $msa_file =~ /^$idch\.(.*)$/;
    my $sprot_id = $1;

    my $string = &extract_aln_type($aln_dir);

    #my $full_path = $aln_dir.$msa_file;
    my $aligned_sequence = &extract_aligned_sequence($aln_dir, $msa_file, $string, $sprot_id);
    printf "for %s found aligned sequence:\n%s.\n", $sprot_id, $aligned_sequence if($v_flag); 

    my %sprot2msa = &sprot2Msa_numbering($aligned_sequence);

    printf "[pdb2msa]For %s%s found %d residues mapped from Sprot to MSA file.\n", $pqs_id, $ch, scalar keys %sprot2msa if ($v_flag);
    
    #combining pdb2sprot and sprot2msa mapping
    foreach my $pdb_res (sort keys %$r_pdb2sprot_hash)
    {
        my ($pdb_id, $chain, $pdb_resnum, $pdb_aa) = split(':', $pdb_res);

        #rest is 'sprot_resnum:sprot_aa'
        my ($ac, $rest) = split(':', $$r_pdb2sprot_hash{$pdb_res}, 2);
        my ($sprot_resnum, $sprot_aa) = split(':', $rest);

        if( ($sprot2msa{$rest}) && ($pdb_aa eq $sprot_aa) )
        {
            print LOG "$sprot2msa{$rest}, pdb_aa $pdb_aa, sprot_aa $sprot_aa\n";
            #losing the res type. not needed hereafter
            chop $pdb_res;
            chop $pdb_res;
            $pdb2Msa{$pdb_res} = $sprot2msa{$rest};
        }
        elsif ($pdb_aa ne $sprot_aa)
        {
            print LOG "For pdb:$pdb_res, PDBSWS has restype $sprot_aa!\n";
        }
        else
        {
            print LOG "PDBSWS mapped pdb_res->$$r_pdb2sprot_hash{$pdb_res}, but $rest not found in MSA file!!!\n";
        }
    }

    return %pdb2Msa;
}


##############################################################
sub extract_aligned_sequence
{
    my ($dir, $file, $msa_type, $id) = @_;
    my $path = "$dir/$file";
    open(MSA, "$path") || die "calc_patch_scorecons (MSA handle): $path not opened\n";
    
    my $in_seq = 0;
    my $sequence = "";

    #in fosta, fasta header of the target sequence is '>P38633|MED31_YEAST', filename is 3fbi_1B.MED31_YEAST
    #in psiblast, fasta header of the target sequence is '>3fbi_1B', filename is 3fbi_1B
    my $seq_header_regex;
    if ($msa_type =~ /fosta/)
    {
        $seq_header_regex = "$id";
    }
    if ($msa_type =~ /blast/)
    {
        $seq_header_regex = "^>$file";
    }

    while (my $line = <MSA>)
    {
        chomp $line;

        if ($in_seq)
        {
            if ($line !~ /^>/)
            {
                $sequence = $sequence.$line;
            }
            else
            {
                $in_seq = 0;
            }
        }
        #works only for the FOSTA sequences...
        if ($line =~ /$seq_header_regex/)
        {
            $in_seq = 1;
        }
    }

    close(MSA);
    return $sequence;
}


##############################################################
sub sprot2Msa_numbering
{
    my $sequence = $_[0];
    my %sprot2msa = ();
    my $msa_pos_count = 0; #counts alignment columns
    my $seq_pos_count = 0; #counts positions in the sequence
    
    my @chars = split(//, $sequence);

    foreach my $c (@chars) 
    {
        $msa_pos_count++;

        if($c ne '-')
        {
            $seq_pos_count++;
            $sprot2msa{"$seq_pos_count:$c"} = $msa_pos_count;
            print "[sprot2msa]$seq_pos_count:$c->$msa_pos_count\n" if ($v_flag);
        }
    }

    return %sprot2msa;
}


##############################################################
sub msa2scorecons_numbering
{
    my ($pqs_id, $ch, $msa_dir) = @_;    
    my %msa2scorecons = ();

    my $idch = $pqs_id.$ch;
    my $msa_file = `ls $aln_dir | grep $idch`;
    chomp $msa_file;

    my $infile = "$msa_dir/$msa_file";
    my $outfile = $pqs_id.$ch.".fosta.scorecons";

    `$Perl_paths::scoreconsExe $infile $Perl_paths::valdar01_params -o $outfile`;
    open(CONS, "$outfile") || die "calc_patch_scorecons (CONS handle): $outfile not opened\n";
    my @scorecons = <CONS>;
    close(CONS);
    `rm -f $outfile`;

    my $line_count = 0;

    foreach my $line (@scorecons)
    {
        chomp $line;

        $line_count++;
        my ($score, $char, $column) = split(/\s+/, $line);
        $msa2scorecons{$line_count} = $score;
    }

    return %msa2scorecons;
}


##############################################################
sub calc_patch_cons_score
{
    my ($pqs_id, $chain, $central, $r_patch_res_array, $r_pdb2Msa_hash, $r_msa2cons_hash) = @_;

    #calculated as averaged scorecons score of all residues in a patch
    my $patch_cons_score = "N/A";
    my $score_sum = 0;
    my $res_count = 0;

    foreach my $res (@$r_patch_res_array)
    {
        #print "res is $res.\n";
        my $key = "$pqs_id:$res";
        my $msa_column = $$r_pdb2Msa_hash{$key};
        
        if (exists $$r_msa2cons_hash{$msa_column})
        {
            $res_count++;
            $score_sum += $$r_msa2cons_hash{$msa_column};
        }
    }

    if (!$res_count)
    {        
        print LOG "for $pqs_id:$chain patch $central, no residue cons. scores were found!!!\n";
    }
    else
    {
        $patch_cons_score =  $score_sum/$res_count;
    }

    return $patch_cons_score;
}


##############################################################
sub Usage
{
    print"\nUSAGE:\n";
    print"\t/acrm/usr/local/bin/perl calc_patch_scorecons.pl -patch_dir dir -aln_dir dir -log_dir dir -out_dir outfile\n\n"; 
    print"i.e.   /acrm/usr/local/bin/perl calc_patch_scorecons.pl -patch_dir /acrm/home/anya/interface_prediction/patches/patches_11/ -aln_dir /acrm/home/anya/interface_prediction/fosta_alignments/ -log_dir log_files/ -out_dir scorecons/fosta/ -v_flag \n\n";
    print "This code will produce a file for every patch-containing pqs_idch\n"; 
    print "for every patch in that pqs_id it will output a line:\n";
    print "\t-<patch A.105> X\n";
    print "\twhere A.105 is the central residue for that patch\n";
    print "\tX is the valdar01 scorecons sequence conservation score, based on provided alignment in aln_dir\n";
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
#returns 'fosta' or 'blast' or exits the script 
#based on the type of alignments in aln_dir
sub extract_aln_type
{
    my $aln_dir = $_[0];
    my $aln_type;
    if($aln_dir =~ /fosta/)
    {
        $aln_type = 'fosta';
    }
    elsif($aln_dir =~ /blast/)
    {
        $aln_type = 'blast';
    }
    else
    {
        die "aln_dir neither of FOSTA nor BLAST type!!! Exiting...\n";
    }

    return $aln_type;
}
