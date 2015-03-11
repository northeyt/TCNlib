#!/acrm/usr/local/bin/perl -w

use strict;
use DBI;
use lib '/home/bsm/anya/perllib'; #sets dir where Perl_paths.pm is
use Perl_paths; #all hard-codes are here
use Getopt::Long;       # for command line options

#### define command line inputs ####
my $in_dir;
my $pdb_dir;
my $log;
my $out_dir;
my $seq_num = 5; #min alignment of 5 seqs
my $cont;        #flag to continue, instead of starting from the scratch
my $help = 0;

GetOptions( 'in_dir=s' => \$in_dir, 'pdb_dir=s' => \$pdb_dir, 'out_dir=s' => \$out_dir, 'log=s' => \$log, 'seq_num' => \$seq_num, 'cont' => \$cont, 'help' => \$help);


#### variables and handles ####
my @files;

my $dbname = "pdbsws";
my $dbserver = 'acrm8';
my $datasource = "dbi:Pg:dbname=$dbname;host=$dbserver";
my $pdbswsdbh;
#my $pdbswssth;

$pdbswsdbh = DBI->connect ($datasource) || die "psiblast_alignments_per_PQSch.pl: Cannot connect to $dbname database.\n";

#### preparing stuff ####
#if command line input is not defined
if ( (!$in_dir) || (!$out_dir) || 
     (!$log) || ! $pdb_dir )
{
    $help = 1;
}

if ($help)
{
    &Usage();
    exit;
}

#list of files to process...
if ($cont)
{
    my @files = &eliminate_done_files($in_dir, $out_dir);

    printf "[after -cont] found %d files to process.\n", scalar @files;
}
else
{
    @files = `ls $in_dir`;
    printf "No -cont flag: found %d files to process.\n", scalar @files;
}

#### main ####
open(RECORD, '>', "$log") || die "psiblast_alignments_per_PQSch.pl: RECORD not opened\n";


#my $patch_file = "4ubpC.patches";
#patch file is like '1a2xA.patches'
foreach my $patch_file (@files)
{
    chomp $patch_file;
    print "[psiblast] $patch_file\n";

    $patch_file =~ /^(\S+)\.patches$/;
    my $pqs_idch = $1;  
    my $pdb_id = substr($pqs_idch, 0, 4);
    
    #extract seqres sequence for the PDB file

    my $single_chain_fasta = &retrieve_chain_seq($pqs_idch, $pdb_dir);

    if (&outfile_not_empty($single_chain_fasta) == 0)
    {                              
        print RECORD "[$pqs_idch] no fasta file for that chain!\n";
    }
    else
    {
        #extract sequence from the fasta format of pqs chain
        my $query_seq = &extract_query_chain_sequence($single_chain_fasta);
                
        #runs blastp
        my $blast_result = &blast_single_chain($single_chain_fasta);
        
        #returns e-val and description 
        my %blast_hits = &extract_blast_hits($blast_result); 
        my $blast_hit_count = scalar keys %blast_hits;
    
        if($blast_hit_count)
        {
            my %seqs_to_align = &seqs_of_filtered_hits(\%blast_hits,$pdbswsdbh);
            my $filt_blast_hits = scalar keys %seqs_to_align;
     
            if ($filt_blast_hits)
            {
                #reduce number of hits to 200
                if ($filt_blast_hits > 200)
                {
                    %seqs_to_align = &top_200_blast_hits(\%seqs_to_align, \%blast_hits);
                }
            
                #adds pqs chain sequence to the hash, it is often already one of the hits
                $seqs_to_align{$pqs_idch} = $query_seq;       
                $filt_blast_hits = scalar keys %seqs_to_align;
            
                if($filt_blast_hits < $seq_num) 
                {
                    print RECORD "[$pqs_idch] Not enough seqs for the alignment ($filt_blast_hits, orig. seq. included)\n";
                }
                else
                {
                    my $fasta_file = &make_fasta_file($pqs_idch, \%seqs_to_align);
                    
                    #creates alignment file in $out_dir
                    my $aln_file = &do_muscle_alignment($fasta_file, $pqs_idch, $out_dir);     
                    
                    #returns 1 if outfile is not empty
                    if (&outfile_not_empty($aln_file) == 0)
                    {
                        print RECORD "$pqs_idch created an empty alignment file!!!\n";
                    }
                }        
            }
            else
            {
                print RECORD "[$pqs_idch] No BLAST hits (after reliability and e-value filtering)!!!\n";
            }
        }
        else
        {
            print RECORD "[$pqs_idch] BLAST returned no hits (before filtering)\n";
        }
    }
}

$pdbswsdbh->disconnect();
close(RECORD);
exit;


########################################################################
#SUBROUTINES
########################################################################
#retrieves seqres sequence from PQS file (one chain only)
#all residues are uppercase (-u), keys are chains, values are sequences per chain
#10.11.2009: will extract sequence per chain - after getchain, 
#            also uses PQS entries, old sub retrieve_seqres was PDB-based
#            outputs sequence in a fasta format
sub retrieve_chain_seq
{
    my ($pqs_id, $pdb_dir) = @_;
    my $seq = "";

    my $chain = chop($pqs_id); 
    
    #need PQS sequences!!!!
    my $in = "$pdb_dir/" . uc( $pqs_id . "_$chain" ) . '.pdb';

    my $ch_pqs_file = $pqs_id.$chain.$Perl_paths::pqsext;

    #extracts single-chain PQS entry
    my $getchain = $Perl_paths::martin_bin."getchain";
    `$getchain -a $chain $in $ch_pqs_file`;

    my $out = $pqs_id.$chain.$Perl_paths::fastaext;
    
    #single chain sequence in $out 
    my $pdb2pir = $Perl_paths::martin_bin."pdb2pir";
    # -i option required to extract only residues from chain
    `$pdb2pir -f -s -u -i $ch_pqs_file $out`;
    `rm -f $ch_pqs_file`;

    return $out;
}


#############################################
#saves a file containing an alignment of a query protein and its FEPs in fosta_alignments directory
sub do_muscle_alignment
{
    my ($fasta_file, $id, $outdir) = @_;
    my $outfile = "$outdir/$id";

    `$Perl_paths::muscleExe -in $fasta_file -out $outfile -maxiters 100 -stable -quiet`;
    #`$Perl_paths::muscleExe -in $fasta_file -out $outfile -maxiters 100 -clw -stable -quiet`;

    `rm $fasta_file`;
    return $outfile;
}


#############################################
#blasts one PQS chain sequence against nr, using LCR filter
#output in XML-like format
sub blast_single_chain
{
    my $in = $_[0];

    $in =~ /(\S+)\.fasta/;
    my $blast_xml = $1.".blast";
 
    #/acrm/data/tmp/fosta/swissprot/uniprot_sprot.fasta is mirrored latest version of SwissProt, and formatdb was performed on it
    #-b 2000 -v 2000 to mimic FOSTA code
    `$Perl_paths::BlastAllExe -p blastp -m 7 -b 2000 -v 2000 -F T -d $Perl_paths::sprot_fasta_mirror -i $in -o $blast_xml`;
    
    #removes single_chain_fasta file
    `rm -f $in`;

    return $blast_xml;
}


#############################################
#extracts every <Hit> tag and its description and e-value
sub extract_blast_hits
{
    my $blast_file = $_[0];
    my $in_hit = 0; 
    my %blast_hits = ();
    my $ac = "N/A";
    my $id = "N/A";
    my $e = "N/A";
    my $def = "N/A";

    open (BL, "$blast_file") || die "extract_blast_hits: Cannot open $blast_file\n";

    while (my $line = <BL>)
    {
        chomp $line;
        #tag open
        if ($line =~ /\s+<Hit>/)
        {
            $in_hit = 1;
        }
        elsif ($line =~ /\s+<\/Hit>/)
        {
            $in_hit = 0;

            #store blast hit details in the hash
            my $key = $ac.":".$id;
            $blast_hits{$key} = "$e->$def";

            #reset labels
            $ac = "N/A";
            $id = "N/A";
            $e = "N/A";
            $def = "N/A";
        }
        else
        {
            if ($in_hit)
            {
                if ($line =~ /\s+<Hit_def>(.*)<\/Hit_def>/)
                {
                    $def = $1;
                    my @begin = split(/\s/, $def);
                    $begin[0] =~ /sp\|(\S{6})\|(\S+)/;
                    $ac = $1;
                    $id = $2;
                }
                elsif ($line =~ /\s+<Hsp_evalue>(.*)<\/Hsp_evalue>/)
                {
                    $e = $1;
                }
                else
                {
                    ;
                }
            }
        }
    }
    close(BL);
    `rm -f $blast_file`;

    return %blast_hits;
}


#############################################
#inputs e-value in 1.3e-13 format
#returns 1 if e-value is <=0.01
sub e_val_filter
{
    my $e = $_[0];

    if ($e =~ /e/) #e-value in exp form
    {
        $e =~ /([0-9\.]+)e\-([0-9]+)/;
        my $f = $1; #factor multiplying
        my $exp = $2; #negative exponent of base 10

        if ( ($exp > 2) ||
             ( ($exp == 2) && ($f <= 1.00000) ) )
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else           #e-value in decimal form
    {
        if ($e <= 0.01000)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }

    return 0;
}


#############################################
#inputs description of a blast hit
#returns 0 if the hit is predicted/putative/hypothetical,
#else returns sprot_ac of that blast hit
sub reliable_filter
{
    my $def = $_[0];

    my ($a, $desc) = split(/\s/, $def, 2);
    $a =~ /^sp\|(\w{6})\|(\w+)$/;
    my $ac = $1;
    my $id = $2;

    $desc = "\L$desc";

    if ( ($desc !~ /putative/) &&
         ($desc !~ /predicted/) &&
         ($desc !~ /hypothetical/) )
    {
        return ($ac.":".$id);
    }        
    
    return 0;
}


#############################################
#keys are sprot_ac:sprot_id, 
#value is sequence (sprot, not segment in blast file)
#returns sequences of blast hits with e-value < 0.01 and
#not containing 'putative/predicted/hypothetical' in description

#changed on 07.06.2010. to extract seqs from pdbsws rather than from indexed fasta file (which was not working at the time)
sub seqs_of_filtered_hits
{
    my ($blast_hits_ref, $dbhandle) = @_;
    #my $blast_hits_ref = $_[0];
    my %seqs_to_align = ();
    
    foreach my $k (keys %$blast_hits_ref)
    {
        my ($e, $def) = split(/\->/, $$blast_hits_ref{$k});
        if (&e_val_filter($e))
        {
            if ( my $ac_id = &reliable_filter($def))
            {           
                my ($ac, $id) = split(':', $ac_id);
     
                #my $fasta = `$Perl_paths::PerlExe $Perl_paths::indexfasta_dir/getfasta.pl $Perl_paths::sprot_fasta_mirror $Perl_paths::sprot_index_file $id`;
                #print "$Perl_paths::PerlExe $Perl_paths::indexfasta_dir/getfasta.pl $Perl_paths::sprot_fasta_mirror $Perl_paths::sprot_index_file $id\n";

                #my @fasta_array = split(/\n/, $fasta);
                #my $sequence = '';

                #foreach my $l (@fasta_array)
                #{
                #    if ($l !~ />sp\|/)
                #    {
                #        $sequence = $sequence.$l;
                #    }
                #} 
                #primary ac-s??????

                my $sql = "SELECT sequence
                           FROM sprot
                           WHERE ac = '$ac';";
                my $sequence = $pdbswsdbh->selectrow_array($sql);
                if ($sequence)
                {
                    $seqs_to_align{$ac_id} = $sequence;
                }
                else
                {
                    print RECORD "No sequence found in PDBSWS for ac:$ac!!!\n";
                }
            }
        }      
    }

    return %seqs_to_align;
}

#############################################
#makes a fasta file for a given sequence
#$seq_id is sprot_ac:sprot_id
sub make_fasta_file
{
    my ($id, $ref_sequences) = @_;

    my $file = "blast.to_align.$id.fasta";
    open (TMP, ">$file") || die "Cannot write to file \'$file\'\n";
   
    while (my ($seq_id, $seq) = each (%$ref_sequences))
    {
        print TMP ">$seq_id\n$seq\n";
    }

    close (TMP);
    return $file;
}

######################################################
sub extract_query_chain_sequence
{
    my $fasta = $_[0];
    open (SEQ, "$fasta") || die "Cannot read $fasta!!!\n";
    my $seq = '';

    #$fasta =~ /(.*)\.fasta/;
    #my $id_chain = $1; #1a3cA format

    while (my $line = <SEQ>)
    {
        chomp $line;

        #if not in header line
        if ($line !~ /^>(.*)/)
        {            
            $seq = $seq.$line;
        }
    }
    close(SEQ);

    return $seq;
}

######################################################
#filtered hash entry is $ac:$id->$sequence
#all hash entry is $id->$e->$definition
#orders entries by e-value and then returns only first 200
sub top_200_blast_hits
{
    my ($filt_blast_hits_ref, $all_blast_hits_ref) = @_;
    my %new_hash = ();
    my %top200 = ();

    #creating a hash with $e for value to be sorted
    foreach my $k (keys %$filt_blast_hits_ref) 
    {
        #put in the e-value from hash#2
        my $a = $$all_blast_hits_ref{$k};
        $a =~ /^(.*)\->sp/;
        my $e = $1;
        $new_hash{$k} = $e;        
    }  

    #sorting %new_hash (ascending) by value, keys in array
    my @sorted = sort {$new_hash{$a} <=> $new_hash{$b}} keys %new_hash;

    #put top 200 into %top200
    for (my $i = 0; $i<200; $i++) 
    {
        $top200{$sorted[$i]} = $$filt_blast_hits_ref{$sorted[$i]};
    }

    return %top200;
}

######################################################
### Prints help message
sub Usage
{
    print "DESCRIPTION:\nIn summary, for each patch file in in_dir creates a muscle alignment of filtered blast hits. First maps a PQS chains to a sprot_ac, extracts sprot sequence, BLASTs against indexed SwissProt fasta file, extract all reliable hits wih e<0.01 (min 4, max. 200 hits) and aligns them using muscle.\n\n";
    print "USAGE:\npsiblast_alignments_per_PQSch.pl -in_dir <dirpath> -pdb_dir <dirpath> -out_dir <dirpath> -log <logfile> [-seq_num] [-cont]\n\n";
    print "Where:\n";
    print "\t-in_dir contains patch files (in '119lA.patches' format),\n";
    print "\t-MSA files are stored in out_dir\n";    
    print "\t-log has messages for all files not outpitting an alignment\n";
    print "\t-seq_num min. size of the alignment to be stored\n"; 
    print "\t(if not specified, default is 5)\n";
    print "\t-cont flag will process only patch files not already processed in out_dir\n"; 
    print "\ni.e. /acrm/usr/local/bin/perl psiblast_alignments_per_PQSch.pl -in_dir /acrm/home/anya/interface_prediction/patches/patches_11/ -out_dir /acrm/home/anya/interface_prediction/psiblast_alignments/ -log psiblast_alignments_per_PQSch.log\n";
}


########################################################################
#returns list of files to process
#indir minus outdir
sub eliminate_done_files
{
    my ($indir, $outdir) = @_;
    my @files2process;

    my @in = `ls $indir`;
    my @out = `ls $outdir`;
    
    foreach my $file (@in)
    {
        chomp $file;
        $file =~ /(.*)\.patches/;
        my $PQSidch = $1;
        my $processed = 0;
        
        foreach my $f (@out)
        {
            chomp $f;
            if ($f eq $file)
            {
                $processed = 1;
            }
        }

        if (!$processed)
        {
            push(@files2process, $file);
        }
        else
        {
            print "[$PQSidch] already processed\n";
        }
    }

    return @files2process;
}


########################################################################
### Returns 1 if outfile is not empty, 0 if it is empty
sub outfile_not_empty
{
    my $file = $_[0];
    open(F, "$file") || die "outfile_not_empty cannot open $file!!!\n";

    my @array = <F>;
    my $size = scalar @array;
    close(F);

    if ($size>0)
    {
        return 1;
    }
    return 0;
}
