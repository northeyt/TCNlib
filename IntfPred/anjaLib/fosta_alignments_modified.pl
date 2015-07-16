#!/acrm/usr/local/bin/perl -w

use strict;
use lib '/home/bsm/anya/perllib'; #sets dir where Perl_paths.pm is
use Perl_paths; #all exe paths are set here
use DBI;
#use POSIX;
use Getopt::Long;       # for command line options

#############

### define command line inputs ####
my $in_dir;          #dir containing all patch files
my $out_dir;         #dir containing MSA files in clustal format
my $log;             #writes the log file (why there is no MSA file for that id) - only for sprot_IDs with no output files
my $out;             #for each patch file (a pdb chain), writes mapped sprot AC and ID, and #FEPS, if any
my $seq_num = 5;     #default minimum 5 sequences to build an alignment
my $help = 0;        #prints usage


#### variables and handles ####
my $dbname2 = "pdbsws";
my $dbserver2 = 'acrm8';
my $datasource2 = "dbi:Pg:dbname=$dbname2;host=$dbserver2";
my $pdbswsdbh;
my $pdbswssth;

my $dbname3 = "fosta";
my $dbserver3 = 'acrm8.biochem.ucl.ac.uk';
my $datasource3 = "dbi:Pg:dbname=$dbname3;host=$dbserver3";
my $fostadbh;
my $fostasth;


#### Checking command line options ####
GetOptions('in_dir=s' => \$in_dir, 'out_dir=s' => \$out_dir, 'out=s' => \$out, 'log=s' => \$log, 'seq_num' => \$seq_num, 'help' => \$help);

if ( (!$in_dir) || (!$out_dir) || 
     (!$out) || (!$log) )
{
    $help = 1;
}

if ($help)
{
    &Usage();
    exit;
}

#### preparing files and handles ####
open(RECORD, ">$log") || die "fosta_alignments.pl: $log not opened\n"; #fosta_alignments.rec or .log
open(OUT, ">$out") || die "fosta_alignments.pl: $out not opened\n"; #fosta_alignments.out

$pdbswsdbh = DBI->connect ($datasource2) || die "fosta_alignments.pl: Cannot connect to $dbname2 database.\n";
$fostadbh = DBI->connect ($datasource3) || die "fosta_alignments.pl: Cannot connect to $dbname3 database.\n";

my $rel_count = 0; #counts unreliable sprot_ids that PDB chains get mapped to

#### main ####
#my @id_list = `ls $Perl_paths::patch_dir`;
my @id_list = `ls $in_dir`;

#my $file = "1a17A.patches";
foreach my $file (@id_list)
{
    chomp $file;
    $file =~ /^(\S+)\.patches$/;

    my $pqs_id = $1;
    my $chain = chop($pqs_id);
    my $pdb_id = substr($pqs_id, 0, 4);
    
    my $sql = "SELECT ac
               FROM pdbsws
               WHERE pdb = '$pdb_id'
               AND chain = '$chain'
               AND valid = 't'
               AND aligned = 't'
               AND ac != 'SHORT'
               AND ac != 'DNA'
               AND ac != 'ERROR';";

    my $sprot_ac = $pdbswsdbh->selectrow_array($sql);#&PDBch_to_sprotAC($pdb_id, $chain, $pdbswsdbh);

    if ($sprot_ac)
    {
        my $sql = "SELECT i.id
                   FROM idac i, acac a
                   WHERE a.altac = '$sprot_ac'
                   AND i.ac = a.ac;";
        my $sprot_id = $pdbswsdbh->selectrow_array($sql);
        #my $sprot_id = &PDBch_to_sprotAC($sprot_ac, $pdbswsdbh); 

        if ($sprot_id)
        { 
            print "PDB: $pdb_id$chain, SwissProt: $sprot_ac, $sprot_id ";
            print OUT "PDB: $pdb_id$chain, SwissProt: $sprot_ac, $sprot_id ";

            #query FOSTA for all FEPs of that sprot ID
            #$sprot_id is included
            my @FEParray = &FindFEP($sprot_id);
            
            #eliminates if no FEPs were found (also eliminates unreliable)
            if ($#FEParray>0) #index of the last element of an array 
            {    
                my %sequences = &retrieve_sequences(@FEParray);
                my $seq_count = scalar keys %sequences;
                
                #at least 5 sequences in alignment...
                if ($seq_count < $seq_num) 
                {
                    print "\n";
                    print OUT "\n";
                    print RECORD "$sprot_id had <5 ($seq_count) sequences in the alignment.\n";
                }
                else
                {
                    my $fasta_file = &make_fasta_file($pqs_id, $chain, \%sequences);
                    my $aln_file = &do_muscle_alignment($fasta_file, $sprot_id, $out_dir);
                    printf "%d seqs in family.\n", scalar keys %sequences;
                    printf OUT "%d seqs in family.\n", scalar keys %sequences;
                }
            }
            else
            {
                print "\n";
                print OUT "\n";
            }
        }
        else
        {
            print RECORD "$sprot_ac had no matching IDs in pdbsws.\n"; 
        }                
    }
    else
    {
        print RECORD "$pdb_id had no matching ACs in pdbsws.\n";
    }       
}

$pdbswsdbh->disconnect();
$fostadbh->disconnect();
close(RECORD);
close(OUT);

exit;


########################################################################
#SUBROUTINES
########################################################################
#checks whether entry is unreliable for orthofind and finds functionally equivalent proteins (FEPs) for reliable query proteins
#returns array, first element is query id, others are its FEPs
sub FindFEP
{
    my $id = $_[0];
    my @array = ();  
    my $fep_id;

    chomp $id;
    if ($id)
    {
        my $sql3 = "SELECT fd.unreliable, f.fosta_family
                    FROM fosta_descriptions fd, feps f
                    WHERE fd.id='$id'
                    AND f.id='$id'";

        my ($unreliable, $family) = $fostadbh->selectrow_array($sql3);
       
        if ($family)
        {
            #$unreliable becomes 0 or 1, although postgres shows t or f
            if (!$unreliable)
            {    
                push (@array, $id);                
                my $found_reliable_feps4reliable_id = 0;

                #ask Lisa is this query OK
                my $sql4 = "SELECT id
                            FROM feps
                            WHERE fosta_family = '$family'
                            AND id != '$id'
                            AND NOT unreliable;";
                
                my $fostasth = $fostadbh->prepare($sql4);
                if($fostasth->execute)
                {
                    while(($fep_id) = $fostasth->fetchrow_array)
                    {
                        push (@array, $fep_id);                        
                        $found_reliable_feps4reliable_id = 1;
                    }
                }
                
                if (!$found_reliable_feps4reliable_id)
                {
                    print RECORD "$id was reliable but had no reliable feps!\n";
                }
            }
            else
            {
                $rel_count++;
                print RECORD "$id is unreliable query protein.\n";
            }
        }
        else
        {
            print RECORD "$id not found in feps table.\n";
        }
    }  
    else 
    {
         print RECORD "FindFEP could not read query id.\n";
    }

    return @array;
}


#############################################
#for given id, puts acc, id and sequence in %seqs
sub retrieve_sequences 
{
    my @seq_ids = @_;
    my %seqs;
    my $seq_id;

    foreach $seq_id (@seq_ids)
    {               
        my $seq_sql = "SELECT sequence
                       FROM fosta_sequences
                       WHERE id='$seq_id'";
        my $seq = $fostadbh -> selectrow_array($seq_sql);

        my $acc_sql = "SELECT a.ac 
                       FROM idac i, acac a
                       WHERE i.id='$seq_id'
                       AND i.ac = a.altac";
        my $acc = $pdbswsdbh -> selectrow_array( $acc_sql );

        if ($seq && $acc)
        {
            $seqs{ "$acc|$seq_id" } = $seq;            
        }    
        else
        {
            print RECORD "$seq_id had no sequence found\n";
        }
    }
   
    return %seqs;
}


#############################################
#makes a fasta file for a given sequence
sub make_fasta_file
{
    my ($id, $ch, $ref_sequences) = @_;
    
    #/tmp/fosta1ek6A.fasta
    my $fasta_file = $Perl_paths::tmp_dir."fosta".$id.$ch.".fasta";

    #my $file = &tmp_file("intf_fosta", ".fasta", "/tmp");
    open (TMP, ">$fasta_file") || die "Cannot write to $fasta_file!!!\n";
   
    while (my ($seq_id, $seq) = each (%$ref_sequences))
    {
        print TMP ">$seq_id\n$seq\n";
    }

    close (TMP);
    return $fasta_file;
}


#############################################
#saves a file containing an alignment of a query protein and its FEPs in fosta_alignments directory
sub do_muscle_alignment
{
    my ($fasta_file, $id, $out_dir) = @_;

    $fasta_file =~ /fosta(.*)\.fasta/;
    my $PQS_idch = $1;
    my $outfile = "$out_dir/$PQS_idch.".".$id";
    #my $outfile = "/home/bsm/anya/intf_prediction/code/alignments/fosta_alignments/".$id;

    `/home/bsm/anya/tools/muscle3.7/muscle -in $fasta_file -out $outfile -maxiters 100 -stable -quiet`;
    #`/home/bsm/anya/tools/muscle3.7/muscle -in $fasta_file -out $outfile -maxiters 100 -clw -stable -quiet`;

    return $outfile;
}


#############################################
#returns a sprot_ac mapping for a pdb chain
sub PDBch_to_sprotAC
{
    my ($PDB_id, $ch, $dbhandle) = @_;
    
    #maps PDB chain (or residue) to a swissprot ac
    my $sql = "SELECT ac
               FROM pdbsws
               WHERE pdb = '$PDB_id'
               AND chain = '$ch'
               AND valid = 't'
               AND aligned = 't'
               AND ac != 'SHORT'
               AND ac != 'DNA'
               AND ac != 'ERROR';";
    if ($dbhandle)
    {
        print "yes\n";
    }
    else                        
    {
        print "no\n";
    }

    my $sprot_ac = $dbhandle->selectrow_array($sql);

    return $sprot_ac;
}

#############################################
#returns a sprot_id mapping for a sprot_ac
sub sprotAC_to_sprotID
{
    my ($ac, $dbhandle) = @_;
 
    #maps swissprot ac to a swissprot id
    my $sql = "SELECT i.id
               FROM idac i, acac a
               WHERE a.altac = '$ac'
               AND i.ac = a.ac;";
    my $sprot_id = $dbhandle->selectrow_array($sql);

    return $sprot_id;
}

#############################################
#prints help message
sub Usage
{
    print "\nUSAGE:\nperl fosta_alignments.pl -in_dir <indir> -out_dir <outdir> -out <outfile> -log <logfile> [-seq_num <int>]\n\n";
    print "\t-out_dir is dir containing patch files (one file per PQS chain)\n";
    print "\t-in_dir is dir containing MSA files in clustal format,\n";
    print "\t-out for each patch file records mapped sprot AC and ID, and #FEPS\n";
    print "\t-log records error messages \n"; 
    print "\t(possible messages: (i)pdb_id had no matching ACs in pdbsws, (ii)sprot_ac had no matching IDs in pdbsws, (iii)sprot_id not found in feps table, (iv)sprot_id had no sequence found)\n"; 
    print "\t-seq_num sets the minimum number of sequences for the alignment, \n\t(if not specified, default is 5) \n\n"; 
}
