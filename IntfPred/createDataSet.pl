#!/usr/bin/perl -w
# createDataSet.pl --- wrapper that runs all required scripts to prepare dataset for WEKA input
# Author: Tom Northey <zcbtfo4@acrm18>
# Created: 13 Jan 2014
# Version: 0.01

use warnings;
use strict;
use IO::CaptureOutput qw ( qxy );
use File::Spec;
use File::chdir;
use Text::Chomped;
use Getopt::Long;
use Carp;
use FindBin qw($Bin);

my $lib = "$Bin/lib/";

@ARGV or Usage();

my($patch_radius, $how_file, $chainIDsFile,
   $run_dir, $class_lab_file, $patches_dir,
   $intfStatFile, $surfStatFile, $newStatFiles,
   $setForTraining, $standardizeAgainst,
   $unsupervised)
    = process_cmdline();

my %input;

# Determine input type: either .how file,
# chainIDs file, or a patches directory
# and class label file
if ($how_file) {
    $input{opt} = "-how_file";
    $input{arg} = $how_file;
}
elsif ($class_lab_file && $patches_dir) {
    $input{opt} = "-chain_file";
    $input{arg} = createpdbIDsAgFile($patches_dir);
}
elsif ($chainIDsFile) {
    $input{opt} = "-chain_file";
    $input{arg} = $chainIDsFile;
}
else {
    Usage();
}

# Get abs file path of .how file and run dir before switching to run directory
# (and optional dirs/files if specified)
processInputFile(\$input{arg});
processInputFile(\$class_lab_file) if $class_lab_file;
processInputDir(\$run_dir);
processInputDir(\$patches_dir) if $patches_dir;

# Use threshold of 50%, as in Anja's work, to determine S, I and U labelling.
# See prepare_csv_for_weka for more details.
$input{threshold} = 50;

# Set unsupervised to be sent to prepare_csv_for_weka
$input{unsupervised} = $unsupervised ? "-u" : "";

# Vars for running processes
my $process;
my $log_file;
my $cmd;

# Vars for capturing script outputs (STDOUT and STDERR) and exit status
my $combined;
my $success;

# Unless new stat files have been requested, skip those processes responsible
# for creating new stat files, rather than ussing Anja's propensity scores by
# inputting her intf/surf.prepare.out files into calc_patch_propensity
my @skip_processes = $newStatFiles ? () : qw(prepare_intf prepare_surf);

# Skip processes responsible for creating a dir of patches, if a dir of patches
# has been specified
if ($patches_dir) {
    my @more_skips = qw(create_patch_files intf_sorting_patches);

    print "The following processes will be skipped because a patch directory "
        ."has been specified: @more_skips\n";
    
    push(@skip_processes, @more_skips);
}
else {
    # Assign patches dir if not specified
    $patches_dir = "$run_dir/patches_$patch_radius" if ! $patches_dir;
    mkdir $patches_dir if ! -e $patches_dir;
    processInputDir(\$patches_dir);
}

my %dr = get_dirs_hash($patches_dir);

my %statFiles = ();
if ($newStatFiles) {
    # These files will be written to during process 
    $statFiles{intf} = "$run_dir/intf.prepare.out";
    $statFiles{surf} = "$run_dir/surf.prepare.out";
}
elsif ($intfStatFile && $surfStatFile) {
    # If pre-prepared stat files have been supplied, use these
    $statFiles{intf} = File::Spec->rel2abs($intfStatFile);
    $statFiles{surf} = File::Spec->rel2abs($surfStatFile);
}
else {
    
    # These files are the output from when these processes were run on Anja's
    # dataset. In most cases we will want to use some of these files because
    # they supply the numbers that some properties are based on e.g. propensity
    # scores for patches should be based on stats from Anja's dataset, not our
    # own.
    $statFiles{intf}
        = File::Spec->rel2abs("$lib/prop_prep_files/intf.prepare.out");
    
    $statFiles{surf}
        = File::Spec->rel2abs("$lib/prop_prep_files/surf.prepare.out");
}

my %hc = (weka_lib => File::Spec->rel2abs("$lib/weka.jar"));

# If a training set is desired, -n opt is sent to create_arff.pl later on.
# Otherwise, if a dataset has been specified to standardize against use that,
# else use the default of Anja's original training set.
$hc{training} = $setForTraining ? "-n"
    : $standardizeAgainst ? "-t " . File::Spec->rel2abs($standardizeAgainst)
    : "-t " . File::Spec->rel2abs("$lib/sets2StandardizeWith/noABpdbs.filtered.NominalToBinary.train.notismissingATT11.r14_intf50.ahpsSH.sS.arff");

# Output file for create_csv step, used as input for proceeding create_arff
# step.
my $csv_file = "dataset.csv";
my $master_log = 'master.log';

my @processes
    = get_processes_array(\%dr, \%hc, \%statFiles,
                          $csv_file, \%input, $class_lab_file);

# get full path names before changing directory
for my $i (0 .. @processes - 1) {
    my $exec = $lib . $processes[$i]->[1];
    $processes[$i]->[1] =  File::Spec->rel2abs($exec);   
}

# Switch to run directory
$CWD = $run_dir;

## MAIN PROCESS ###############################################################
###############################################################################

my($MASTER, %complete_processes) = process_master_log($master_log);

for my $i ( 0 .. @processes - 1 ) {

    my($process, $exec, $input_string) = @{ $processes[$i] };

    # Check if current process has been completed, or is set to be skipped
    if ( $complete_processes{$process} ) {
        print "Process '$process' logged as complete in master.log, "
            . "skipping\n";
        next;
    }
    elsif ( grep {/^$process$/} @skip_processes ) {
        print "Process '$process' has been set to be skipped, skipping\n";
        next;
    }
    
    print "Running process $process ... ";

    my $log_file = $process . '.log';
    
    open(my $LOG, '>', $log_file)
        or die "Cannot open log file '$log_file', $!";

    my $cmd = "$exec $input_string";

    print "$cmd\n";
    
    ($combined, $success) = qxy($cmd);

    print {$LOG} $combined;
    close $LOG;
    
    if ( ! $success ) {
        die "\nERROR: $process failed. See log '$log_file'\n";
    }
    else {
        print "done\n";
        print {$MASTER} "$process\n";
    }
}

close $MASTER;

print "$0 finished\n";
exit;


## SUBROUTINES ################################################################
###############################################################################

sub Usage {
    print <<EOF;
$0 USAGE : [-c FILE -h FILE ] -u [-p DIR] -i FILE -s FILE -n [-t -T FILE] patchRadius OutputDir

 -p : Specify a directory containing patch files, where line format for each
      line is like this example (i.e. patch summaries):
       <patch A.151> A:146  A:147  A:148  A:149  A:150  A:151
      The files must be named like this example:
       1afvA.patches

 -c : Specify a file containing class labels for patches. The line format is:
       pdbCode(lowercase):chainID(Uppercase):CentralResidueResSeq:label(I or S)

 -i : Specify an interface residue propensity stats file.

 -s : Specify a surface residue propensity stats file.

 -h : Specify a .how file to be used as input.

 -f : Specify a file containing white-space separated pdb chain ids to be used
      as input

 -n : Specify that new residue propensity stat files be made from input data,
      rather than using stat files derived from Anja's original test set, or
      user specified stat files.

 -r : Create training set. This standardizes the final data set against itself.

 -e : Create test set. This standardizes the final set against the .arff file
      specified.

 -u : Unsupervised. If set, then all patches will be missing class labels

If option -p is specified, then option -c must be specfified

If neither -r nor -e are set, then the dataset will be standardized against
Anja's original training set.

This script is a wrapper that runs all required scripts to prepare a dataset
for WEKA input. The script can either be supplied with a .how file (e.g. a
DiscoTope dataset), or directly supplied with a set of patches: see options
for details.

EOF
    exit;
}

sub process_cmdline {

    my $class_label_file = "";
    my $patch_file_dir = "";
    my $how_file = "";
    my $chainIDsFile = "";
    my $newStatFiles;
    my $setForTraining;
    my $standardizeAgainst;
    my $intfStatFile = "";
    my $surfStatFile = "";
    my $unsupervised;
    
    GetOptions("c=s" => \$class_label_file,
               "p=s" => \$patch_file_dir,
               "i=s" => \$intfStatFile,
               "s=s" => \$surfStatFile,
               "h=s" => \$how_file,
               "f=s" => \$chainIDsFile,
               "n"   => \$newStatFiles,
               "r"   => \$setForTraining,
               "e=s" => \$standardizeAgainst,
               "u"   => \$unsupervised);

    # Check option selection is valid
    if ($class_label_file | $patch_file_dir) {
        Usage() unless
            $class_label_file & $patch_file_dir;
    }

    if ($chainIDsFile) {
        croak "Currently can't run a chainIDs file!";
    }
    
    my $patch_radius = shift @ARGV or Usage();
    $patch_radius == int $patch_radius
        or die "Patch radius must be a positive integer";
    
    my $run_dir  = shift @ARGV or Usage();

    return($patch_radius, $how_file, $chainIDsFile,
           $run_dir, $class_label_file, $patch_file_dir,
           $intfStatFile, $surfStatFile,
           $newStatFiles, $setForTraining, $standardizeAgainst,
           $unsupervised);
}

# Given a ref to a file name, this sub checks the file exists and is readable,
# then returns the absolue file path
sub processInputFile {
    my $ref2FileName = shift;

    my $fileName = ${$ref2FileName};
    
    -e -r ($fileName) or croak "$fileName does not exist";
    ${$ref2FileName} = File::Spec->rel2abs($fileName);
}

# Given a ref to a dir name, this sub checks the dir exists,
# then returns the absolue file path
sub processInputDir {
    my $ref2DirName = shift;

    my $dirName = ${$ref2DirName};
    
    -d $dirName or croak "$dirName does not exist";
    ${$ref2DirName} = File::Spec->rel2abs($dirName);
}


# Returns hash of dirs in form:  file type => directory path
# This hash is used in defining cmd-line args for processes to be run
sub get_dirs_hash {
    my $patches_dir = shift
        or die "get_dirs_hash must be passed a patches directory name";

    
    my @dirs = qw(pdb xmas patches_dir log_files absASA rASA
                  intf_sorting_patches intf_data_dir propensities
                  hydrophobicity planarity second_str SSbonds Hbonds
                  psiblast_alignments psiblast_scorecons fosta_alignments
                  fosta_scorecons);
    
    my %dirs_hash = ();
    
    # Create any non-existing directories
    foreach my $file_type (@dirs) {

        my $dir = "";
        
        # Name patches dir after set radius
        if ($file_type eq 'patches_dir') {
            $dir = $patches_dir;
        }
        else {
            $dir = "$run_dir/$file_type";     
        }
       
        if (! -d $dir) {
            mkdir($dir) or die "Cannot make $dir, $!\n";
        }
        $dirs_hash{$file_type} = $dir;
    }

    return %dirs_hash;
}

# Return an array of processes to be run
# Must be passed ref to a hash that is returned by get_dirs_hash
# and a ref to a hash that defines hardcoded file names (see main script)
# Array elements: 0 = process name, 1 = exec path, 2 = input string
sub get_processes_array {
    my $dirs_hashr = shift
        or die "get_processes_array must be passed a ref to a dirs hash";
    my %dr = %{$dirs_hashr};
    
    my $hardcode_hashr = shift
        or die "get_processes_array must be passed ref to  a 'hardcoded' hash";
    my %hr = %{$hardcode_hashr};

    my $statFileHref = shift
        or die "get_processes array must be a ref to a 'statsFiles' hash";
    my %statFiles = %{$statFileHref};
    
    my $csv_file = shift
        or die "get_processes_array must be passed an out csv filename";

    my $inputHref = shift
        or die "get_processes_array must be passed ref to an input hash";
    my %input = %{$inputHref};
    
    my $class_label_file = shift;
    my $classOptString = "";
    if ($class_label_file) {
        $classOptString = "-c $class_label_file";
    }
    
    my @processes
        = ( ['create_process_files', 'create_process_files.pl',
             "$input{opt} $input{arg} -pdb_dir $dr{pdb} -xmas_dir $dr{xmas}"],
            [ 'neisr', 'new.extract.intf.surf.residues_modified.pl',
              "prots.PQS.mimic $dr{xmas}", ],
            [ 'create_patch_files', 'create_patch_files.pl',
              "$patch_radius $dr{pdb} $dr{xmas} patch.centres", ],
            [ 'intf_sorting_patches', 'intf_sorting_patches_modified.pl',
              "-xmas_dir $dr{xmas} -patch_dir $dr{patches_dir} -intf 'intf.residues' -log_dir $dr{log_files} -sort_dir $dr{intf_sorting_patches} -absasa_dir $dr{absASA} -rasa_dir $dr{rASA}", ],
            [ 'prepare_intf','prepare_for_propensities_modified.pl',
              "-in intf.residues -xmas_dir $dr{xmas} -out $statFiles{intf}"],
            [ 'prepare_surf', 'prepare_for_propensities_modified.pl',
              "-in patch.residues -xmas_dir $dr{xmas} -out $statFiles{surf}"],
            [ 'calc_patch_propensity', 'calc_patch_propensity_modified.pl',
              "-patch_dir $dr{patches_dir} -xmas_dir $dr{xmas} -intf_summary $statFiles{intf} -surf_summary $statFiles{surf} -log_dir $dr{log_files} -out_dir $dr{propensities} -v_flag"],
            [  'calc_patch_hydrophobicity', 'calc_patch_hydrophobicity_modified.pl',
               "-patch_dir $dr{patches_dir} -xmas_dir $dr{xmas} -log_dir $dr{log_files} -out_dir $dr{hydrophobicity}" ],
            [ 'calc_patch_planarity', 'calc_patch_planarity_modified.pl',
              "-patch_dir $dr{patches_dir} -pdb_dir $dr{pdb} -log_dir $dr{log_files} -out_dir $dr{planarity}" ],
            [ 'calc_SS_Hb_2Dstr', 'calc_SS_Hb_2Dstr_modified.pl',
              "-patch_dir $dr{patches_dir} -xmas_dir $dr{xmas} -SS $dr{SSbonds} -Hb $dr{Hbonds} -secstr $dr{second_str}", ],
            [ 'psiblast_alignments', 'psiblast_alignments_per_PQSch_modified.pl',
              "-in_dir $dr{patches_dir} -pdb_dir $dr{pdb} -out_dir $dr{psiblast_alignments} -log $dr{log_files}/psiblast_alignments.log" ],
            [ 'psiblast_patch_scorecons', 'calc_patch_blast_scorecons_new.pl',
              "-pdb_dir $dr{pdb} -patch_dir $dr{patches_dir} -aln_dir $dr{psiblast_alignments} -log_dir $dr{log_files} -out_dir $dr{psiblast_scorecons}"],
            [ 'fosta_alignments', 'fosta_alignments_modified.pl',
              "-in_dir $dr{patches_dir} -out_dir $dr{fosta_alignments} -out fosta_alignments.out -log $dr{log_files}" . "/fosta_alignments.log"  ],
            [ 'fosta_scorecons', 'calc_patch_scorecons_modified.pl',
              "-patch_dir $dr{patches_dir} -aln_dir $dr{fosta_alignments} -log_dir $dr{log_files} -out_dir $dr{fosta_scorecons}" ],
            [ 'create_csv', 'prepare_csv_for_weka_modified.pl',
              "-patch_dir $dr{patches_dir} $classOptString -prop_dir $dr{propensities} -hydr_dir $dr{hydrophobicity} -plan_dir $dr{planarity} -sstr_dir $dr{second_str} -rASA_dir $dr{rASA} -aASA_dir $dr{absASA} -SS_dir $dr{SSbonds} -Hb_dir $dr{Hbonds} -class_dir $dr{intf_sorting_patches} -fsc_dir $dr{fosta_scorecons} -bsc_dir $dr{psiblast_scorecons} -single ahpsSH -msa F -intf $input{threshold} -v_flag $input{unsupervised} -out $csv_file",
          ],
            [ 'create_arff', 'create_arff.pl',
              "-c $csv_file -w $hc{weka_lib} $hc{training} > dataset.arff"],
        );
    
    return @processes;
}

# Returns a fh to master log and a hash of complete processes
# (parsed from master log)
sub process_master_log {
    my $master_log = shift
        or die "process_master_log must be passed a master log filename";
    
    if ( ! -e $master_log ) {
        `touch $master_log`;
    }

    open(my $MASTER, '<', $master_log)
        or die "Canot open log '$master_log', $!";

    my %complete_processes = map { chomped $_ => 1 } <$MASTER>;
    
    close $MASTER;
    
    open($MASTER, '>>', $master_log)
        or die "Cannot open log '$master_log' in append mode, $!\n";

    return($MASTER, %complete_processes);
}

sub createpdbIDsAgFile {
    my $patches_dir = shift;

    my $fName = 'pdbIDs.antigen';
    
    open(my $OUT, ">", $fName) or die "Cannot open file $fName, $!";
    opendir(my $DH, $patches_dir) or die "Cannot open dir $patches_dir, $!"; 
    my @pdbIDs;

    while (my $dirFile = readdir($DH)) {
        next if $dirFile =~ /^\./;

        $dirFile =~ s/\.patches//;
        my $pdbID = $dirFile;
        push(@pdbIDs, $pdbID);
    }

    print {$OUT} "$_\n" foreach @pdbIDs;

    return $fName;   
}

__END__
