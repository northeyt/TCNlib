#!/acrm/usr/local/bin/perl -w

# October 2011: 
# 1) added relative and absolute per-patch ASA values 
# and patch size (number of residues in patch yielding ASA values)
# to intf_sorting_patches outputs
# 2) modified this script accordingly to calculate average absASA and rASA

use strict;
use DBI;
use Getopt::Long;  
use lib '/home/bsm/anya/perllib'; #sets dir where Perl_paths.pm is
use Perl_paths; 
use Carp;

####   Define command line inputs   ####
########################################
my $patch_dir;
my $single = '';
my $msa = '';
my $intf_threshold;
my $help = '';
my $v_flag = '';

my $prop_dir;
my $hydr_dir;
my $plan_dir;
my $sstr_dir;
my $rASA_dir;
my $aASA_dir;
my $SS_dir;
my $Hb_dir;
my $class_dir;
my $out_fname;

my $fsc_dir;
my $bsc_dir;

my $class_labels_file;
my $unsupervised;
my $keepUnlabelled;

####   Checking command line options   ####
###########################################
GetOptions('patch_dir=s' => \$patch_dir,
           'prop_dir=s'  => \$prop_dir,
           'hydr_dir=s'  => \$hydr_dir,
           'plan_dir=s'  => \$plan_dir,
           'sstr_dir=s'   => \$sstr_dir,
           'rASA_dir=s'  => \$rASA_dir,
           'aASA_dir=s'  => \$aASA_dir,
           'SS_dir=s'    => \$SS_dir,
           'Hb_dir=s'    => \$Hb_dir,
           'class_dir=s' => \$class_dir,
           'fsc_dir=s'   => \$fsc_dir,
           'bsc_dir=s'   => \$bsc_dir,
           'single=s'    => \$single,
           'msa=s'       => \$msa,
           'intf=i'      => \$intf_threshold,
           'help'        => \$help,
           'v_flag'      => \$v_flag,
           'out=s'       => \$out_fname,
           'c=s'         => \$class_labels_file,
           'u'           => \$unsupervised,
           'k'           => \$keepUnlabelled,
       );
#i.e. /acrm/usr/local/bin/perl prepare_csv_for_weka.pl -patch_dir /acrm/home/anya/interface_prediction/patches/patches_11/ [-single [a][h][p][s][r][w][S][H][F]] [-msa [A][B][s][S][F]] -intf 50 -v_flag

my $OK = &command_line($patch_dir, $help, $single, $msa, $intf_threshold);

if(!$OK)
{
    &print_usage();
    exit;
}

if ($single =~ /F/)
{
    $single = 'ahpsrwSH';
}

if ($msa =~ /F/)
{
    # Removed A and B options (FOSTA and BLAST-based propensities - don't
    # think Anja used these properties in the final WEKA input, it at all
    $msa = 'sS';
}


####   variables and handles   ####
###################################
my %single_dir;
$single_dir{"propensity"}     = $prop_dir;
$single_dir{"hydrophobicity"} = $hydr_dir;
$single_dir{"planarity"}      = $plan_dir;
$single_dir{"secondary_str"}  = $sstr_dir;
$single_dir{"rASA"}           = $rASA_dir;
$single_dir{"absASA"}         = $aASA_dir;
$single_dir{"SSbonds"}        = $SS_dir; #Andrew's SSbonds
$single_dir{"Hbonds"}         = $Hb_dir; 
$single_dir{"class"}          = $class_dir; 
  
my %single_ext;
$single_ext{"propensity"}     = ".propensity.patches";
$single_ext{"hydrophobicity"} = ".hydr.patches";
$single_ext{"planarity"}      = ".patch.rms";
$single_ext{"secondary_str"}  = "..2D";
$single_ext{"rASA"}           = ".rASA";
$single_ext{"absASA"}         = ".absASA";
$single_ext{"SSbonds"}        = ".SSbonds";
$single_ext{"Hbonds"}         = "..Hbbonds";
#class is always processed
$single_ext{"class"}          = ".sorted.patches";

my %msa_dir;
#$msa_dir{"fosta_propensity"} = $Perl_paths::intf_data."";
#$msa_dir{"blast_propensity"} = $Perl_paths::intf_data."";
$msa_dir{"fosta_scorecons"}   = $fsc_dir;
$msa_dir{"blast_scorecons"}   = $bsc_dir;
#$msa_dir{"blast_impact"}       = $Perl_paths::intf_data."";
#$msa_dir{"fosta_impact"}       = $Perl_paths::intf_data."";

my %msa_ext;
$msa_ext{"fosta_scorecons"}   = ".patch.scorecons";
$msa_ext{"blast_scorecons"}   = ".patch.scorecons";
#quite a few missing...

my $cat_num = length($single) + length($msa) + 1; #+1 is for class 

#### Main ####
##############

# If out filename was defined by user, use this. Otherwise, create out fname
my $outfile
    = $out_fname ? $out_fname : &outfile_name($single, $msa, $intf_threshold);

open(OUT, ">$outfile") || die "prepare_csv_for_weka.pl: $outfile not opened\n";

my $header = &header($single, $msa);

#print CSV header
print OUT "$header\n";
print "$header\n" if($v_flag);

#extracting attribute values for patches
my @p_filelist = `ls $patch_dir`;
#my @p_filelist = ("119lA.patches", "1b62A.patches");

# If class_label_file was defined by user, generate hashref for use by get_class
my $classLabelsHref;
if ($class_labels_file) {
    $classLabelsHref = getClassLabelsFromFile($class_labels_file);        
}


foreach my $patchfile (@p_filelist)
{
    chomp $patchfile;
    my $current_column = 0;       

    my %patch_properties = &patch_ids($patch_dir, $patchfile, $cat_num);
    printf "\n%d patches in $patchfile.\n", scalar keys %patch_properties if ($v_flag);

    #single-seq-based data
    if ($single)
    {
        if ($single =~ /a/)
        {
            my $valuecount = &get_value_per_patch(\%patch_properties, $patchfile, $single_dir{"propensity"}, $single_ext{"propensity"}, $current_column);
            printf "Found single AA propensity value for %d patches.\n", $valuecount if ($v_flag);
            $current_column++;
        }
        
        if ($single =~ /h/)
        {
            my $valuecount = &get_value_per_patch(\%patch_properties, $patchfile, $single_dir{"hydrophobicity"}, $single_ext{"hydrophobicity"}, $current_column);            
            printf "Found single AA hydrophobicity value for %d patches.\n", $valuecount if ($v_flag);
            $current_column++;
        }
        
        if ($single =~ /p/)
        {            
            my $valuecount = &get_value_per_patch(\%patch_properties, $patchfile, $single_dir{"planarity"}, $single_ext{"planarity"}, $current_column);
            
            printf "Found single AA planarity value for %d patches.\n", $valuecount if ($v_flag);
            $current_column++;
        } 

        if ($single =~ /s/)
        {            
            my $valuecount = &get_value_per_patch(\%patch_properties, $patchfile, $single_dir{"secondary_str"}, $single_ext{"secondary_str"}, $current_column);
            
            printf "Found single 2D value for %d patches.\n", $valuecount if ($v_flag);
            $current_column++;
        }

        if ($single =~ /r/)
        {
            my $valuecount = &get_value_per_patch(\%patch_properties, $patchfile, $single_dir{"rASA"}, $single_ext{"rASA"}, $current_column);

            printf "Found single rASA value for %d patches.\n", $valuecount if ($v_flag);
            $current_column++;
        }

        if ($single =~ /w/)
        {
            my $valuecount = &get_value_per_patch(\%patch_properties, $patchfile, $single_dir{"absASA"}, $single_ext{"absASA"}, $current_column);

            printf "Found single absASA value for %d patches.\n", $valuecount if ($v_flag);
            $current_column++;
        }

        if ($single =~ /S/)
        {            
            my $valuecount = &get_value_per_patch(\%patch_properties, $patchfile, $single_dir{"SSbonds"}, $single_ext{"SSbonds"}, $current_column);
            
            printf "Found single SSbonds value for %d patches.\n", $valuecount if ($v_flag);
            $current_column++;
        }

        if ($single =~ /H/)
        {            
            my $valuecount = &get_value_per_patch(\%patch_properties, $patchfile, $single_dir{"Hbonds"}, $single_ext{"Hbonds"}, $current_column);
            
            printf "Found single H-bonds value for %d patches.\n", $valuecount if ($v_flag);
            $current_column++;
        }
    }
         
    #MSA-based data                        
    if($msa)
    {
        if ($msa =~ /S/)
        {
            my $valuecount;
            my $ret
                = eval {
                    $valuecount = &get_value_per_patch(\%patch_properties, $patchfile, $msa_dir{"fosta_scorecons"}, $msa_ext{"fosta_scorecons"}, $current_column);
                    1;
                };

            if (! $ret ) {
                print "Found no FOSTA scorecons file\n"
            }
            else {
                printf "Found FOSTA scorecons value for %d patches.\n",
                    $valuecount if ($v_flag);
            }
            $current_column++;
        } 
        
        if ($msa =~ /s/)
        {
            my $valuecount;
            my $ret
                = eval {
                    $valuecount = &get_value_per_patch(\%patch_properties, $patchfile, $msa_dir{"blast_scorecons"}, $msa_ext{"blast_scorecons"}, $current_column);
                    1;
                };
            
            if (! $ret) {
                print "Found no BLAST scorecons file\n";
            }
            else {
            printf "Found BLAST scorecons value for %d patches.\n",
                $valuecount if ($v_flag);
            }
            $current_column++;
        }                          
    }

    #adding class values (S/U/I) to the last column
    my $classcount
        = &get_class(\%patch_properties, $patchfile, $single_dir{"class"},
                     $single_ext{"class"}, $current_column, $intf_threshold,
                     $classLabelsHref);
    
    printf "Found CLASS value for %d patches.\n", $classcount if ($v_flag);

    #use Data::Dumper;
    #print Dumper \%patch_properties;
    
    #printing CSV instances to outfile
    foreach my $k (sort keys %patch_properties)
    {
        my $value = $patch_properties{$k};

        #added 25.8.2010.
        my $class = chop($value);

        if ($unsupervised || ($class eq "U" && $keepUnlabelled)) {
            # If unsupervised is set, or instance is unlabelled and
            # keepUnabelledthen has been set, then set class label to missing
            print OUT "$k,$value" . "?\n";
        }
        elsif ($class eq "I" || $class eq "S")
        {
            # print I and S instances to file, as well as U if set.
            print OUT "$k,$value$class\n";
        }
    }
}
close(OUT);
exit;

#### Subroutines ####
#####################
sub get_value_per_patch
{
    my ($r_hash, $patchfile, $dir, $ext, $colnum) = @_;

    my $value_count = 0;
    $patchfile =~ /(.*)\.patches/;
    my $pqs_id = $1;

    my @paths
        = ( "$dir/$pqs_id" . "$ext",
            "$dir/" . substr($pqs_id, 0, 4) . '.' . substr($pqs_id, 4, 1)
            . "$ext"
        );
    
    my $path
        = -e $paths[0] ? $paths[0] :
          -e $paths[1] ? $paths[1] :
              croak "Could not find a valid path. Attempted paths: @paths";
    
   
    my $chain = chop($pqs_id);
    
    if(-e $path)
    {
        open(V, "$path");
        while (my $line = <V>)
        {
            chomp $line;
            if (length($line))
            {
                $line =~ /patch\s(\S+)>\s(.*)$/;
                my $central_res = $1;            
                my $curr_value = $2;
                
                $central_res =~ s/\./:/g;

                if($curr_value =~ /N\/A/)
                {
                    $curr_value = '?';
                }

                my $key = $pqs_id.":".$central_res;
                my $csv_value = $$r_hash{$key};

                #changes the appropriate column value
                my $new_value = &save_value($csv_value, $curr_value, $key, $colnum, $ext);
                $$r_hash{$key} = $new_value;                
                $value_count++;
            }
        }
        close(V);
    }
    else {
        print "path $path does not exist!\n";
    }
    
    return $value_count;
}


##############################################################
sub get_class
{
    my ($r_hash, $patchfile, $dir, $ext,
        $colnum, $threshold, $classLabelsHref) = @_;
    
    my $value_count = 0;
    $patchfile =~ /(.*)\.patches/;
    my $pqs_id = $1;
    
    my $path = "$dir/$pqs_id$ext";
    
    my $chain = chop($pqs_id);
    
    if(-e $path){
        open(V, "$path");
        while (my $line = <V>){
            chomp $line;
            
            if (length($line)){
                $line =~ /patch\s(\S+)>\s(.*)$/;
                my $central_res = $1;            
                my $curr_value = $2;
                my $class = '?';
                
                $central_res =~ s/\./:/g;

                if($curr_value =~ /N\/A/){
                    $class = '?';
                }
                #difference from get_value_per_patch:
                else{
                    my ($central, $perc_patch, $perc_intf, $abs_patcharea,
                        $rel_patcharea, $rescount_in_patch) = split(' ', $curr_value);
                    
                    if ($perc_patch == 0){
                        $class = 'S';
                    }
                    elsif ($perc_patch <= $threshold){
                        $class = 'U';
                    }
                    else{
                        $class = 'I';
                    }
                }
                
                my $key = $pqs_id.":".$central_res;
                my $csv_value = $$r_hash{$key};

                #print "TEST: old $csv_value\n";
                
                #changes the appropriate column value
                my $new_value = &save_value($csv_value, $class, $key, $colnum, $ext);
                $$r_hash{$key} = $new_value;
                #print "TEST: new $new_value\n";
                $value_count++;
            }
        }
        close(V);
    }
    elsif ($classLabelsHref) {

        my $firstKey = join(":", ($pqs_id,$chain));
        
        foreach my $patchKey (keys %{$classLabelsHref->{$firstKey}}) {
            my $csv_value = $$r_hash{$patchKey};

            my $class = $classLabelsHref->{$firstKey}->{$patchKey};            
            
            my $new_value = &save_value($csv_value, $class, $patchKey, $colnum, $ext);
            $$r_hash{$patchKey} = $new_value;
            ++$value_count;
        }
    }
    
    return $value_count;
}

##############################################################
sub outfile_name
{
    my ($single, $msa, $threshold) = @_;
    my $outfile = "intf$threshold.";

    #creating outfile name
    if(length($single))
    {
        $outfile = $outfile.$single;
    }

    $outfile = $outfile.".";
    
    if(length($msa))
    {
        $outfile = $outfile.$msa;
    }

    $outfile = $outfile.".csv";
    $outfile =~ s/\.\./\./g;
    print "Outfile: $outfile\n";

    return $outfile;
}


##############################################################
sub header
{
    my ($single, $msa) = @_;

    my $header = 'patchID,';
    
    if ($single)
    {
        if ($single =~ /a/)
        {
            $header = $header."propensity,";
        }
        if ($single =~ /h/)
        {
            $header = $header."hydrophobicity,";
        }
        if ($single =~ /p/)
        {
            $header = $header."planarity,";
        }
        if ($single =~ /s/)
        {
            $header = $header."secondary_str,";
        }
        if ($single =~ /r/)
        {
            $header = $header."rASA,";
        }
        if ($single =~ /w/)
        {
            $header = $header."absASA,";
        }
        if ($single =~ /S/)
        {
            $header = $header."SSbonds,";
        }
        if ($single =~ /H/)
        {
            $header = $header."Hbonds,";
        }
    }
    if($msa)
    {
        if ($msa =~ /S/)
        {
            $header = $header."fosta_scorecons,";
        }
        if ($msa =~ /s/)
        {
            $header = $header."blast_scorecons,";
        }
    }
    $header = $header."intf_class";

    return $header;
}


##############################################################
sub save_value
{
    my ($old, $col_value, $key, $column, $ext) = @_;
    my $new = '';

    if($old)
    {
        my @old_in_columns = split(',', $old);
        if ($old_in_columns[$column] ne '?')
        {
            print "already have value for $key???\n";
        }
        else
        {
            $old_in_columns[$column] = $col_value;             
        }

        foreach my $c (@old_in_columns)
        {
            $new = $new.$c.",";
        }
        chop($new); #removing last comma
    }
    else
    {
        print "for $ext found patch $key, but does not exist in patch_dir?\n";
    }

    #print "$col_value:$column [$old][$new]\n";

    return $new;
}


##############################################################
sub patch_ids
{
    my ($dir, $file, $categories) = @_;
    my %patches = ();

    #has '?' for every category
    my $cat_string = &create_cat_string($categories);

    $file =~ /(.*)\.patches/;
    my $pqs_id = $1; #id:chain, actually
    my $chain = chop($pqs_id);

    open(P, "$dir/$file") || die "prepare_csv_for_weka.pl: (P handle) $dir$file not opened\n";

    while (my $line = <P>)
    {
        chomp $line;

        if (length($line))
        {
            $line =~ /patch\s(\S+)>\s(.*)$/;
            my $central_res = $1;
            $central_res =~ s/\./:/g;
            
            my $key = $pqs_id.":".$central_res;
            
            $patches{$key} = $cat_string;
        }
    }
    close(P);

    return %patches;
}

sub getClassLabelsFromFile {
    my $file = shift;

    open(my $fh, "<", $file) or die "Cannot open file $file, $!";

    my %class = ();
    
    while (my $line = <$fh>) {
        chomp $line;
        my ($pdbCode, $chainID, $centralResSeq, $label) = split(":", $line);
        my $firstKey = join(":", ($pdbCode, $chainID));
        my $secndKey = join(":", ($firstKey, $centralResSeq));
        $class{$firstKey} = {} if ! exists $class{$firstKey};
        
        $class{$firstKey}->{$secndKey} = $label;
    }
    return \%class;
}


##############################################################
#will return comma-separated list of '?', $num of them
sub create_cat_string
{
    my $num = $_[0];
    my $string = '';

    for(my $i=0; $i<$num; $i++)
    {
        $string = $string."?,";
    }
    chop($string); # removing the last ,
    
    return $string;
}

##############################################################
sub print_usage
{
    print <<EOF;
USAGE:
$0 -patch_dir DIR -intf INT [-single STRING] [-msa STRING] [-v_flag] [-c FILE]
i.e. $0 -patch_dir /acrm/home/anya/interface_prediction/patches/patches_11/ [-single [a][h][p][s][r][w][S][H][F]] [-msa [A][B][s][S][F]] [-c FILE]

This code will produce a CSV file containing header with the list of attributes.
Last attribute is the class.
Every subsequent row is an instance of attribute values.
    -patch_dir is directory where patch files are
    -intf is threshold (in %) above which patches are classified as interface
     patches I, below intf are unclassified U, if overlap is 0 then surface S
    -single extracts single-sequence-based attributes:
        a=amino acid propensities
        h=hydrophobicity
        p=planarity (to add protrusion)
        s=secondary structure
        r=relative solvent-accessibility
        w=absolute solvent-accessibility
        S=SS-bonds
        H=H-bonds
    -msa extracts mustiple-sequence-alignment-based attributes:
        A=amino acid propensities (FOSTA-based)
        B=amino acid propensities (BLAST-based)
        S=sequence conservation (FOSTA-based)
        s=sequence conservation (BLAST-based)(to add ImPACT)

    -c Specify a file to overwrite patch class labels.
       Line must follow this format:
       pdbCode(lowercase):chainID(Uppercase):CentralResidueResSeq:label(I or S0)

Needs at least one attribute to be selected to run.
Running -single F and -msa F will run full version, including all options.
Outfile will be intf[threshold].[options].csv\n\n";
EOF
}


##############################################################
sub command_line
{
    my ($patch_dir, $help, $single, $msa, $intf) = @_;
    my $categories_count = length($single) + length($msa);
    
    if( ($help) || (!$patch_dir) )#|| (!$intf) )
    {
        return 0;
    }
    elsif (!$categories_count)
    {
        print "ERROR!!! Select at least one attribute!\n";
        return 0;
    }
    else
    {
        return 1;
    }
}
