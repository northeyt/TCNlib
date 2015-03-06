#!/acrm/usr/local/bin/perl
# This script tests the speed of automatic_patches under a number of different
# conditions, primarily to assess the efficacy of forking during repeated
# make_patch calls
use strict;
use warnings;
use pdb::pdb;
use pdb::automatic_patches;
use Benchmark;

# 1afv is a large PDB file containing 9000 atom lines
my $testPdbCode = "1afv";
my $testPdbFile = "$testPdbCode.pdb";
my $testXmasFile = "$testPdbCode.xmas";

my $testPdb = pdb->new(pdb_code => $testPdbCode,
                       pdb_file => $testPdbFile,
                       xmas_file => $testXmasFile,
                   );

# Read ASAc values for pdb
$testPdb->read_ASA();

# Iterate through 10 pdb objects derived from the test pdb object, each with a
# proportion of atoms from the original, from 1:10 to 10:10.
# This will allow us to vary the size of the pdb object
for (my $i = 1 ; $i < 11 ; ++$i) {

    my $partition = int(scalar @{$testPdb->atom_array()} / 100 * $i* 10) - 1;
    
    my @paritionAtoms = @{$testPdb->atom_array()}[0..$partition];

    # Remove any atoms from residues that no longer have C-alpha
    # after partitioning
    my $cAlphaFlag = 0;
    until ($cAlphaFlag) {
        if ($paritionAtoms[-1]->name eq 'CA') {
            $cAlphaFlag = 1;
        }
        else {            
            pop @paritionAtoms;
        }
    }

    print "Initializing parition object with " . @paritionAtoms
        . " atoms\n";
    
    my $partitionObj = pdb->new(pdb_code => $testPdbCode,
                                pdb_file => $testPdbFile,
                                atom_array => \@paritionAtoms,
                                has_read_ASA => 1,
                                resid2RelASAHref => $testPdb->resid2RelASAHref);
    
    foreach my $flag (0, 1) {
        my $aP = automatic_patches->new(pdb_object => $partitionObj,
                                        patch_type => 'normal',
                                        radius => 8,
                                        ASA_type => 'ASAc',
                                        forkFlag => $flag,
                                        pdb_code => $testPdbCode,
                                );
        
        
        my $t0 = Benchmark->new();
        
        $aP->get_patches();
        my $t1 = Benchmark->new();
        my $td = timediff($t1, $t0);
        print "Flag $flag, TIME: ", timestr($td), "\n";
    }
}
