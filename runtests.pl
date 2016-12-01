#!/usr/bin/env perl
use TAP::Harness;
use Cwd;

# Set these so that tests can find modules and external packages
$ENV{'TCNlib'} = getcwd();
$ENV{'PERL5LIB'} = $ENV{'TCNlib'} . "/lib:" . $ENV{'PERL5LIB'};

my $harness = TAP::Harness->new(\%args);
my @tests 
    = qw(ARFF.t  cdhit.t  FOSTA.t  MSA.t  scorecons.t  VectorCalcs.t  WEKA.t  write2tmp.t
pdb/automatic_patches.t        pdb/hbondFinder.t  pdb/pdbFunctions.t        pdb/solv.t
pdb/BLAST.t                    pdb/pdbsws.t       pdb/ss.t
pdb/chaincontacts.t            pdb/pdb.t          pdb/ViewPatch.t
pdb/file_cache.t               pdb/makepatch.t    pdb/RadiusFinder.t
pdb/fixGetInterfaceResidues.t  pdb/multiChain.t   pdb/remote_file_getter.t
pdb/get_files.t                pdb/patch_desc.t   pdb/secstrCalculator.t);

@tests = map {'t/' . $_} @tests;

$harness->runtests(@tests);
