#!/usr/bin/perl -Iblib/lib -Iblib/arch -I../blib/lib -I../blib/arch
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl fixGetInterfaceResidues.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl fixGetInterfaceResidues.t'
#   Tom Northey <zcbtfo4@acrm18>     2015/02/24 15:07:43

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;
use Test::More qw( no_plan );
use pdb::pdb;

my $errChain = chain->new(pdb_code => '2iff', chain_id => 'Y',
                          pdb_file => '2iff.pdb');

# This pdb was previously proving tricky because solvent residues were being
# processed. These are now correctly skipped.
ok($errChain->getInterfaceResidues([$errChain->createOtherChains()]), "getInterfaceResidues works for problem case");
