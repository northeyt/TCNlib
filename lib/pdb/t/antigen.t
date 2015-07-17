#!/usr/bin/env perl -Iblib/lib -Iblib/arch -I../blib/lib -I../blib/arch
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl antigen.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl antigen.t'
#   Tom Northey <zcbtfo4@acrm18>     2013/09/23 16:32:31

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use Data::Dumper;

use lib ( '..' );
use Test::More qw( no_plan );
use pdb::pdb;

#BEGIN { use_ok( 'antigen' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.


my $antigen = pdb->new( pdb_code => '1djs',
                        pdb_file => '1djs.pdb',
                        antigen_chain_ids => 'A',
                        epitope_residue_array => [ qw( A147 A148 A151 ) ],
               );

print Dumper $antigen->epitope_atom_index;

# Problems with 'resSeq not recognised'

$antigen = chain->new( pdb_code => '2p8l',
                       pdb_file => '/acrm/data/pdb/pdb2p8l.ent',
                       chain_id => 'C',
                       'epitope_residue_array'
                           => [ qw(C0 C2 C3 C4 C5 C6 C7 C8 C10) ],
                );

print Dumper $antigen->atom_array;

$antigen->epitope_atom_index;
