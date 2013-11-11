#!/acrm/usr/local/bin/perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl pdb.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl pdb.t'
#   Tom Northey <zcbtfo4@acrm18>     2013/09/11 15:07:05

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use Data::Dumper;

use lib ( '..' );
use Test::More qw( no_plan );
use Test::Exception;

BEGIN { use_ok('pdb'); }

use xmas2pdb;

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.


# Atom object tests
my $ATOM_line
    = "ATOM     13  CG  LEU A 148      13.227  45.000  23.178  1.00 47.53           C";

my $atom = atom->new( ATOM_line => $ATOM_line );

isa_ok($atom, 'atom', "atom object created ok");

# pdb object tests

my $pdb_file = '1djs.pdb';

my $pdb = pdb->new( pdb_file => $pdb_file, pdb_code => '1djs' );
isa_ok($pdb, 'pdb', "pdb object created ok");

# fh builder test
is( ref $pdb->pdb_data, 'ARRAY', "fh builder okay" );

## _parse_ATOM_lines

my @ATOM_lines = $pdb->_parse_ATOM_lines;

## _parse_atoms

$pdb->_parse_atoms();

$pdb->atom_array();

# atom_index

$pdb->atom_index();

ok($pdb->resid_index, "resid_index ok" );

# chain object

my $chain_id = 'A';
my $chain
    = chain->new( pdb_file => $pdb_file,
                  pdb_code => '1djs',
                  chain_id => 'A' );

isa_ok($chain, 'chain', "chain object created ok");

## are around modifiers for _parse_ATOM_lines working?

$chain->_parse_atoms();

# get accession codes using pdbsws?

is( $chain->accession_codes()->[0], 'P21802',
    "Accession codes from pdbsws" );


# reading ASAs and radii from xmas2pdb object

my $xmas2pdb_file = 'xmas2pdb';
my $xmas_file     = '1djs.xmas';
my $radii_file    = 'radii.dat';
my $form          = 'multimer';

my %arg = (xmas2pdb_file => $xmas2pdb_file,
           xmas_file     => $xmas_file,
           radii_file    => $radii_file,
           form          => $form,
       );

my $x2p = xmas2pdb->new(%arg);

my $test_atom = $pdb->atom_array->[1];

$pdb->read_ASA($x2p);

ok($test_atom->radius(), "Radius read from xmas2pdb object" );
ok($test_atom->ASAc(), "Multimer ASA read from xmas2pdb object" );

$arg{form} = 'monomer';

my $mono_x2p = xmas2pdb->new(%arg);

$pdb->read_ASA($mono_x2p);

ok($test_atom->ASAm(), "Monomer ASA read from xmas2pdb object" );

# patch_centres

my %pc_arg = ( ASA_threshold => 25 );

my ($errors, $patch_centres) = $pdb->patch_centres( %pc_arg );

#foreach my $index (1360, 1357, 1361, 1364, 1359, 1363, 1362, 1358) {
#    print "$index " .  $pdb->atom_array->[$index]->{serial} . ' '
#         . $pdb->atom_array->[$index]->ASAc() . "\n";
#}

($errors, $patch_centres) = $chain->patch_centres( %pc_arg);

ok( @{ $errors }, "Errors returned if ASAm has not been set for chain" );

$chain->read_ASA($mono_x2p);

ok( $chain->patch_centres( %pc_arg),
    "patch_centres modified for chain object" );
