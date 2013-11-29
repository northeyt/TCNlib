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
use Test::Deep;
use Test::Exception;

BEGIN { use_ok('pdb'); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.


# Atom object tests
my $ATOM_line
    = "ATOM     13  CG  LEU A 148      13.227  45.000  23.178  1.00 47.53           C";

my $atom = atom->new( ATOM_line => $ATOM_line );

isa_ok($atom, 'atom', "atom object created ok");

# HETATM

my $HETATM_line
    = "HETATM 2105  N   MSE B  67      68.539  35.469  16.161  1.60  0.00           N  \n";

my $het_atom = atom->new( ATOM_line => $HETATM_line );

is( $het_atom->is_het_atom, 1, "het atom identified" );

my $het_atom_string = "$het_atom";

is( $het_atom_string, $HETATM_line, "stringify okay for het atoms" );

# pdb object tests

my $pdb_file = '1djs.pdb';

my $pdb = pdb->new( pdb_file => $pdb_file, pdb_code => '1djs' );
isa_ok($pdb, 'pdb', "pdb object created ok");

# BUILD method gets attributes from getresol objet

cmp_deeply( [ $pdb->experimental_method(), $pdb->resolution(),
              $pdb->r_value ],
            [ 'crystal', '2.40', '0.227' ],
            "BUILD method gets attributes from getresol object okay" );


# fh builder test
is( ref $pdb->pdb_data, 'ARRAY', "fh builder okay" );

## _parse_ATOM_lines

my @ATOM_lines = $pdb->_parse_ATOM_lines;

ok( (grep /^HETATM/, @ATOM_lines),
    "HETATM lines included by _parse_ATOM_lines");

## _parse_atoms

$pdb->_parse_atoms();

is( scalar @{ $pdb->atom_array() }, 3066,
    "all ATOM and HETATM lines parsed and included in atom_array" );

# _parse_ter

my $ter_line = "TER    2725      SO4 A 407 \n";

my($serial, $chainID) = pdb::_parse_ter($ter_line);

is( $serial, '2725', "_parse_ter parses serial ok"  );
is( $chainID, 'A'  , "_parse_ter parses chainID ok" );

# pdb file containing multiple models - for the moment, exception needs
# to be thrown (untill I write in capacity to deal with multi-model pdbs)

my $multimodel_file = '/acrm/data/pdb/pdb2q1z.ent';

dies_ok(
    sub { pdb->new( pdb_file => $multimodel_file, pdb_code => '2q1z' )
      },  "Croak  if pdb is multimodel" );

# atom_ and resid_ index

ok($pdb->atom_index(), "atom_index ok");
ok($pdb->resid_index, "resid_index ok" );

# test out chain solvent determination

my $multi_term_pdb_code = '2we8';
my $multi_term_file = '/acrm/data/pdb/pdb2we8.ent';

my $mterm_pdb = pdb->new( pdb_code => $multi_term_pdb_code,
                          pdb_file => $multi_term_file, );

$mterm_pdb->atom_array();

# is_terminal atom attribute

my $term_atom  = $pdb->atom_array->[ $pdb->resid_index->{A362}->{CB}  ];
my $term_atom2 = $pdb->atom_array->[ $pdb->resid_index->{B140}->{OD2} ];

cmp_deeply( [ $term_atom->is_terminal, $term_atom2->is_terminal ], [ 1, 1 ],
            "terminal atom labelled");

my $term_string = $term_atom->stringify_ter();
my $exp_string  = "TER    1628  CB  ALA A 362 \n\n";
is($term_string, $exp_string, "stringify_ter returns TER string ok" );

# terminal_atom_index

print "$_\n" foreach  @{ $pdb->terminal_atom_index };

# chain object

my $chain_id = 'A';
my $chain
    = chain->new( pdb_file => $pdb_file,
                  pdb_code => '1djs',
                  chain_id => 'A' );

isa_ok($chain, 'chain', "chain object created ok");

## are around modifiers for _parse_ATOM_lines working?

$chain->_parse_atoms();

# Chain length working?

is($chain->chain_length, 206, "chain length determined ok");


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

($errors, $patch_centres) = $chain->patch_centres( %pc_arg);

ok( @{ $errors }, "Errors returned if ASAm has not been set for chain" );

$chain->read_ASA($mono_x2p);

ok( $chain->patch_centres( %pc_arg),
    "patch_centres modified for chain object" );

# test multi_resName_resid

my $bad_chain = chain->new( pdb_code => '3u5e',
                            chain_id => 'O',
                            pdb_file => '3u5e.pdb'
                        );

$bad_chain->atom_array();

is( scalar keys %{ $bad_chain->multi_resName_resid() }, '175',
    "multi_resName_resid captures resids with multi resName-atoms" );
