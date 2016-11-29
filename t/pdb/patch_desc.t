#!/usr/bin/env perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl patch_desc.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl patch_desc.t'
#   Tom Northey <zcbtfo4@acrm18>     2013/10/10 10:23:54

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use lib ( '..' );
use pdb;

use Data::Dumper;
use Math::Trig;
use Math::VectorReal;
use Test::More qw( no_plan );
use Test::Deep;

BEGIN { use_ok( 'pdb::patch_desc' ); }
use FindBin qw($RealBin);
chdir($RealBin); # So test data files can be found

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.


# Test _true_angle_from_x
print "Testing _true_angle_from_x\n";

my @v = ( vector(1,1,0), vector(-1,1,0), vector(-1,-1,0), vector(1,-1,0) );

is(pdb::patch_desc->_true_angle_from_x($v[0]) *  180 / pi, 45,
   "1st quad angles okay");

is(pdb::patch_desc->_true_angle_from_x($v[1]) *  180 / pi, 135,
   "2nd quad angles okay");

is(pdb::patch_desc->_true_angle_from_x($v[2]) * 180 / pi, 225,
   "3rd quad angles okay");

is(pdb::patch_desc->_true_angle_from_x($v[3]) * 180 / pi, 315,
   "4th quad angles okay");


my $chain = chain->new(pdb_code => '4hou',
                       chain_id => 'A',
                       pdb_file => 'pdb4hou.ent',
                       xmas_file => 'pdb4hou.xmas',
                   );

my $dump = '4houA_patch112.dump';

open(my $fh, '<', $dump) or die "Cannot open patch object dump $dump";

my $data;

{
    local $/;
    $data = <$fh>;
    close $fh;
}

my $patch = eval $data;

my $patch_desc = pdb::patch_desc->new( patch => $patch,
                                  parent => $chain );

isa_ok( $patch_desc, 'pdb::patch_desc' );

# Hashes created internally by patch order method
my ( $surface_error, $surface_atoms ) = $patch_desc->_surface_atoms;

my %all_atom = map { $_->serial => $_ } @{ $patch_desc->parent->atom_array };

my %surf_atom = map { $_->serial => $_ } @{ $surface_atoms };

my %patch_atom = map { $_->serial => $all_atom{$_->serial} }
    @{ $patch_desc->patch->atom_array };

my %patch_surf_atom
    = map { $_->serial => $_ }
    grep ( defined $patch_atom{$_->serial}, @{ $surface_atoms } );


my %nonpatch_atom
    = map { $all_atom{$_}->serial => $all_atom{$_} }
    grep( ! defined $patch_atom{$_}, keys %all_atom );


is( [$patch_desc->patch_order]->[0], 'ALRETQ', "patch order works okay" )
    or diag explain $patch_desc->patch_order;

my $file = 'fullrun.pdb';

open($fh, '>', $file) or die "Cannot open $file, $!";

print {$fh} sort { $a->serial <=> $b->serial } values %all_atom; 
