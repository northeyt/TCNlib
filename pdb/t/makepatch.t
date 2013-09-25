#!/usr/bin/perl -Iblib/lib -Iblib/arch -I../blib/lib -I../blib/arch
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl makepatch.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl makepatch.t'
#   Tom Northey <zcbtfo4@acrm18>     2013/09/12 13:38:30

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use warnings;
use strict;
use Data::Dumper;

use lib ( '..' );
use pdb;

use Test::More qw( no_plan );
BEGIN { use_ok( 'makepatch' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my $makepatch_file = `which makepatch`;
chomp $makepatch_file;

my $pdb_file       = '1djs.pdb';
my $radius         = 8;
my $patch_type     = 'contact';

my $atom_line = "ATOM     31  CD  PRO A 150      16.450  43.163  16.346  1.00 44.51           C";

my $atom = atom->new( ATOM_line => $atom_line );

my $makepatch = makepatch->new( makepatch_file => $makepatch_file,
                                pdb_file       => $pdb_file,
                                patch_type     => $patch_type,
                                radius         => $radius,
                                central_atom   => $atom, );


my $output = $makepatch->output;

print Dumper $output;

my $patch = new_ok( 'patch' => [ central_atom => $makepatch->central_atom,
                                 pdb_file => $makepatch->output ] );

print "patch->new() okay when passed a makepatch object directly?\n";

my $dir_patch = new_ok( 'patch' => [ $makepatch ] );

print Dumper $dir_patch->atom_array();

print $dir_patch->summary();

$dir_patch->run_PatchOrder;
