#!/usr/bin/env perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl xmas2pdb.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl xmas2pdb.t'
#   Tom Northey <zcbtfo4@acrm18>     2013/09/13 14:57:41

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use Data::Dumper;

use lib ( '..' );

use Test::More qw( no_plan );
BEGIN { use_ok( 'pdb::xmas2pdb' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my $xmas2pdb_file = 'xmas2pdb';
my $xmas_file     = '1djs.xmas';
my $radii_file    = 'radii.dat';
my $form          = 'multimeric';

    
my %arg = (
    radii_file    => $radii_file,
    xmas2pdb_file => $xmas2pdb_file,
    xmas_file     => $xmas_file,
    form          => $form,
);

my @arg = %arg;

my $xmas2pdb = new_ok( 'xmas2pdb' => \@arg );

can_ok( $xmas2pdb, ('output') );
can_ok( $xmas2pdb, ('output_file') );

# test chain selection

%arg = (
    radii_file    => $radii_file,
    xmas2pdb_file => $xmas2pdb_file,
    xmas_file     => $xmas_file,
    form          => $form,
    chain_ids => ['A'],
);

$xmas2pdb = xmas2pdb->new(%arg);
