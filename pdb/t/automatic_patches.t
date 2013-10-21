#!/usr/bin/perl -Iblib/lib -Iblib/arch -I../blib/lib -I../blib/arch
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl automatc_patches.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl automatc_patches.t'
#   Tom Northey <zcbtfo4@acrm18>     2013/09/16 17:35:35

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use Data::Dumper;

use lib ( '..' );

use Test::More qw( no_plan );
BEGIN { use_ok( 'automatic_patches' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.


my %arg
    = ( radius => 8,
        patch_type => 'normal',
        pdb_code => '1djs',
    );

my $auto = automatic_patches->new(%arg);
$auto->get_patches();

# Test to see if tmp xmas file write works
$arg{pdb_code} = '1nox';
my $noxmas = automatic_patches->new(%arg);

# Directly access pdb_file and pdb code  attr for testing purposes
$noxmas->{pdb_code} = '1noxm'; # Will not be found in xmas dir
$noxmas->{pdb_file} = '1nox.pdb';

ok($noxmas->xmas_file, "xmas file created when not found in xmas dir");

