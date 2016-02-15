#!/usr/bin/env perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl pdb2xmas.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl pdb2xmas.t'
#   Tom Northey <zcbtfo4@acrm18>     2014/02/14 16:57:19

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use Test::More qw( no_plan );
use Data::Dumper;

use lib ('..');

BEGIN { use_ok( 'pdb::pdb2xmas' ); }
use FindBin qw($RealBin);
chdir($RealBin); # So test data files can be found

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my $pdb_file = '1djs.pdb';

my $test;

ok($test = pdb::pdb2xmas->new(pdb_file => $pdb_file), "new pdb2xmas ok");

ok($test->output(), "pdb2xmas->output() works ok");
ok($test->last_output(), "pdb2xmas->last_output() works ok");
ok($test->log(), "pdb2xmas->log() works ok" );
