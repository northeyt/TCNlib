#!/acrm/usr/local/bin/perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl pdbsws.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl pdbsws.t'
#   Tom Northey <zcbtfo4@acrm18>     2013/11/11 10:32:01

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More qw( no_plan );
BEGIN { use_ok( 'pdb::pdbsws' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my $pdbsws = pdb::pdbsws->new();

is( [ $pdbsws->get_ac('4houA') ]->[0], 'P09914', "get_ac returns ac okay" );

is( ref $pdbsws->get_ac('4houZ'), 'local::error', 
    "get_ac returns error when no ac found" );
