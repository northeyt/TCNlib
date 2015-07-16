#!/acrm/usr/local/bin/perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl princip.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl princip.t'
#   Tom Northey <zcbtfo4@acrm18>     2015/04/15 09:28:31

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More qw( no_plan );
BEGIN { use_ok( 'pdb::princip' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my $testObj = new_ok("pdb::princip");
$testObj->input("4houA_patch112.pdb");
ok($testObj->run, "run ok");

is($testObj->getPlanarity, 7.0838456, "getPlanarity works ok");

