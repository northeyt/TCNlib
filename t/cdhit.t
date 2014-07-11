#!/acrm/usr/local/bin/perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl cdhit.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl cdhit.t'
#   Tom Northey <zcbtfo4@acrm18>     2014/06/30 14:54:15

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;
use Test::More qw( no_plan );
BEGIN { use_ok( 'cdhit' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my $testObj = new_ok(cdhit => [input => 'test.fasta']);
ok($testObj->getClusters(), "getClusters() works okay");
