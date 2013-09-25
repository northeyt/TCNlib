#!/usr/bin/perl -Iblib/lib -Iblib/arch -I../blib/lib -I../blib/arch
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl pcentres.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl pcentres.t'
#   Tom Northey <zcbtfo4@acrm18>     2013/09/12 16:25:35

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use Test::More qw( no_plan );
BEGIN { use_ok( pcentres ); }


#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.


