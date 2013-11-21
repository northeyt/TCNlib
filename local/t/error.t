#!/usr/bin/perl -Iblib/lib -Iblib/arch -I../blib/lib -I../blib/arch
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl error.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl error.t'
#   Tom Northey <zcbtfo4@acrm18>     2013/11/21 12:54:29

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use Test::More qw( no_plan );
BEGIN { use_ok( 'local::error' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my @args = ( 'message', "Some message", 'type', 'an_error' );
my $error = new_ok( 'local::error' , \@args );

is($error->id(), 1, "error id assignment works okay");
