#!/usr/bin/env perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl RadiusFinder.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl RadiusFinder.t'
#   Tom Northey <zcbtfo4@acrm18>     2015/08/26 16:34:51

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More qw( no_plan );
use pdb;
BEGIN { use_ok( 'pdb::RadiusFinder' ); }
use FindBin qw($RealBin);
chdir($RealBin); # So test data files can be found

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

subtest 'constructor' => sub {
    new_ok('pdb::RadiusFinder');
};

subtest 'findRadiusOfAtomFromNames' => sub {
    my $testObj  = pdb::RadiusFinder->new(radiiFile => 'radii.dat');
    my $resName  = 'ALA';
    my $atomName = 'N';
    is($testObj->findRadiusOfAtomFromNames($resName, $atomName), '1.65',
       "expected radius returned for ALA, N");

    my $resName  = 'ARG';
    my $atomName = 'CA';
    is($testObj->findRadiusOfAtomFromNames($resName, $atomName), '1.87',
       "expected radius returned for ARG, CA");    
};
