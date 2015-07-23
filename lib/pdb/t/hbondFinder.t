#!/usr/bin/env perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl hbondFinder.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl hbondFinder.t'
#   Tom Northey <zcbtfo4@acrm18>     2015/07/22 16:57:13

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;
use Test::More qw( no_plan );
use Test::Deep;
BEGIN { use_ok( 'pdb::hbondFinder' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.


new_ok('pdb::hbondFinder');
new_ok('pdb::Hb', [donorSerial => 1, acceptorSerial => 2]);

subtest "getHbonds" => sub {
    my $testFile = "1qok.pdb";
    my @gotArray
        = pdb::hbondFinder->new(input => $testFile)->getHbonds();

    # TOMORROW: Test all elements of array are pdb::Hb objects
    #           Test correct number of elements returned
    #           Test a Hb object for acceptor and donor serials
    cmp_deeply(\@gotArray, array_each(isa("pdb::Hb")),
               "getHbonds returns array of Hb");

    is(scalar @gotArray, "137", "that has correct number of Hbs");
    
    my ($testHb) = grep {$_->donorSerial() == 17} @gotArray;
    is($testHb->acceptorSerial(), 176,
       "donor and acceptor correctly parsed in one test case");
};
