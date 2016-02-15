#!/usr/bin/perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl ss.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl ss.t'
#   Tom Northey <zcbtfo4@acrm18>     2015/07/22 13:36:36

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;
use Test::More qw( no_plan );
use Test::Deep;
BEGIN { use_ok( 'pdb::ss' ); }
use FindBin qw($RealBin);
chdir($RealBin); # So test data files can be found

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.


new_ok('pdb::ssFinder');
new_ok('pdb::ss', [atomSerialPairAref => [1, 2], resIDPairAref => ["A1", "B1"],
              distance => "1.2"]);

subtest "parsing sslist output lines" => sub {
    my $testLine = "  A23  Atom   164 :   A88  Atom   682 : 2.076";
    my $ss       = pdb::ssFinder->_getssFromLine($testLine);

    cmp_deeply($ss->atomSerialPairAref, [164, 682],    "atom serials parsed ok");
    cmp_deeply($ss->resIDPairAref,      [qw(A23 A88)], "resIDs parsed ok");
    is($ss->distance, 2.076, "distance parsed ok");
};

subtest "running getssArray" => sub {
    my $testFile = "1qok.pdb";
    my @ssArray   = pdb::ssFinder->new(input => $testFile)->getssArray();

    cmp_deeply(\@ssArray, array_each(isa("pdb::ss")), "getssArray returns an array of ss objects")
        or explain \@ssArray;

    is(scalar @ssArray, 2, "and array is expected length");
}
