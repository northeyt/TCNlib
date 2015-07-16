#!/usr/bin/perl -Iblib/lib -Iblib/arch -I../blib/lib -I../blib/arch
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl parseXMAS.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl parseXMAS.t'
#   Tom Northey <zcbtfo4@acrm18>     2015/04/15 13:33:58

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More qw( no_plan );
use Test::Deep;
BEGIN { use_ok( 'pdb::parseXMAS' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

$testObj = new_ok("pdb::parseXMAS", [xmasData => getXmasData()]);

$testObj->parseXMAS;

cmp_deeply($testObj->ppHbDonor2AcceptorHref, expppHbDonor2AcceptorHref(),
       "Hydrogen bonds parsed ok");

cmp_deeply($testObj->SSbondHref, expSSbondHref(),
           "SS bonds parsed ok");

testResID2secStructHref($testObj);

sub testResID2secStructHref {
    my $testObj = shift;

    is(ref $testObj->resID2secStructHref, 'HASH',
       "testResID2secStructHref returns Href");
    
    my @resIDs        = keys   %{$testObj->resID2secStructHref};
    my @sStructLabels = values %{$testObj->resID2secStructHref};

    my $resIDRe        = re('\w+\.\w+'); # e.g. B.123
    my $sStructLabelRe = re('[IiHhBbCEeTtSGg]'); # Valid sect struct flags
    
    cmp_deeply(\@resIDs, array_each($resIDRe),
               "resID hash keys ok");
    cmp_deeply(\@sStructLabels, array_each($sStructLabelRe),
               "sStruct labels ok");
}

sub getXmasData {
    my $testFile = "1djs.xmas";
    open(my $IN, "<", $testFile) or die "Cannot open file $testFile, $!";

    my @data = <$IN>;

    return \@data;
}

sub expppHbDonor2AcceptorHref {
    return {735 => 39, 41 => 289, 57 => 281, 83 => 269, 114 => 97,
            123 => 101, 2497 => 113, 149 => 134, 148 => 1724};
}

sub expSSbondHref {
    return {265 => 683, 683 => 265, 1041 => 1482, 1482 => 1041};
}
