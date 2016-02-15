#!/usr/bin/env perl
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl pdbFunctions.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl pdbFunctions.t'
#   Tom Northey <zcbtfo4@acrm18>     2014/06/24 17:07:15

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;
use Test::More qw( no_plan );
use Test::Deep;

BEGIN { use_ok( 'pdb::pdbFunctions' ); }
use FindBin qw($RealBin);
chdir($RealBin); # So test data files can be found

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.


# test compareResSeqs
my @test = qw(1 1A);
my %order = ("1A" => "1", "1" => "2");
my @exp  = qw(1A 1);
cmp_deeply([sort {pdb::pdbFunctions::compare_resSeqs($a, $b, \%order)} @test],
           \@exp, "compareResSeqs sorts 1, 1A => 1A, 1 (using order hash)");

@test = qw(1A 1B);
%order = ("1A" => "2", "1B" => "1");
@exp  = qw(1B 1A);
cmp_deeply([sort {pdb::pdbFunctions::compare_resSeqs($a, $b, \%order)} @test],
           \@exp, "compareResSeqs sorts 1A, 1B => 1B, 1A (using order hash)");

@test = qw(52A 52);
%order = ("52" => "1", "52A" => "2");
@exp  = qw(52 52A);
cmp_deeply([sort {pdb::pdbFunctions::compare_resSeqs($a, $b)} @test],
           \@exp, "compareResSeqs sorts 52A, 52 => 52, 52A (using order hash)");

@test = qw(52 52A);
cmp_deeply([sort {pdb::pdbFunctions::compare_resSeqs($a, $b)} @test],
           \@exp, "compareResSeqs keeps 52, 52A the same");

@test = qw(52B 52A);
@exp  = qw(52A 52B);

cmp_deeply([sort {pdb::pdbFunctions::compare_resSeqs($a, $b)} @test],
           \@exp, "compareResSeqs sorts 52B, 52A => 52A, 52B");
 
my $riA = "C121";
my $riB = "C27";

is(pdb::pdbFunctions::compare_resids($riA, $riB), 1,
   "compare_resids works okay");
 
