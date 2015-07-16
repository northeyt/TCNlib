#!/usr/bin/perl -Iblib/lib -Iblib/arch -I../blib/lib -I../blib/arch
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl cluster.t'

# Test file created outside of h2xs framework.
# Run this like so: `perl cluster.t'
#   THOMAS NORTHEY <tcn@THOMASs-MacBook-Pro-2.local>     2014/02/06 11:40:46

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;
use Test::More qw( no_plan );
use Test::Deep;

use Data::Dumper;

use lib ('..');
BEGIN { use_ok( 'cluster' ); }

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my @array = ( 1, 5, -4, 9, 124, 15, 300 );

my $cluster = cluster->new(array => \@array, condition => \&test_condition,);

my @expected = ( [1,5,9], [-4], [124], [15], [300],  );

my @clusters = $cluster->array_of_clusters();

cmp_deeply( \@clusters, \@expected, "returned expected clusters" );

sub test_condition {
    my($cluster, $j) = @_;

    my @cluster = @{ $cluster };

    foreach my $i (@cluster) {
        if ( abs ( $i - $j )  < 5 ) {
            return 1;
        }
    }
    return 0;
}
